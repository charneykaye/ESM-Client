import concurrent.futures
import subprocess
import json
import tqdm
from itertools import cycle

sequence_names = []
fasta_sequences = []

with open('inputs/input.fasta', 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith('>'):  # Start of a new sequence
            if 'current_sequence_name' in locals():
                joined_sequence = ''.join(current_sequence_parts).replace('*', '').replace('X', 'G')
                if ':' in joined_sequence:
                    parts = joined_sequence.split(':')
                    for j, part in enumerate(parts):
                        fasta_sequences.append(part)
                        complex_name = (current_sequence_name + f"|Complex_{j+1}")[:70]  # Truncate
                        sequence_names.append(complex_name)
                else:
                    fasta_sequences.append(joined_sequence)
                    sequence_names.append(current_sequence_name[:70])  # Truncate
            current_sequence_name = line
            current_sequence_parts = []
        else:  # Non-empty line that is part of the current sequence
            current_sequence_parts.append(line)
    # Handle the last sequence in the file
    if 'current_sequence_name' in locals():
        joined_sequence = ''.join(current_sequence_parts).replace('*', '').replace('X', 'G')
        if ':' in joined_sequence:
            parts = joined_sequence.split(':')
            for j, part in enumerate(parts):
                fasta_sequences.append(part)
                complex_name = (current_sequence_name + f"|Complex_{j+1}")[:70]  # Truncate
                sequence_names.append(complex_name)
        else:
            fasta_sequences.append(joined_sequence)
            sequence_names.append(current_sequence_name[:70])  # Truncate

# # batch fasta sequences into groups and write each batch into a .fasta file
# from itertools import zip_longest
# from itertools import cycle
# def grouper(iterable, n, fillvalue=None):
#     "Collect data into fixed-length chunks or blocks"
#     # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
#     args = [iter(iterable)] * n
#     return zip_longest(*args, fillvalue=fillvalue)

# max_batches = 4
# max_gpus = 4
# from math import ceil
# batch_size = ceil(len(fasta_sequences) / max_batches)
# batch_files = []

# for i, (batch, batch_names) in enumerate(zip(grouper(fasta_sequences, batch_size, ''), grouper(sequence_names, batch_size, ''))):
#     batch_file_name = f'batch_{i+1}.fasta'
#     batch_files.append(batch_file_name)
#     with open(f'batches/{batch_file_name}', 'w') as batch_file:
#         for name, sequence in zip(batch_names, batch):
#             if name and sequence:  # Avoid writing empty sequences or names
#                 batch_file.write(f'{name}\n{sequence}\n')



# batch fasta sequences into groups and write each batch into a .fasta file
# balancing the number of characters in fasta_sequences between all batches
max_batches = 1
max_gpus = 1
batch_files = []

total_chars = sum(len(seq) for seq in fasta_sequences)
avg_chars_per_batch = total_chars / (max_batches )  # Calculate average for one less batch to allow the last batch to go over

# Helper function to write a batch to a file and print stats
def write_batch(batch_num, names, sequences):
    batch_file_name = f'batch_{batch_num}.fasta'
    batch_files.append(batch_file_name)
    num_sequences = len(names)
    total_length = sum(len(seq) for seq in sequences)
    with open(f'batches/{batch_file_name}', 'w') as batch_file:
        for name, sequence in zip(names, sequences):
            batch_file.write(f'>{name}\n{sequence}\n')
    print(f'Batch {batch_num} written with {num_sequences} sequences and {total_length} total characters.')

# Create batches
current_batch_num = 1
current_batch_chars = 0
current_batch_names = []
current_batch_sequences = []
for name, sequence in zip(sequence_names, fasta_sequences):
    current_batch_names.append(name)
    current_batch_sequences.append(sequence)
    current_batch_chars += len(sequence)
    # If the current batch has reached the average number of characters, write it to a file
    # Allow the last batch to exceed the average
    if current_batch_chars >= avg_chars_per_batch and current_batch_num < max_batches:
        write_batch(current_batch_num, current_batch_names, current_batch_sequences)
        current_batch_num += 1
        current_batch_chars = 0
        current_batch_names = []
        current_batch_sequences = []

# Write the last batch regardless of its size
write_batch(current_batch_num, current_batch_names, current_batch_sequences)


def generate_pdb(batch_file, gpu_id):
    """Function to run a Docker command for a given batch file of FASTA sequences on a specific GPU."""
    command = f'sudo docker run --rm --gpus device={gpu_id} -v $(pwd)/batches:/input -v $(pwd)/outputs:/output rexpository/esmfold -i /input/{batch_file} -o /output/{batch_file}.pdb > ./logs/{batch_file}_pred.log 2>./logs/{batch_file}_pred.err'
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.stdout.decode().strip()
    except subprocess.CalledProcessError as e:
        return e.stderr.decode().strip()

# Using ThreadPoolExecutor to run the commands in parallel on batch files across multiple GPUs
with concurrent.futures.ThreadPoolExecutor() as executor:
    # Map the function to the batch files and assign them to GPUs in a round-robin fashion
    gpu_ids = cycle(range(max_gpus)) 
    futures = {executor.submit(generate_pdb, batch_file, next(gpu_ids)): batch_file for batch_file in batch_files}
    for future in tqdm.tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Processing batch files"):
        result = future.result()  # Ensure the container has finished executing
        # Process Complete Notification
        batch_file_name = futures[future]
        total_chars = sum(len(seq) for seq in fasta_sequences if seq)
        # gpu_ids is a cycle object and not directly printable in a meaningful way
        # We need to keep track of the GPU ID used for each batch file
        gpu_id_for_batch = next((gpu_id for gpu_id, file in zip(cycle(range(max_gpus)), batch_files) if file == batch_file_name), None)
        print(f"Batch file {batch_file_name} finished processing on GPU {gpu_id_for_batch} with {total_chars} characters in fasta sequences.")
