import concurrent.futures
import subprocess
import json
import tqdm
from itertools import cycle
from math import exp

sequence_names = []
fasta_sequences = []

# Parse the input FASTA file
with open('inputs/input_small.fasta', 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith('>'):  # Start of a new sequence
            if 'current_sequence_name' in locals():
                joined_sequence = ''.join(current_sequence_parts).replace('*', '').replace('X', 'G') # Remove Stop Codons (*) and replae unknown amino acids (X) with G
                # Split the sequence into multiple sequence complexes if it contains a colon
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


# Batch fasta sequences into groups and write each batch into a .fasta file, balancing the run time of each batch
max_batches = 4
max_gpus = 4
batch_files = []

# Calculate the weight of each sequence using the derived equation
weights = [2.3499 * exp(0.0050 * len(seq)) for seq in fasta_sequences]
total_weight = sum(weights)

# Helper function to write a batch to a file and print stats
def write_batch(batch_num, batch_data):
    batch_file_name = f'batch_{batch_num}.fasta'
    batch_files.append(batch_file_name)
    num_sequences = len(batch_data)
    total_length = sum(len(seq) for _, seq, _ in batch_data)
    total_weight = sum(weight for _, _, weight in batch_data)
    with open(f'batches/{batch_file_name}', 'w') as batch_file:
        for name, sequence, _ in batch_data:
            batch_file.write(f'>{name}\n{sequence}\n')
    print(f'Batch {batch_num} written with {num_sequences} sequences, {total_length} total characters, and {total_weight:.2f} total weight.')

# Greedy algorithm to create batches with balanced weights
# Sort sequences by weight in descending order
sorted_sequences = sorted(zip(sequence_names, fasta_sequences, weights), key=lambda x: x[2], reverse=True)
batches = [[] for _ in range(max_batches)]
batch_weights = [0] * max_batches

for name, sequence, weight in sorted_sequences:
    # Find the batch with the minimum weight and add the current sequence to it
    min_batch_index = batch_weights.index(min(batch_weights))
    batches[min_batch_index].append((name, sequence, weight))
    batch_weights[min_batch_index] += weight

# Write all batches to files
for i, batch_data in enumerate(batches, start=1):
    write_batch(i, batch_data)


# 
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
    # Map the function to the batch files and assign them to GPUs
    gpu_ids = cycle(range(max_gpus)) 
    futures = {executor.submit(generate_pdb, batch_file, next(gpu_ids)): batch_file for batch_file in batch_files}
    for future in tqdm.tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Processing batch files"):
        result = future.result()  # Ensure the container has finished executing
        # Process Complete Alert
        batch_file_name = futures[future]
        # Calculate the total characters for the current batch
        total_chars_in_batch = sum(len(seq) for name, seq, _ in batches[batch_files.index(batch_file_name)])
        # Keep track of the GPU ID used for each batch file
        gpu_id_for_batch = next((gpu_id for gpu_id, file in zip(cycle(range(max_gpus)), batch_files) if file == batch_file_name), None)
        
        print(f"Batch file {batch_file_name} finished processing on GPU {gpu_id_for_batch} with {total_chars_in_batch} characters in fasta sequences.")