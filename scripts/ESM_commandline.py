import concurrent.futures
import subprocess
import json
import tqdm

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
            
# batch fasta sequences into groups of 5 and write each batch into a .fasta file
from itertools import zip_longest

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

max_batches = 1
from math import ceil
batch_size = ceil(len(fasta_sequences) / max_batches)
batch_files = []

for i, (batch, batch_names) in enumerate(zip(grouper(fasta_sequences, batch_size, ''), grouper(sequence_names, batch_size, ''))):
    batch_file_name = f'batch_{i+1}.fasta'
    batch_files.append(batch_file_name)
    with open(f'batches/{batch_file_name}', 'w') as batch_file:
        for name, sequence in zip(batch_names, batch):
            if name and sequence:  # Avoid writing empty sequences or names
                batch_file.write(f'{name}\n{sequence}\n')

def generate_pdb(batch_file):
    """Function to run a Docker command for a given batch file of FASTA sequences."""
    command = f'docker run --rm --gpus all -v $(pwd)/batches:/input -v $(pwd)/outputs:/output rexpository/esmfold -i /input/{batch_file} -o /output/{batch_file}.pdb > ./logs/{batch_file}_pred.log 2>./logs/{batch_file}_pred.err'
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.stdout.decode().strip()
    except subprocess.CalledProcessError as e:
        return e.stderr.decode().strip()

# Using ThreadPoolExecutor to run the commands in parallel on batch files
with concurrent.futures.ThreadPoolExecutor() as executor:
    # Map the function to the batch files and execute them
    futures = {executor.submit(generate_pdb, batch_file): batch_file for batch_file in batch_files}
    for future in tqdm.tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Processing batch files"):
        result = future.result()  # Ensure the container has finished executing
        print(f"Batch file {futures[future]} finished procesing")