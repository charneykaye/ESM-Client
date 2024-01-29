import concurrent.futures
import subprocess
import json
import tqdm

sequence_names = []
fasta_sequences = []

with open('inputs/input.fasta', 'r') as file:
    for line in file:
        line = line.strip()
        if not line:  # Skip empty lines
            continue
        if line.startswith('>'):
            if ':' in line:
                sequence_names.append(line + '|Complex_1')  # Append |Complex_1 to the first sequence name
            else:
                sequence_names.append(line)
        else:
            # Replace unwanted characters and split if it's a protein complex
            parts = line.replace('*', '').replace('X', 'G').split(':')
            for j, part in enumerate(parts):
                fasta_sequences.append(part)
                if len(parts) > 1:
                    # Ensure the complex number increments without keeping the previous complex name
                    complex_name = sequence_names[-1].split('|Complex_')[0] + f"|Complex_{j+2}"
                    sequence_names.append(complex_name)
                else:
                    # Append the sequence name without incrementing the complex number
                    sequence_names.append(sequence_names[-1].split('|Complex_')[0] + '|Complex_1')
            # Remove the last name added in the sequence_names list because it's a duplicate
            sequence_names.pop()
            
# batch fasta sequences into groups of 5 and write each batch into a .fasta file
from itertools import zip_longest

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

batch_size = 5
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
    command = f'docker run --rm --gpus device={} -v $(pwd)/batches:/input -v $(pwd)/outputs:/output rexpository/esmfold -i /input/{batch_file} -o /output/{batch_file}.pdb > ./logs/{batch_file}_pred.log 2>./logs/{batch_file}_pred.err'
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.stdout.decode().strip()
    except subprocess.CalledProcessError as e:
        return e.stderr.decode().strip()

# Using ThreadPoolExecutor to run the commands in parallel on batch files
with concurrent.futures.ThreadPoolExecutor() as executor:
    # Map the function to the batch files and execute them
    results = list(tqdm.tqdm(executor.map(generate_pdb, batch_files), total=len(batch_files), desc="Processing batch files"))

    # # Process the results and save them in the outputs folder
    # for batch_file_name, pdb_result in zip(batch_files, results):
    #     # Adjust the output file name to match the batch file name
    #     with open(f"outputs/{batch_file_name}.txt", 'w') as file:
    #         file.write(pdb_result)