import concurrent.futures
import subprocess
import json
import tqdm

sequence_names = []
fasta_sequences = []

with open('input.fasta', 'r') as file:
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

def generate_pdb(sequence):
    """Function to run a cURL command for a given FASTA sequence in a Docker container."""
    command = f'docker run --rm rexpository/curl curl --insecure -X POST --data "{sequence}" https://ap i.esmatlas.com/foldSequence/v1/pdb/'
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.stdout.decode().strip()
    except subprocess.CalledProcessError as e:
        return e.stderr.decode().strip()

# Using ThreadPoolExecutor to run the commands in parallel
with concurrent.futures.ThreadPoolExecutor() as executor:
    # Map the function to the FASTA sequences and execute them
    results = list(tqdm.tqdm(executor.map(generate_pdb, fasta_sequences), total=len(fasta_sequences), desc="Processing sequences"))

    # Process the results and save them in the outputs folder
    for sequence_name, pdb_result in zip(sequence_names, results):
        with open(f"outputs/{sequence_name}.txt", 'w') as file:
            file.write(pdb_result)
