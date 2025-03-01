
import sys
import os
import pandas as pd
import numpy as np
from Bio import SeqIO  # Biopython for handling DNA sequences

# Retrieve command-line arguments
args = sys.argv
chrom = int(args[1]) - 1  # Convert 1-based index to 0-based
fasta_file = args[2]  # Path to the input FASTA file

# Ensure the FASTA file exists
if not os.path.isfile(fasta_file):
    print("Error: FASTA file not found.")
    sys.exit(1)

# Read the genome sequence from the FASTA file
records = list(SeqIO.parse(fasta_file, "fasta"))
seq = records[chrom].seq  # Extract chromosome sequence
seq_length = len(seq)

# Create output directories if they don't exist
microhomology_dir = "microhomology_output"
summary_dir = "summary_microhomology"
os.makedirs(microhomology_dir, exist_ok=True)
os.makedirs(summary_dir, exist_ok=True)

# Detect microhomology regions for maxLength = 50, 100, 150, 200
for maxLength in range(50, 201, 50):
    for kmer in range(6, 11):
        output = []  # Store detected positions

        for i in range(0, seq_length - maxLength):
            short = seq[i:(i + kmer)]  # Extract k-mer
            if "N" not in short:  # Ignore sequences with 'N'
                shortRevC = short.reverse_complement()  # Reverse complement
                fragment = seq[i:(i + maxLength - 1)]  # Extract fragment
                start_pos = min(kmer + 11, len(fragment))  # Ensure valid search start
                pos = fragment.find(shortRevC, start=start_pos)  # Search reverse complement

                if pos > 0:  # If found
                    output.append([i, pos + kmer])  # Store start and end positions

        # Save detected microhomology regions
        microhomology_file = f"{microhomology_dir}/kmer{kmer:02d}_chrom{chrom+1:02d}_maxLength{maxLength:03d}.csv"
        pd.DataFrame(output, columns=["sp", "end"]).to_csv(microhomology_file, index=False)

# Summarize microhomology coverage
window_size = 1000  # Define the window size
for maxLength in range(50, 201, 50):  # Loop over maxLength values (50, 100, 150, 200)
    pass_array = np.zeros(seq_length, dtype=int)  # Initialize array

    for kmer in range(6, 11):
        microhomology_file = f"{microhomology_dir}/kmer{kmer:02d}_chrom{chrom+1:02d}_maxLength{maxLength:03d}.csv"

        # Read the detected microhomology data
        try:
            csv_data = pd.read_csv(microhomology_file)
        except FileNotFoundError:
            print(f"Warning: {microhomology_file} not found. Skipping...")
            continue

        # Adjust the end positions
        csv_data["end"] = csv_data["sp"] + csv_data["end"] - 1

        # Mark detected positions
        for _, row in csv_data.iterrows():
            sp = max(row["end"] - window_size, 0)  # Ensure non-negative start
            ep = row["sp"]  # End position
            pass_array[int(sp):int(ep)] = 1  # Mark detected range

        # Save summary results
        summary_file = f"{summary_dir}/kmer{kmer:02d}_chrom{chrom+1:02d}_maxLength{maxLength:03d}.csv"
        summary_df = pd.DataFrame({"pass": [np.count_nonzero(pass_array)], "total": [seq_length]})
        summary_df.to_csv(summary_file, index=False)

