# # Code Explanation
# # This script processes a FASTA file containing genomic sequences and identifies microhomology regions by searching for reverse complement k-mer matches.
# #Example Usage
# #python script.py 1 500 GCF_000002035.5_GRCz10_genomic.fna
# #This runs the script on chromosome 1, with a maximum fragment length of 500.
# #The script searches for microhomology regions in the given FASTA file.
# #The results are saved in microhomology_output/.
# # The results are stored in a Pandas DataFrame with two columns:
# # "sp" → Start position of the match.
# # "end" → End position of the match.
# # ex)
# # sp,end
# # 5,161
# # 31,169
# # 83,191
# 
# import sys
# import os
# from Bio import SeqIO  # Import SeqIO from Biopython to read FASTA files
# from Bio.Seq import Seq  # Import Seq to handle nucleotide sequences
# import pandas as pd  # Import pandas for handling and saving data
# 
# # Retrieve command-line arguments
# args = sys.argv
# chrom = int(args[1]) - 1  # Convert chromosome number from 1-based to 0-based index
# maxLength = int(args[2])  # Maximum sequence length to analyze
# fasta_file = args[3]  # Path to the input FASTA file
# 
# # Check if the provided FASTA file exists
# if not os.path.isfile(fasta_file):
#     print("Error: FASTA file not found.")
#     sys.exit(1)
# 
# # Read the FASTA file and store sequences in a list
# records = list(SeqIO.parse(fasta_file, "fasta"))
# 
# # Extract the sequence of the specified chromosome
# seq = records[chrom].seq
# 
# # Create an output directory if it doesn't already exist
# if not os.path.exists("microhomology_output"):
#     os.mkdir("microhomology_output")
# 
# # Iterate through k-mer sizes from 6 to 10
# for kmer in range(6, 11):
#     output = []  # List to store detected microhomology positions
# 
#     # Slide through the sequence within the valid range
#     for i in range(0, len(seq) - maxLength):
#         short = seq[i:(i + kmer)]  # Extract a k-mer sized substring
#         contains_N = short.find("N")  # Check if the sequence contains "N" (unknown base)
# 
#         if contains_N == -1:  # Continue only if "N" is not present
#             shortRevC = short.reverse_complement()  # Get the reverse complement of the k-mer
#             fragment = seq[i:(i + maxLength - 1)]  # Extract a fragment of maxLength
#             
#             # Ensure that the starting position for searching remains within bounds
#             start_pos = min(kmer + 11, len(fragment))  
#             pos = fragment.find(shortRevC, start=start_pos)  # Find reverse complement within the fragment
# 
#             if pos > 0:  # If a match is found
#                 match = [i, pos + kmer]  # Store the start and end positions
#                 output.append(match)
# 
#     # Convert the detected matches into a Pandas DataFrame
#     df = pd.DataFrame(output, columns=["sp", "end"])
# 
#     # Save results as a CSV file in the output directory
#     df.to_csv(
#         "microhomology_output/kmer%02d_chrom%02d_maxLength%03d.csv"
#         % (kmer, chrom + 1, maxLength),
#         index=False
#     )
# 

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
            sp = max(row["end"] - window_size + 1, 0)  # Ensure non-negative start
            ep = row["sp"]  # End position
            pass_array[int(sp):int(ep)] = 1  # Mark detected range

            # Save summary results
        summary_file = f"{summary_dir}/kmer{kmer:02d}_chrom{chrom+1:02d}_maxLength{maxLength:03d}.csv"
        summary_df = pd.DataFrame({"pass": [np.count_nonzero(pass_array)], "total": [seq_length]})
        summary_df.to_csv(summary_file, index=False)

