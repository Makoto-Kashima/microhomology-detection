# 🧬 Microhomology Detection and Summary Script

This Python script detects **microhomology regions** in a **FASTA genome sequence** by identifying **reverse complement k-mers (6-10 bp)**.  
It then summarizes the distribution of detected microhomology across different fragment lengths.

---

##  Features
✔ Reads a **FASTA file** and extracts the sequence of a given chromosome  
✔ Searches for **reverse complement k-mers (6-10 bp)** within a given sequence length (`maxLength = 50, 100, 150, 200`)  
✔ **Saves detected microhomology regions** as CSV files  
✔ **Summarizes microhomology coverage** over different fragment lengths  
✔ Supports **multiple `maxLength` values in a single execution**  

---

## Installation
### Clone this repository
```bash
git clone https://github.com/Makoto-Kashima/microhomology-detection
cd microhomology-detection
```

### Install required dependencies
This script requires **Python 3** and the following Python libraries:
```bash
pip install biopython pandas numpy
```

## Usage
To run the script, use:
```bash
python script.py <chromosome_number> <fasta_file>
```

## Example
```bash
python script.py 1 GCF_000002035.5_GRCz10_genomic.fna
```
Arguments:
<chromosome_number> → Chromosome index (1-based)
<fasta_file> → Path to the FASTA genome file
The script will automatically process maxLength values: 50, 100, 150, and 200.
