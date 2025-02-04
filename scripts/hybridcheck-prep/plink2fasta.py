import argparse

def tped_to_fasta(tped_file, tfam_file, output_fasta):
    # Read sample IDs from .tfam
    with open(tfam_file, 'r') as f:
        samples = [line.split()[1] for line in f]  # Extract individual IDs

    # Initialize sequences for each individual
    sequences = {sample: [] for sample in samples}

    # Read .tped and construct sequences
    with open(tped_file, 'r') as f:
        for line in f:
            cols = line.strip().split()
            alleles = cols[4:]  # Genotypes start from the 5th column
            for i, sample in enumerate(samples):
                genotype = alleles[i * 2:(i * 2) + 2]  # Two alleles per individual
                # Convert to a single character (e.g., A, T, C, G, or N for missing data)
                if genotype[0] == "0" or genotype[1] == "0":
                    sequences[sample].append("N")  # Missing data
                elif genotype[0] == genotype[1]:
                    sequences[sample].append(genotype[0])  # Homozygous
                else:
                    sequences[sample].append("N")  # Heterozygous sites as ambiguous

    # Write to FASTA
    with open(output_fasta, 'w') as f:
        for sample, seq in sequences.items():
            f.write(f">{sample}\n")
            f.write("".join(seq) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert PLINK .tped and .tfam files to a FASTA alignment.")
    parser.add_argument("--tped", required=True, help="Path to the input .tped file.")
    parser.add_argument("--tfam", required=True, help="Path to the input .tfam file.")
    parser.add_argument("--output", required=True, help="Path to the output FASTA file.")

    args = parser.parse_args()

    tped_to_fasta(tped_file=args.tped, tfam_file=args.tfam, output_fasta=args.output)