#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --job-name=index_ref
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END
#SBATCH --mail-user=abigail.williams@biology.ox.ac.uk
#SBATCH --output=logs/index_%A.log
#SBATCH --error=logs/index_%A.error

# Bwa-mem2 path
export PATH=/data/biol-silvereye/ball6625/software/bwa-mem2-2.2.1_x64-linux:$PATH

# Paths
REF=/data/biol-silvereye/ref_genome/Zlat_2_Tgut_pseudochromosomes.shortChromNames.fasta.gz

# Index reference (only needed once)
if [ ! -f "$REF.bwt.2bit.64" ]; then
    bwa-mem2 index $REF
fi