#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=32G
#SBATCH --job-name=align_stats
#SBATCH --partition=short
#SBATCH --mail-type=END
#SBATCH --output=logs/stats_%A.log
#SBATCH --error=logs/stats_%A.error
#SBATCH --mail-user=abigail.williams@biology.ox.ac.uk

# Load samtools
ml SAMtools/1.16.1-GCC-11.3.0

# Set path to bams
BAM_PATH=/data/biol-silvereye/ball6625/norfolk-island/outputs/bams

# Get sample names
mapfile -t SAMPLES < samples.txt

# Run stats in series
for SAMPLE in "${SAMPLES[@]}"; do

  echo "Calculating stats for $SAMPLE..."

  samtools stats "${BAM_PATH}/${SAMPLE}.sorted.bam" > "${BAM_PATH}/${SAMPLE}.stats"

  echo "Stats for $SAMPLE complete."


done
