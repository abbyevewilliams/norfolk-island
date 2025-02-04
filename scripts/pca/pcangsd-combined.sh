#!/bin/bash
#SBATCH --time=2-00:00:00 
#SBATCH --job-name=PCAngsd
#SBATCH --partition=long
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/PCAngsd_%A.log
#SBATCH --error=logs/PCAngsd_%A.error
#SBATCH --mail-type=END
#SBATCH --mail-user=abigail.williams@biology.ox.ac.uk

ml Anaconda3

# Path to environment
source activate /data/biol-silvereye/ball6625/norfolk-island/scripts/pca/pca-env

# Paths
OUT_DIR=/data/biol-silvereye/ball6625/norfolk-island/outputs/pca
mkdir -p $OUT_DIR

# Concatenate museum and contemporary beagle files
zcat $OUT_DIR/museum.beagle.gz $$OUT_DIR/contemporary.beagle.gz | gzip > $OUT_DIR/museum-contemporary.beagle.gz

# Run pcangsd
pcangsd -b $OUT_DIR/museum-contemporary.beagle.gz \
    --iter 1000 \
    -o $OUT_DIR/museum-contemporary

