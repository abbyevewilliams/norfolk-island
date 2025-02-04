#!/bin/bash
#SBATCH --time=12:00:00 
#SBATCH --job-name=PCAngsd
#SBATCH --partition=short
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
IN_DIR=/data/biol-silvereye/ball6625/norfolk-island/outputs/genotype-likelihoods
OUT_DIR=/data/biol-silvereye/ball6625/norfolk-island/outputs/pca
mkdir -p $OUT_DIR

# Concatenate all beagle files
# We can do this without removing headers as these will be ignored by PCA anyway

if [ ! -f "$OUT_DIR/museum.beagle.gz" ]; then
    echo "Combining beagle files..."
    zcat $IN_DIR/chr*/*.beagle.gz | gzip > $OUT_DIR/museum.beagle.gz
fi

# Run pcangsd
pcangsd -b $OUT_DIR/museum.beagle.gz \
    --iter 1000 \
    -o $OUT_DIR/museum

