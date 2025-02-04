#!/bin/bash
#SBATCH --time=12:00:00               
#SBATCH --job-name=ngsadmix     
#SBATCH --partition=short
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4        
#SBATCH --output=logs/ngsadmix_%A.log      
#SBATCH --error=logs/ngsadmix_%A.error
#SBATCH --mail-type=END
#SBATCH --mail-user=abigail.williams@biology.ox.ac.uk

# Set paths
IN_DIR=/data/biol-silvereye/ball6625/norfolk-island/outputs/genotype-likelihoods
OUT_DIR=/data/biol-silvereye/ball6625/norfolk-island/outputs/ngsadmix
mkdir -p $OUT_DIR

# Concatenate all beagle files
# Ignore headers

INPUT_BEAGLE=$OUT_DIR/museum.beagle.gz

if [ ! -f $INPUT_BEAGLE ]; then
    echo "Combining beagle files..."
    (zcat $IN_DIR/chr*/*.beagle.gz | grep -m 1 "^marker"; zcat $IN_DIR/chr*/*.beagle.gz | grep -v "^marker") | gzip > $OUT_DIR/museum.beagle.gz
fi

#Run NGSadmix
/data/biol-silvereye/ball6625/software/NGSadmix \
    -likes $INPUT_BEAGLE \
    -K 3 \
    -P 4 \
    -misTol 0.1 \
    -o $OUT_DIR/results