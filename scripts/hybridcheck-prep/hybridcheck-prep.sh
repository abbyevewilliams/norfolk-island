#!/bin/bash
#SBATCH --time=12:00:00               
#SBATCH --job-name=hybridcheck_prep    
#SBATCH --partition=short
#SBATCH --mem=32G      
#SBATCH --error=logs/hybridcheck_prep_%A.error
#SBATCH --output=logs/hybridcheck_prep_%A.log

# Load plink
ml PLINK/1.9b_6.21-x86_64

# Set paths
IN_DIR=/data/biol-silvereye/ball6625/norfolk-island/outputs/genotype-likelihoods
OUT_DIR=/data/biol-silvereye/ball6625/norfolk-island/outputs/hybridcheck-prep
mkdir -p $OUT_DIR

# Concat beagle files
if [ ! -f "$OUT_DIR/museum.BeagleVersion.vcf.gz" ]; then
    echo "Combining beagle files..."
    (zcat $IN_DIR/chr*/*.BeagleVersion.vcf.gz | grep -m 1 "^marker"; zcat $IN_DIR/chr*/*.BeagleVersion.vcf.gz | grep -v "^marker") | gzip > $OUT_DIR/museum.BeagleVersion.vcf.gz
fi

echo "Converting beagle to plink..."

# Convert beagle file to plink format
plink --vcf $OUT_DIR/museum.BeagleVersion.vcf.gz --make-bed --out $OUT_DIR/plink_format

# Extract SNPs for hybridcheck
plink --bfile $OUT_DIR/plink_format --recode transpose --out $OUT_DIR/plink_transpose

#echo "Activating python env for conversion..."

# Load  python env for conversion (already created)
ml Anaconda3
export CONPREFIX=/data/biol-silvereye/ball6625/norfolk-island/scripts/hybridcheck-prep/myenv
source activate $CONPREFIX

#echo "Converting plink to fasta..."

# Run python conversion script
python3 plink2fasta.py --tped $OUT_DIR/plink_transpose.tped --tfam $OUT_DIR/plink_transpose.tfam --output $OUT_DIR/output.fasta
