#!/bin/bash
#SBATCH --time=24:00:00               
#SBATCH --job-name=dsuite     
#SBATCH --partition=medium
#SBATCH --mem=24G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8        
#SBATCH --output=logs/dsuite_%A.log      
#SBATCH --error=logs/dsuite_%A.error
#SBATCH --mail-type=END
#SBATCH --mail-user=abigail.williams@biology.ox.ac.uk

# Set paths
export PATH=/data/biol-silvereye/ball6625/software/Dsuite/Build:$PATH
IN_PATH=/data/biol-silvereye/ball6625/norfolk-island/outputs/genotype-likelihoods
OUT_PATH=/data/biol-silvereye/ball6625/norfolk-island/outputs/dsuite
CHROMS=$(cat chrom.list)
mkdir -p $OUT_PATH

# Loop over chromosomes
for CHROM in $CHROMS
do
	# Logging
	echo "Starting analysis for chromosome: ${CHROM}"

	# Mkdir and move into it
	mkdir -p ${OUT_PATH}/${CHROM}
	
	# Assign VCF
	VCF=${IN_PATH}/${CHROM}/*.BeagleVersion.vcf.gz

	# Run Dtrios
	Dsuite Dtrios \
		${VCF} \
		SETS.txt \
		-o ${OUT_PATH}/${CHROM}/Dsuite_${CHROM}
	
	# Run Dinvestigate
	Dsuite Dinvestigate \
    ${VCF} \
    SETS.txt \
	test_trios.txt \
	-n ${OUT_PATH}/${CHROM}/Dsuite_${CHROM}

done

echo "Individual chromosomes analysed."

# Once complete, combine the results across chromosomes

RESULT_FILES=$(find ${OUT_PATH} -type d -name 'chr*' -exec find {} -type f -name '*.Dmin.txt' \; | sort)
Dsuite Dcombine ${OUT_PATH}/combined_results ${RESULT_FILES}

echo "Results combined."