#!/bin/bash
#SBATCH --partition=ycga
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50000

# Define the input files and paths
VCF_DIR="../phased"
OUT_DIR="out"

# Create the output directory if it doesn't exist
mkdir -p $OUT_DIR

# Process the fourth population file (e.g., GBR.txt)
POP_FILE="GBR.txt"
POP_NAME=$(basename $POP_FILE .txt)

# Create a list of individuals for the current population
awk 'NR>1 {print $1"_"$1}' $POP_FILE > ${POP_NAME}.iids.txt

# Loop through chromosomes 1 to 22
for CHR in {1..22}; do
    VCF_FILE="${VCF_DIR}/chr${CHR}.all.100kga.phased_output.vcf.gz"

    # Extract the VCF for the current population and chromosome
    bcftools view -S ${POP_NAME}.iids.txt --force-samples $VCF_FILE -Oz -o ${POP_NAME}.chr${CHR}.vcf.gz

    # Index the VCF file using tabix
    tabix -p vcf ${POP_NAME}.chr${CHR}.vcf.gz

    # Define the output file for vcf2smc
    OUT_FILE="${OUT_DIR}/${POP_NAME}.chr${CHR}.smc.gz"

    # Extract individual IDs from the VCF header
    IND_LIST=$(bcftools query -l ${POP_NAME}.chr${CHR}.vcf.gz | paste -sd, -)

    # Run vcf2smc with a missing cutoff
    smc++ vcf2smc -c 1000 ${POP_NAME}.chr${CHR}.vcf.gz $OUT_FILE ${CHR} $POP_NAME:$IND_LIST
done

# Clean up the temporary file
rm ${POP_NAME}.iids.txt

echo "VCF extraction and vcf2smc conversion completed for ${POP_NAME}."
