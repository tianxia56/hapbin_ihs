#!/bin/bash

# Define populations and their corresponding individual numbers ["sasi"]=150 ["PJL"]=96 
declare -A populations=( ["SASI"]=90 )

# Loop through each population
for pop in "${!populations[@]}"; do
  ind=${populations[$pop]}
  n=$((ind * 2))
  N=$(echo "($n * 1.25)/1" | bc)

  # Submit a single job for each population
  sbatch <<EOT
#!/bin/bash
#SBATCH --partition=week
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000

# Load BCFtools module
#module load BCFtools/1.11-GCC-10.2.0

# Create the list of individuals in the required format
#awk '{print $1"_"$1}' /home/tx56/palmer_scratch/100kga/PCA/${pop}.filtered.txt > ${pop}_individuals.txt

# Create the directory for the new VCF files
mkdir -p sasi90

# Loop through each chromosome
for i in {1..22}; do
  # Debugging output for each chromosome
  echo "Processing chromosome: \$i"

  # Extract the subset of individuals using bcftools
  bcftools view -S ${pop}_individuals.txt -Oz -o sasi90/${pop}.90.chr\${i}.vcf.gz ${pop}.chr\${i}.vcf.gz

  # Create rmap file
  pyrho optimize --tablefile ${pop}_n_${n}_N_${N}_lookuptable.hdf --vcffile sasi90/${pop}.90.chr\${i}.vcf.gz --outfile ${pop}.chr\${i}.rmap --blockpenalty 50 --windowsize 50 --logfile . --numthreads 12
done
EOT
done
