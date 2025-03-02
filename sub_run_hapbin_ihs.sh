#!/bin/bash

# Define populations 
populations=( "IND" "MGN" "PJL" "SASP" "SASI" "KOR" )

# Loop through each population
for pop in "${populations[@]}"; do

  # Submit a single job for each population
  sbatch <<EOT
#!/bin/bash
#SBATCH --partition=ycga
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000
#SBATCH --job-name=ihs_${pop}
#SBATCH --output=/home/tx56/palmer_scratch/100kga/ihs/logs/${pop}_%j.log
#SBATCH --error=/home/tx56/palmer_scratch/100kga/ihs/logs/${pop}_%j.err

# Create output and log directories if they don't exist
mkdir -p /home/tx56/palmer_scratch/100kga/ihs/output
mkdir -p /home/tx56/palmer_scratch/100kga/ihs/logs

# Loop through each chromosome
for i in {1..22}
do
  /home/tx56/hapbin/build/ihsbin --hap /home/tx56/palmer_scratch/100kga/ihs/hap/${pop}.\${i}.hap --map /home/tx56/palmer_scratch/100kga/ihs/hap/${pop}.\${i}.map --out /home/tx56/palmer_scratch/100kga/ihs/output/${pop}.\${i}.ihs --cutoff 0.05
done
EOT

done
