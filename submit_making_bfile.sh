#!/bin/bash
#SBATCH --partition=ycga
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50000

mkdir -p /home/tx56/palmer_scratch/100kga/ihs/bfile/

# List of populations to run selscan
target_pops1=("IND" "MGN" "SASI" "SASP")

# Loop through each population and chromosome
for pop in ${target_pops1[@]}
do
  for i in {1..22}
  do
    plink --vcf /home/tx56/palmer_scratch/100kga/cluster/${pop}.chr${i}.vcf.gz --make-bed --biallelic-only --exclude /home/tx56/1KGP/exclude_rsid.txt --out /home/tx56/palmer_scratch/100kga/ihs/bfile/${pop}.${i}
  done
done

# List of populations to run selscan
target_pops2=("PJL" "KOR")

# Loop through each population and chromosome
for pop in ${target_pops2[@]}
do
  for i in {1..22}
  do
    plink --vcf /home/tx56/palmer_scratch/100kga/reference/${pop}.chr${i}.vcf.gz --make-bed --biallelic-only --exclude /home/tx56/1KGP/exclude_rsid.txt --out /home/tx56/palmer_scratch/100kga/ihs/bfile/${pop}.${i}
  done
done
