#!/bin/bash

# Define populations and their corresponding individual numbers ["IND"]=89 ["SASP"]=97 ["MGN"]=102 ["SASI"]=149 
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
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=12000

# Debugging output
echo "Running job for population: $pop"
echo "VCF directory: $vcf_dir"
echo "n: $n, N: $N"

# Create lookup table
pyrho make_table -n $n -N $N --mu 1.25e-8 --logfile . --outfile ${pop}_n_${n}_N_${N}_lookuptable.hdf --approx --smcpp_file analysis_${pop}/${pop}.csv --decimate_rel_tol 0.1 --numthreads 12

# Create hyper-parameter files
pyrho hyperparam -n $n --mu 1.25e-8 --blockpenalty 50,100 --windowsize 25,50 --logfile . --tablefile ${pop}_n_${n}_N_${N}_lookuptable.hdf --num_sims 5 --smcpp_file analysis_${pop}/${pop}.csv --outfile ${pop}_hyperparam_results.txt --numthreads 12

# Loop through each chromosome
for i in {1..22}; do
  # Debugging output for each chromosome
  echo "Processing chromosome: \$i"
  echo "VCF file: ${vcf_dir}/${pop}.chr\${i}.vcf.gz"

  # Create rmap file
  pyrho optimize --tablefile ${pop}_n_${n}_N_${N}_lookuptable.hdf --vcffile ${pop}.chr\${i}.vcf.gz --outfile ${pop}.chr\${i}.rmap --blockpenalty 50 --windowsize 50 --logfile . --numthreads 12
done
EOT
done
