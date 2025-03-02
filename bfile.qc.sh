#!/bin/bash

base_dir="/home/tx56/palmer_scratch/100kga/ihs"
qc_dir="${base_dir}/qc_bfile"
tmp_dir="${base_dir}/tmp"

# Create the necessary directories if they don't exist
mkdir -p ${qc_dir}
mkdir -p ${tmp_dir}

pops=("IND" "MGN" "SASI" "SASP" "PJL" "KOR")
chr=( {1..22} )

# Loop through each population and chromosome
for pop in "${pops[@]}"; do
  for i in "${chr[@]}"; do
    echo "Processing ${pop}.${i}..."

    raw_bfile="${base_dir}/bfile/${pop}.${i}"
    qc_bfile="${qc_dir}/${pop}.${i}.qc1"
    duplist_file="${tmp_dir}/${pop}.${i}.duplist"
    rs_exclude_list="${tmp_dir}/${pop}.${i}.exclude.list"
    bim_tmp_file="${tmp_dir}/${pop}.${i}.bim.tmp"

    if [ ! -f "${raw_bfile}.bed" ]; then
      echo "ERROR: Raw bfile for chromosome ${i} not found!"
      continue
    fi

    # Remove RSIDs that do not start with "rs" and genotypes that are not A, T, C, G
    awk '($2 ~ /^rs/) && ($5 ~ /^[ATCG]$/) && ($6 ~ /^[ATCG]$/) {print $0}' ${raw_bfile}.bim > ${bim_tmp_file}

    # Ensure unique RSIDs
    awk '{if (!seen[$2]++) print}' ${bim_tmp_file} > ${qc_bfile}.bim

    # Generate the cleaned bfile using PLINK
    plink --bfile ${raw_bfile} --extract ${qc_bfile}.bim --make-bed --out ${qc_bfile}

    # Double-check for duplicate RSIDs
    awk '{count[$2]++} END {for (rsid in count) if (count[rsid] > 1) print rsid}' ${qc_bfile}.bim > ${duplist_file}

    # Check if duplist is not empty and handle it
    if [ -s ${duplist_file} ]; then
      echo "Duplicates found in ${pop}.${i}. Writing to ${duplist_file}"
      plink --bfile ${qc_bfile} \
            --exclude ${duplist_file} \
            --make-bed \
            --out ${qc_bfile}.clean

      # Replace the QC file with the cleaned version
      mv ${qc_bfile}.clean.bim ${qc_bfile}.bim
      mv ${qc_bfile}.clean.bed ${qc_bfile}.bed
      mv ${qc_bfile}.clean.fam ${qc_bfile}.fam
      rm ${qc_bfile}.clean.log
    else
      echo "No duplicates found in ${pop}.${i}."
      rm ${duplist_file}
    fi

    # Cleaning up temp files
    rm ${bim_tmp_file}
    rm ${rs_exclude_list}
  done
done
