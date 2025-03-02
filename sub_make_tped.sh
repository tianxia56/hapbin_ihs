#!/bin/bash

base_dir="/home/tx56/palmer_scratch/100kga/ihs"
pops=("IND" "MGN" "SASI" "SASP" "PJL" "KOR")
mkdir -p clean_bfile
mkdir -p clean_tped

# Loop through each population
for pop in "${pops[@]}"; do

  # Submit a single job for each population
  sbatch <<EOT
#!/bin/bash
#SBATCH --partition=ycga
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10000

base_dir="${base_dir}"
pop="${pop}"

for i in {1..22}; do
    echo "Processing chromosome \${i}..."
    clean_dir="\${base_dir}/clean_bfile"
    raw_bfile="\${base_dir}/qc_bfile/\${pop}.\${i}.qc1"
    
    if [ ! -f "\${raw_bfile}.bed" ]; then
        echo "ERROR: Raw bfile for chromosome \${i} not found!"
        exit 1
    fi
    
    # =====去重和移除. RSID=====
    # 提取出现多于1次的rsid和rsid为"."的条目
    awk '(\$2 == ".") || count[\$2]++ {print \$2}' \${raw_bfile}.bim > \${clean_dir}/tmp.exclude.list
    
    # 排除重复的rsid和.条目
    plink --bfile \${raw_bfile} \
          --exclude \${clean_dir}/tmp.exclude.list \
          --make-bed \
          --out \${clean_dir}/\${pop}.\${i}.dedup
    
    # ===== 验证并清理 dedup.bim =====
    dedup_bim="\${clean_dir}/\${pop}.\${i}.dedup.bim"
    awk '(\$2 == ".") || count[\$2]++ {print \$2}' \${dedup_bim} > \${clean_dir}/tmp.dedupcheck.exclude.list

    # If tmp.dedupcheck.exclude.list is not empty, exclude entries and redo the dedup.bim
    if [ -s \${clean_dir}/tmp.dedupcheck.exclude.list ]; then
        echo "Removing duplicates and . RSIDs from deduplicated file for chromosome \${i}..."
        plink --bfile \${clean_dir}/\${pop}.\${i}.dedup \
              --exclude \${clean_dir}/tmp.dedupcheck.exclude.list \
              --make-bed \
              --out \${clean_dir}/\${pop}.\${i}.dedup.clean
        dedup_bim="\${clean_dir}/\${pop}.\${i}.dedup.clean.bim"
    fi

    # ===== 参考等位更新 =====
    aa_bim="\${base_dir}/AA_bfile/\${pop}.\${i}.AA.bim"
    aa_ref="\${base_dir}/AA_bfile/\${pop}.\${i}.ref_alleles.txt"
    
    # 生成参考等位文件并排除. RSID
    awk 'BEGIN {OFS="\t"} NR==FNR {if(\$2 != ".") a[\$2]=\$5; next} (\$2 in a) {print \$2, a[\$2]}' \${aa_bim} \${dedup_bim} | sort -u > \${aa_ref}
    
    # 严格检查参考文件
    if [ ! -s "\${aa_ref}" ]; then
        echo "ERROR: Empty reference allele file for chr \${i}"
        exit 1
    fi
    
    plink --bfile \${clean_dir}/\${pop}.\${i}.dedup${CLEAN_SUFFIX} \
          --a1-allele \${aa_ref} \
          --make-bed \
          --out \${clean_dir}/\${pop}.\${i}.AA

    # ===== 转换为tped =====
    plink --bfile \${clean_dir}/\${pop}.\${i}.AA \
          --recode transpose \
          --cm-map \${base_dir}/rmap/\${pop}.\${i}.map.txt \${i} \
          --biallelic-only \
          --exclude /home/tx56/1KGP/exclude_rsid.txt \
          --out \${base_dir}/clean_tped/\${pop}.\${i}.AA
    
    # 清理中间文件
    rm \${clean_dir}/tmp.*
    rm \${clean_dir}/tmp.dedupcheck.exclude.list

done
EOT
done
