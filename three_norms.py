import pandas as pd
import numpy as np
import os
import zipfile
from glob import glob

# ========================
# 配置参数
# ========================
NUM_DAF_BINS = 50
MIN_SNPS_PER_BIN = 2
COMPRESS_LEVEL = 9
DEFAULT_FILL_VALUE = 0.0
INPUT_DIR = "output"
BIM_DIR = "hap"
OUTPUT_DIR = "normed_ihs"  # mkdir -p

IHS_COLUMN_NAMES = [
    'Index', 'ID', 'Freq',
    'iHH_0', 'iHH_1', 'iHS',
    'Std iHS'
]

# ========================
# 数据加载函数
# ========================
def load_bim_positions(pop, chr_num):
    """加载.bim文件获取物理位置和重组率"""
    bim_path = os.path.join(BIM_DIR, f"{pop}.{chr_num}.map")
    return pd.read_csv(
        bim_path,
        sep=r'\s+',
        header=None,
        usecols=[0, 1, 2, 3],
        names=['chr', 'ID', 'cM', 'pos'],
        dtype={'chr': 'uint8', 'ID': 'string', 'cM': 'float32', 'pos': 'uint32'},
        engine='python'
    )

# ========================
# 标准化计算函数（修复版）
# ========================
def calculate_global_norm(full_df):
    """全基因组DAF分箱标准化"""
    bins = np.linspace(0, 1, NUM_DAF_BINS + 1)
    full_df['DAF_bin_global'] = pd.cut(
        full_df['Freq'], bins, right=False, include_lowest=True
    )
    
    # 修正聚合语法
    global_stats = full_df.groupby('DAF_bin_global', observed=False)['iHS'].agg(
        global_mean=lambda x: np.nanmean(x),
        global_std=lambda x: np.nanstd(x, ddof=1),
        count='size'
    ).reset_index()
    
    # 处理小样本分箱
    global_stats['adj_global_mean'] = np.where(
        global_stats['count'] >= MIN_SNPS_PER_BIN,
        global_stats['global_mean'],
        full_df['iHS'].mean()
    )
    global_stats['adj_global_std'] = np.where(
        global_stats['count'] >= MIN_SNPS_PER_BIN,
        global_stats['global_std'],
        full_df['iHS'].std(ddof=1)
    )
    
    merged = pd.merge(
        full_df,
        global_stats[['DAF_bin_global', 'adj_global_mean', 'adj_global_std']],
        on='DAF_bin_global',
        how='left'
    )
    
    merged['GW_DAF_bin_iHS'] = (
        (merged['iHS'] - merged['adj_global_mean']) / merged['adj_global_std']
    ).fillna(DEFAULT_FILL_VALUE)
    
    return merged.drop(columns=['DAF_bin_global'])

# ========================
# JV500bin 标准化计算函数
# ========================
def JV500bin(full_df):
    """新的标准化计算方法"""
    # Initial Standardization by Derived Allele Frequency
    bins = np.linspace(0, 1, NUM_DAF_BINS + 1)
    full_df['DAF_bin'] = pd.cut(
        full_df['Freq'],
        bins=bins,
        right=False,
        include_lowest=True
    )
    
    # Further Binning by Local Recombination Rate
    recomb_bins = np.linspace(full_df['cM'].min(), full_df['cM'].max(), 11)
    full_df['recomb_bin'] = pd.cut(
        full_df['cM'],
        bins=recomb_bins,
        right=False
    )
    
    # Create 500 bins (50 allele frequency bins × 10 recombination rate bins)
    full_df['combined_bin'] = full_df['DAF_bin'].astype(str) + '_' + full_df['recomb_bin'].astype(str)
    
    # Standardize within each of the 500 bins
    bin_stats = full_df.groupby('combined_bin')['iHS'].agg(
        mean_adj=lambda x: np.nanmean(x),
        std_adj=lambda x: np.nanstd(x, ddof=1)
    ).reset_index()
    
    global_mean = full_df['iHS'].mean()
    global_std = full_df['iHS'].std(ddof=1)
    bin_stats['mean_adj'] = bin_stats['mean_adj'].fillna(global_mean)
    bin_stats['std_adj'] = bin_stats['std_adj'].fillna(global_std)
    
    merged = pd.merge(
        full_df,
        bin_stats,
        on='combined_bin',
        how='left'
    )
    
    merged['JV_rmap_adj_norm_iHS'] = (
        (merged['iHS'] - merged['mean_adj']) / merged['std_adj']
    ).fillna(DEFAULT_FILL_VALUE)
    
    return merged.drop(columns=['DAF_bin', 'recomb_bin', 'combined_bin'])

# ========================
# 主处理函数
# ========================
def process_population(pop):
    try:
        all_chr_dfs = []
        chr_files = glob(os.path.join(INPUT_DIR, f"{pop}.?.ihs")) + \
                    glob(os.path.join(INPUT_DIR, f"{pop}.??.ihs"))
        
        for fpath in chr_files:
            base_name = os.path.basename(fpath)
            chr_num = int(base_name.split('.')[-2])
            
            ihs_df = pd.read_csv(
                fpath,
                sep=r'\s+',
                names=IHS_COLUMN_NAMES,
                header=0,
                dtype={'ID': 'string', 'pos': 'uint32'},
                engine='python'
            )
            
            bim_df = load_bim_positions(pop, chr_num)
            merged_df = pd.merge(ihs_df, bim_df, on='ID', how='inner')
            merged_df['chr'] = chr_num
            all_chr_dfs.append(merged_df)

        full_df = pd.concat(all_chr_dfs, ignore_index=True)
        full_df = calculate_global_norm(full_df)
        full_df = JV500bin(full_df)
        
        output_cols = [
            'chr', 'pos', 'ID', 'Freq',
            'iHH_0', 'iHH_1', 'iHS',
            'Std iHS', 'GW_DAF_bin_iHS', 'JV_rmap_adj_norm_iHS'
        ]
        output_df = full_df[output_cols].rename(columns={
            'Std iHS': 'hapbin_std_iHS',
            'pos': 'Position'
        })
        
        zip_path = os.path.join(OUTPUT_DIR, f"{pop}.norm.ihs.zip")
        with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=COMPRESS_LEVEL) as zf:
            with zf.open(f"{pop}.norm.ihs", 'w') as zf_out:
                output_df.to_csv(
                    zf_out,
                    sep='\t',
                    index=False,
                    float_format='%.6f',
                    encoding='utf-8'
                )
        
        print(f"✅ {pop} processed: {len(output_df):,} SNPs")
        
    except Exception as e:
        print(f"❌ Error processing {pop}: {str(e)}")
        raise

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    all_pops = sorted({os.path.basename(f).split('.')[0] for f in glob(os.path.join(INPUT_DIR, "*.ihs"))})
    
    print(f"🔍 Detected {len(all_pops)} populations")
    for i, pop in enumerate(all_pops, 1):
        print(f"🚀 Processing ({i}/{len(all_pops)}) {pop}")
        process_population(pop)
    
    print(f"\n🎉 All results saved to: {os.path.abspath(OUTPUT_DIR)}")
