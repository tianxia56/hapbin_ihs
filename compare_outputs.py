import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import zipfile
import os

# 配置参数
INPUT_DIR = "/vast/palmer/pi/reilly/tian/100kga"  # 更新后的ZIP文件目录
OUTPUT_DIR = "./figs"  # 结果输出目录
POPULATIONS = ["IND", "KOR", "MGN", "PJL", "SASI", "SASP"]  # 需要分析的群体

def load_ihs_data(pop):
    """从ZIP文件加载标准化iHS数据"""
    zip_path = os.path.join(INPUT_DIR, f"{pop}.norm.ihs.zip")
    with zipfile.ZipFile(zip_path) as z:
        with z.open(f"{pop}.norm.ihs") as f:
            return pd.read_csv(f, sep='\t')

def calculate_stats(df):
    """计算关键统计指标"""
    stats_dict = {}
    for col in ['hapbin_std_iHS', 'GW_DAF_bin_iHS', 'JV_rmap_adj_norm_iHS']:
        stats_dict[col] = {
            'mean': np.mean(df[col]),
            'std': np.std(df[col], ddof=1),
            'skewness': stats.skew(df[col]),
            'kurtosis': stats.kurtosis(df[col]),
            '>2sd (%)': np.mean(np.abs(df[col]) > 2)*100
        }
    return pd.DataFrame(stats_dict).T

def plot_comparison(df, pop, output_dir):
    """生成群体特异的对比可视化"""
    plt.figure(figsize=(18, 12))
    plt.suptitle(f"iHS Standardization Comparison - {pop}", y=1.02, fontsize=40)
    
    # 直方图分布对比
    colors = ['blue', 'green', 'red']
    for col, color in zip(df.columns[-3:], colors):
        plt.hist(df[col], bins=50, density=True, alpha=0.5, label=col, color=color)
    plt.xlabel('iHS Value', fontsize=40)
    plt.ylabel('Density', fontsize=40)
    plt.title('Distribution Comparison', fontsize=40)
    plt.legend(fontsize=20)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{pop}_comparison_hist.png"), dpi=300, bbox_inches='tight')
    plt.close()

def plot_binned_qqplots(df, pop, output_dir, n_cols=4):
    """按DAF分箱绘制QQ图对比"""
    # 创建分箱，确保包含0到1范围内的所有分箱
    daf_bins = np.arange(0, 1.05, 0.1)
    df['DAF_bin'] = pd.cut(df['Freq'], bins=daf_bins, right=False, include_lowest=True)
    
    # 分箱绘图
    n_bins = len(daf_bins) - 1
    n_rows = (n_bins + n_cols - 1) // n_cols  # 计算行数

    plt.figure(figsize=(6 * n_cols, 6 * n_rows))
    
    for i, (daf_bin, sub_df) in enumerate(df.groupby('DAF_bin', observed=True), 1):
        plt.subplot(n_rows, n_cols, i)
        
        # 绘制三组QQ线
        colors = ['blue', 'green', 'red']
        for col, color in zip(['hapbin_std_iHS', 'GW_DAF_bin_iHS', 'JV_rmap_adj_norm_iHS'], colors):
            if sub_df[col].count() > 0:
                (osm, osr), _ = stats.probplot(sub_df[col], dist="norm")
                plt.plot(osm, osr, 'o', markersize=3, alpha=0.5, color=color, label=col)
        
        plt.plot([-3,3], [-3,3], 'k--', alpha=0.5)
        plt.title(f"DAF: {daf_bin}\n(n={len(sub_df):,})", fontsize=20)
        plt.xlabel('Theoretical Quantiles', fontsize=20)
        plt.ylabel('Ordered Values', fontsize=20)
        plt.legend(fontsize=10)
    
    plt.suptitle(f"DAF-Binned QQ Plots - {pop}", y=1.02, fontsize=20)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{pop}_daf_binned_qq.png"), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)  # Ensure the output directory exists
    
    for pop in POPULATIONS:
        print(f"Processing {pop}...")
        
        # Load data
        pop_df = load_ihs_data(pop)
        
        # Calculate statistics
        stats_df = calculate_stats(pop_df)
        stats_df.to_csv(os.path.join(OUTPUT_DIR, f"{pop}_statistics.tsv"), sep='\t')
        
        # Generate visualizations
        plot_comparison(pop_df, pop, OUTPUT_DIR)
        plot_binned_qqplots(pop_df, pop, OUTPUT_DIR, n_cols=4)
        
        # Save correlation matrices
        corr_methods = ['pearson', 'spearman']
        for method in corr_methods:
            corr_df = pop_df[['hapbin_std_iHS', 'GW_DAF_bin_iHS', 'JV_rmap_adj_norm_iHS']].corr(method=method)
            corr_df.to_csv(os.path.join(OUTPUT_DIR, f"{pop}_{method}_corr.tsv"), sep='\t')

if __name__ == "__main__":
    main()
