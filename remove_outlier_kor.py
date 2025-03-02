import pandas as pd
import numpy as np
from scipy import stats

# Load the data from KOR_HAN_JPN.pcs.txt
data = pd.read_csv('KOR_HAN_JPN.pcs.txt', sep='\t')

# Filter to keep only KOR population
kor_data = data[data['POP'] == 'KOR']

# Calculate Z-scores for PC1 and PC2
z_scores = np.abs(stats.zscore(kor_data[['PC1', 'PC2']]))

# Keep rows where Z-scores are less than a threshold (e.g., 3)
kor_data_no_outliers = kor_data[(z_scores < 3).all(axis=1)]

# Ensure there are at least 90 rows
if len(kor_data_no_outliers) < 90:
    raise ValueError("Not enough KOR population data to keep 90 rows after removing outliers.")

# Randomly sample 90 rows from KOR population without outliers
kor_filtered = kor_data_no_outliers.sample(90, random_state=42)

# Write the filtered data to KOR.filtered.txt
kor_filtered.to_csv('KOR.filtered.txt', sep='\t', index=False)

print("Filtered KOR data saved to 'KOR.filtered.txt'.")
