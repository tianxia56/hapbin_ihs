import pandas as pd
import numpy as np
from scipy import stats

# Load the data from SASI.pcs.txt
data = pd.read_csv('SASI.pcs.txt', sep='\t')

# Ensure each population has at least 2 individuals
pop_counts = data['POP'].value_counts()
valid_pops = pop_counts[pop_counts >= 2].index
filtered_data = data[data['POP'].isin(valid_pops)]

# Calculate Z-scores for PC1 and PC2
z_scores = np.abs(stats.zscore(filtered_data[['PC1', 'PC2']]))

# Keep rows where Z-scores are less than a threshold (e.g., 3)
filtered_data_no_outliers = filtered_data[(z_scores < 3).all(axis=1)]

# Calculate the number of rows to keep from each population proportionally
total_rows = 92
pop_proportions = filtered_data_no_outliers['POP'].value_counts(normalize=True)
rows_to_keep = (pop_proportions * total_rows).astype(int)

# Ensure at least 2 individuals are kept for each population
rows_to_keep[rows_to_keep < 2] = 2

# Adjust the total number of rows to keep to be exactly 90
while rows_to_keep.sum() > total_rows:
    rows_to_keep[rows_to_keep.idxmax()] -= 1

while rows_to_keep.sum() < total_rows:
    rows_to_keep[rows_to_keep.idxmin()] += 1

# Proportionally remove outliers from each population
filtered_data_list = []
for pop in rows_to_keep.index:
    pop_data = filtered_data_no_outliers[filtered_data_no_outliers['POP'] == pop]
    if len(pop_data) > rows_to_keep[pop]:
        pop_data = pop_data.sample(rows_to_keep[pop], random_state=42)
    filtered_data_list.append(pop_data)

final_filtered_data = pd.concat(filtered_data_list)

# Write the filtered data to SASI.filtered.txt without header
final_filtered_data.to_csv('SASI.filtered.txt', sep='\t', index=False, header=False)

print("Filtered data saved to 'SASI.filtered.txt' without header.")
