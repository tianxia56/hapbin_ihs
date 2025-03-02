import pandas as pd

# Load the CSV file
file_path = "/home/tx56/palmer_scratch/100kga/ADMIXTURE/extrapop/admixture_data_k6_1kgp_and_100kga.csv"
df = pd.read_csv(file_path)

# Extract IID and POP columns
iid_pop_df = df[['IID', 'POP']]

# Filter and save the data for each specified population
populations = ['YRI', 'GBR', 'KOR', 'PJL']
for pop in populations:
    pop_df = iid_pop_df[iid_pop_df['POP'] == pop]
    pop_df.to_csv(f"{pop}.txt", index=False, header=False, sep='\t')

print("Files saved successfully.")
