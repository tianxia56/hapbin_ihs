import pandas as pd
import os

def convert_rmap_to_genetic_map(rmap_file, output_file, chromosome):
    # Read the RMAP file into a DataFrame
    rmap_df = pd.read_csv(rmap_file, sep='\t', header=None, names=['start', 'end', 'rate'])
    
    # Initialize the genetic map DataFrame
    genetic_map_df = pd.DataFrame(columns=['pos', 'chr', 'cM'])
    
    # Calculate the cumulative genetic distance
    cumulative_cM = 0
    rows = []
    for index, row in rmap_df.iterrows():
        start_pos = int(row['start'])
        end_pos = row['end']
        rate = row['rate']
        
        # Calculate the genetic distance for the current segment
        segment_length = end_pos - start_pos
        genetic_distance = segment_length * rate * 100
        
        # Update the cumulative genetic distance
        cumulative_cM += genetic_distance
        
        # Append the current position to the list of rows
        rows.append({'pos': start_pos, 'chr': chromosome, 'cM': cumulative_cM})
    
    # Convert the list of rows to a DataFrame
    genetic_map_df = pd.DataFrame(rows)
    
    # Ensure pos column is integer
    genetic_map_df['pos'] = genetic_map_df['pos'].astype(int)
    
    # Save the genetic map DataFrame to the output file
    genetic_map_df.to_csv(output_file, sep='\t', index=False)

# List of populations "IND", "MGN", "SASP"-cluster, "KOR"-reference
populations = ["KOR"]

# Ensure the output directory exists
output_dir = '/home/tx56/palmer_scratch/100kga/ihs/rmap/'
os.makedirs(output_dir, exist_ok=True)

# Loop through each population and chromosome
for pop in populations:
    for i in range(1, 23):
        rmap_file = f'/home/tx56/palmer_scratch/100kga/reference/{pop}.chr{i}.rmap'
        output_file = f'{output_dir}{pop}.{i}.map.txt'
        chromosome = i
        
        convert_rmap_to_genetic_map(rmap_file, output_file, chromosome)
