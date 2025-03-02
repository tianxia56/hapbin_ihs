import os

def extract_and_clean_columns(input_file, output_file_map, output_file_hap, allele_map):
    with open(input_file, 'r') as file, open(output_file_map, 'w') as map_file, open(output_file_hap, 'w') as hap_file:
        for line in file:
            columns = line.split()
            if len(columns) < 5:
                continue

            rsid = columns[1]
            ref_allele = allele_map.get(rsid)

            # Extract marker information (first four columns)
            map_file.write(f"{columns[0]}\t{columns[1]}\t{columns[2]}\t{columns[3]}\n")

            # Clean genotype columns
            cleaned_columns = ['0' if col.strip() == ref_allele else '1' if col.strip() in 'ATCG' else col for col in columns[4:]]
            hap_file.write(' '.join(cleaned_columns) + '\n')

def create_map_file(original_map_file):
    temp_file = original_map_file + ".tmp"
    with open(original_map_file, 'r') as map_file, open(temp_file, 'w') as new_file:
        for line in map_file:
            columns = line.split()
            if len(columns) != 4:
                continue
            # Remove any non-printable characters and ensure proper formatting
            formatted_columns = [col.encode('ascii', 'ignore').decode('ascii') for col in columns]
            cleaned_line = ' '.join(formatted_columns).strip()
            new_file.write(cleaned_line + '\n')
    os.replace(temp_file, original_map_file)

def main():
    target_pops = ["IND", "MGN", "PJL", "SASP", "SASI", "KOR"]
    path = "/home/tx56/palmer_scratch/100kga/ihs/clean_tped"

    for pop in target_pops:
        for i in range(1, 23):
            input_file = f"{path}/{pop}.{i}.AA.tped"
            
            if not os.path.exists(input_file):
                print(f"Skipping {pop} chromosome {i} as {input_file} does not exist.")
                continue

            ref_allele_file = f"/home/tx56/palmer_scratch/100kga/ihs/AA_bfile/{pop}.{i}.ref_alleles.txt"
            if not os.path.exists(ref_allele_file):
                print(f"Reference allele file {ref_allele_file} does not exist.")
                continue

            allele_map = {}
            with open(ref_allele_file, 'r') as ref_file:
                for line in ref_file:
                    rsid, allele = line.split()
                    allele_map[rsid] = allele

            os.makedirs("/home/tx56/palmer_scratch/100kga/ihs/hap", exist_ok=True)
            output_file_hap = f"/home/tx56/palmer_scratch/100kga/ihs/hap/{pop}.{i}.hap"
            output_file_map = f"/home/tx56/palmer_scratch/100kga/ihs/hap/{pop}.{i}.map"
            extract_and_clean_columns(input_file, output_file_map, output_file_hap, allele_map)
            create_map_file(output_file_map)
            print(f"Extracted columns saved to {output_file_map} and {output_file_hap}")

if __name__ == "__main__":
    main()
