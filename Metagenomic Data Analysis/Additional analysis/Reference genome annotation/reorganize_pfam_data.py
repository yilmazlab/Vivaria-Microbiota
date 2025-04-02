import re
import pandas as pd

# Define the file paths
input_file_path = '/storage/homefs/ib21q318/data/wild_mice/genomes/Bacteroides_acidifaciens/GUT_GENOME000221_reps.genes.faa_vs_Pfam.hmm'  # replace with your actual file path
output_file_path = '/storage/homefs/ib21q318/data/wild_mice/genomes/Bacteroides_acidifaciens/GUT_GENOME000221_reps.genes.faa_vs_Pfam.hmm.organized.csv'  # this will be the output file

# Function to extract data from a line and organize it into columns
def extract_columns(line):
    # Split the line into parts using spaces, but respect the embedded data fields by splitting manually
    parts = line.split()
    
    # Extracting the locations and ID from the end of the line using regex
    location_data = re.search(r'# (\d+) # (\d+) # ([-\d]+) # ID=(.*?);', line)
    
    if location_data:
        start_location = location_data.group(1)
        end_location = location_data.group(2)
        strand = location_data.group(3)
        gene_id = location_data.group(4)
        
        # Reorganize into a new line with extracted columns
        new_line = parts[:22] + [start_location, end_location, strand, gene_id]
        return new_line
    else:
        return parts

# Read and process the file
with open(input_file_path, 'r') as file:
    content = file.readlines()

processed_lines = []
for line in content:
    if not line.startswith('#'):  # Ignore header and comment lines
        processed_lines.append(extract_columns(line))

# Define the column names
columns = ['target_name', 'accession', 'tlen', 'query_name', 'query_accession', 'qlen', 
           'E-value', 'score', 'bias', 'of', 'c-Evalue', 'i-Evalue', 'score_domain', 
           'bias_domain', 'hmm_from', 'hmm_to', 'ali_from', 'ali_to', 'env_from', 
           'env_to', 'acc', 'description', 'start_location', 'end_location', 'strand', 'gene_id']

# Convert the processed lines into a DataFrame
df = pd.DataFrame(processed_lines, columns=columns)

# Save the DataFrame to a CSV file
df.to_csv(output_file_path, index=False)

print(f"Data has been reorganized and saved to {output_file_path}")