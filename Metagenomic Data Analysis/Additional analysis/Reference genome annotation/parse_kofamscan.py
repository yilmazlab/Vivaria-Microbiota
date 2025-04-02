import os
import pandas as pd
from collections import defaultdict

def parse_kofamscan(floc):
    """
    v1.1: 5/1/2023
    - update "KO_definition" to be 5

    v1.0: 1/6/2023

    Parse kofamscan results. Only save results where KO score > threshold

    Returns:
        Adb: DataFrame with KOfam results
    """
    table = defaultdict(list)

    with open(floc, 'r') as o:

        # This block is for RAM efficiency
        while True:
            line = o.readline()
            if not line:
                break

            line = line.strip()
            if line[0] == '#':
                continue

            lw = line.split()
            if lw[0] == '*':
                del lw[0]

            if lw[2] == '-':
                lw[2] = 0

            try:
                if float(lw[3]) >= float(lw[2]):
                    g = lw[0]
                    k = lw[1]

                    table['gene'].append(g)
                    table['KO'].append(k)
                    table['thrshld'].append(float(lw[2]))
                    table['score'].append(float(lw[3]))
                    table['e_value'].append(float(lw[4]))
                    table['KO_definition'].append(' '.join(lw[5:]))
            except:
                print(line)
                assert False
    o.close()

    Adb = pd.DataFrame(table)
    return Adb

# Define the input directory
input_dir = "/storage/homefs/ib21q318/data/wild_mice/genomes/Bacteroides_acidifaciens/"
input_file = "genes.faa.kofamscan"
floc = os.path.join(input_dir, input_file)

# Parse the KofamScan file
Adb = parse_kofamscan(floc)

# Save the DataFrame to a CSV file
output_csv = "/storage/homefs/ib21q318/data/wild_mice/genomes/Bacteroides_acidifaciens/kofamscan_results.csv"
Adb.to_csv(output_csv, index=False)
