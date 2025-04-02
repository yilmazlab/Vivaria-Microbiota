import pandas as pd
from collections import defaultdict

def parse_Pfam(floc, aloc):
    """
    Parse Pfam results.

    Args:
        floc (str): Path to the filtered Pfam hits file.
        aloc (str): Path to the Pfam annotation info file.

    Returns:
        pd.DataFrame: DataFrame with annotated Pfam results.
    """
    # Read the filtered Pfam results
    PFdb = pd.read_csv(floc, header=1, sep='\s+', 
                       names=['query-id', 'match-id', 'score', 'boundaries', 
                              'resolved', 'cond-evalue', 'indp-evalue', 'junk'], 
                       index_col=None)
    # Remove unnecessary column
    del PFdb['junk']

    # Compress the results by removing duplicates
    Ofdb = PFdb.rename(columns={'query-id':'gene', 'match-id':'pFam'})
    Ofdb = Ofdb.sort_values('score').drop_duplicates(subset=['gene', 'pFam'], keep='last').sort_values('gene')
    Ofdb['pFam'] = Ofdb['pFam'].astype('category')

    # Add more annotation data
    PFDdb = pd.read_csv(aloc)
    TDdb = pd.merge(Ofdb, PFDdb[['NAME', 'DESC', 'ACC']].rename(columns={'NAME':'pFam'}), on='pFam')
    TDdb = TDdb.drop_duplicates()

    return TDdb

def make_Pfam_info(LOCATION):
    """
    Parse the Pfam-A.hmm file to extract annotations.

    Args:
        LOCATION (str): Path to the Pfam-A.hmm file.

    Returns:
        pd.DataFrame: DataFrame with Pfam annotations.
    """
    a2c = {}
    grab_next = False
    with open(LOCATION) as f:
        for line in f.readlines():
            line = line.strip()
            if grab_next:
                assert line.split()[0] == 'ACC', line
                a2c[acc] = ' '.join(line.split()[1:])
                grab_next = False

            if line[:4] == 'NAME':
                acc = line.split()[1]
                grab_next = True
    Ddb = pd.DataFrame(list(a2c.items()))
    Ddb.rename(columns={0:'NAME', 1:'ACC'}, inplace=True)

    a2c = {}
    grab_next = False
    with open(LOCATION) as f:
        for line in f.readlines():
            line = line.strip()
            if grab_next:
                if line.split()[0] == 'NC':
                    a2c[acc] = float(line.split()[1])
                    grab_next = False

            if line[:3] == 'ACC':
                acc = line.split()[1]
                grab_next = True
    Ndb = pd.DataFrame(list(a2c.items()))
    Ndb.rename(columns={0:'ACC', 1:'NC'}, inplace=True)

    a2c = {}
    grab_next = False
    with open(LOCATION) as f:
        for line in f.readlines():
            line = line.strip()
            if grab_next:
                assert line.split()[0] == 'DESC', line
                a2c[acc] = ' '.join(line.split()[1:])
                grab_next = False

            if line[:3] == 'ACC':
                acc = line.split()[1]
                grab_next = True
    DEdb = pd.DataFrame(list(a2c.items()))
    DEdb.rename(columns={0:'ACC', 1:'DESC'}, inplace=True)

    PFdb = pd.merge(Ddb, Ndb).merge(DEdb)
    return PFdb

if __name__ == "__main__":
    # Define file locations
    hmm_file = "/storage/homefs/ib21q318/Pfam-A.hmm"
    filtered_hits_file = "/storage/homefs/ib21q318/data/wild_mice/genomes/Bacteroides_acidifaciens/GUT_GENOME000221_reps.genes.faa_vs_Pfam.hmm.filtered.txt"  # Updated file name
    annotation_output_file = "/storage/homefs/ib21q318/data/wild_mice/genomes/Bacteroides_acidifaciens/Pfam-A.info.csv"
    final_output_file = "/storage/homefs/ib21q318/data/wild_mice/genomes/Bacteroides_acidifaciens/Pfam_results_annotated.csv"

    # Create the Pfam annotation file
    PIdb = make_Pfam_info(hmm_file)
    PIdb.to_csv(annotation_output_file, index=False)

    # Parse Pfam results and annotate
    Hdb = parse_Pfam(filtered_hits_file, annotation_output_file)

    # Save the final annotated results
    Hdb.to_csv(final_output_file, index=False)

    print(f"Annotated Pfam results saved to {final_output_file}")