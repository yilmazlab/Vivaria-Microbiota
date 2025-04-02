import os
import pandas as pd
import json
from collections import defaultdict

def parse_card(floc, jloc=None):
    """
    v1.0 - 1/6/2023

    Parse CARD

    Returns:
        Rdb: DataFrame with CARD results
    """
    h = ['gene', 'target', 'percentID', 'alignment_length', 'mm', 'gaps',
         'querry_start', 'querry_end', 'target_start', 'target_end', 'e-value', 'bit_score',
         'extra']
    db = pd.read_csv(floc, sep='\t', names=h)
    del db['extra']

    db['protein_seq_accession'] = [t.split('|')[1] for t in db['target']]
    db['ARO'] = [t.split('|')[2].split(':')[-1] for t in db['target']]
    db['CARD_short_name'] = [t.split('|')[3].split(':')[-1] for t in db['target']]

    # Reorder
    header = ['gene', 'CARD_short_name', 'ARO', 'target']
    db = db[header + [x for x in list(db.columns) if x not in header]]

    if jloc is None:
        return db

    # Parse additional JSON data
    j = json.load(open(jloc))

    aro2name = {}
    aro2categories = {}

    for n, m2t in j.items():
        if type(m2t) != dict:
            continue

        if 'ARO_description' in m2t:
            aro2name[m2t['ARO_accession']] = m2t['ARO_description']

        if 'ARO_category' in m2t:
            cats = []
            for cat, c2t in m2t['ARO_category'].items():
                if 'category_aro_accession' in c2t:
                    cats.append(c2t['category_aro_accession'])
            aro2categories[m2t['ARO_accession']] = '|'.join(cats)

    db['ARO_description'] = db['ARO'].map(aro2name)
    db['ARO_category_accessions'] = db['ARO'].map(aro2categories)

    header = ['gene', 'CARD_short_name', 'ARO', 'ARO_description', 'ARO_category_accessions', 'target']
    db = db[header + [x for x in list(db.columns) if x not in header]]

    return db


floc = os.path.join("/storage/homefs/ib21q318/data/wild_mice/genomes/Bacteroides_acidifaciens/", "GUT_GENOME000221genes.faa_vs_CARD.dm")
jloc = os.path.join("/storage/homefs/ib21q318/data/wild_mice/scripts/", "card.json")

# Parse the CARD file
Rdb = parse_card(floc, jloc=jloc)

# Save the DataFrame to a CSV file
output_csv = os.path.join("/storage/homefs/ib21q318/data/wild_mice/genomes/Bacteroides_acidifaciens/", "card_results.csv")
Rdb.to_csv(output_csv, index=False)
