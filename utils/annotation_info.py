import os
import pandas as pd

file_dir = os.path.dirname(os.path.abspath(__file__))

main_dir = '../TGNE/'

go_df = pd.read_csv(os.path.join(file_dir, main_dir, '../active_files/go_annotations.csv'))
kegg_df = pd.read_csv(os.path.join(file_dir, main_dir, '../active_files/kegg_annotations.csv'))
ec_df = pd.read_csv(os.path.join(file_dir, main_dir, '../active_files/ec_annotations.csv'))
pfam_df = pd.read_csv(os.path.join(file_dir, main_dir, '../active_files/pfam_annotations.csv'))
interpro_df = pd.read_csv(os.path.join(file_dir, main_dir, '../active_files/interpro_annotations.csv'))


def get_GO_info(go_term):
    try:
        name = go_df['GOs_description'].loc[go_df['GOs'] == go_term].values[0]
        definition = go_df['GOs_definition'].loc[go_df['GOs'] == go_term].values[0]
        obsolete = go_df['GOs_obsolete'].loc[go_df['GOs'] == go_term].values[0]
        
        return name, definition, obsolete
    except:
        raise ValueError(f'THE GO TERM {go_term} DOES NOT HAVE A DESCRIPTION ENTRY.')
    
def get_KEGG_info(term):
    try:
        return kegg_df['KEGG_ko_description'].loc[kegg_df['KEGG_ko'] == term].values[0]
    except:
        raise ValueError(f'THE KEGG TERM {term} DOES NOT HAVE A DESCRIPTION ENTRY.')
    
def get_EC_info(term):
    try:
        return ec_df['EC_description'].loc[ec_df['EC'] == term].values[0]
    except:
        raise ValueError(f'THE EC TERM {term} DOES NOT HAVE A DESCRIPTION ENTRY.')

def get_PFAM_info(term):
    try:
        return pfam_df['PFAMs_description'].loc[pfam_df['PFAMs'] == term].values[0]
    except:
        raise ValueError(f'THE PFAM TERM {term} DOES NOT HAVE A DESCRIPTION ENTRY.')
    
def get_InterPro_info(term):
    try:
        return interpro_df['InterPro_description'].loc[interpro_df['InterPro'] == term].values[0]
    except:
        raise ValueError(f'THE InterPro TERM {term} DOES NOT HAVE A DESCRIPTION ENTRY.')

# As of 2020 https://www.ncbi.nlm.nih.gov/research/cog/
COG_dict = {
    "A" : "RNA processing and modification",
    "B" : "Chromatin structure and dynamics",
    "C" : "Energy production and conversion",
    "D" : "Cell cycle control, cell division, chromosome partitioning",
    "E" : "Amino acid transport and metabolism",
    "F" : "Nucleotide transport and metabolism",
    "G" : "Carbohydrate transport and metabolism",
    "H" : "Coenzyme transport and metabolism",
    "I" : "Lipid transport and metabolism",
    "J" : "Translation, ribosomal structure and biogenesis",
    "K" : "Transcription",
    "L" : "Replication, recombination, and repair",
    "M" : "Cell wall/membrane/envelope biogenesis",
    "N" : "Cell motility",
    "O" : "Posttranslational modification, protein turnover, chaperones",
    "P" : "Inorganic ion transport and metabolism",
    "Q" : "Secondary metabolites biosynthesis, transport and catabolism",
    "T" : "Signal transduction mechanisms",
    "U" : "Intracellular trafficking, secretion, and vesicular transport",
    "V" : "Defense mechanisms",
    "W" : "Extracellular structures",
    "X" : "Mobilome: prophages, transposons",
    "Y" : "Nuclear structure",
    "Z" : "Cytoskeleton",
    "R" : "General function prediction only",
    "S" : "Function unknown",
}