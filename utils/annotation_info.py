import os
import pandas as pd

file_dir = os.path.dirname(os.path.abspath(__file__))

main_dir = '../TGNE/'

go_df = pd.read_csv(os.path.join(file_dir, main_dir, '../active_fastas/go_annotations.csv'))
kegg_df = pd.read_csv(os.path.join(file_dir, main_dir, '../active_fastas/kegg_annotations.csv'))
ec_df = pd.read_csv(os.path.join(file_dir, main_dir, '../active_fastas/ec_annotations.csv'))


def get_GO_info(go_term):
    
    name = go_df['GOs_description'].loc[go_df['GOs'] == go_term].values[0]
    
    definition = go_df['GOs_definition'].loc[go_df['GOs'] == go_term].values[0]
    
    obsolete = go_df['GOs_obsolete'].loc[go_df['GOs'] == go_term].values[0]
    
    return name, definition, obsolete

def get_KEGG_info(term):
    return kegg_df['KEGG_ko_description'].loc[kegg_df['KEGG_ko'] == term].values[0]

def get_EC_info(term):
    return ec_df['EC_description'].loc[ec_df['EC'] == term].values[0]

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