import requests
import re
import pandas as pd
import time
import numpy as np
from bs4 import BeautifulSoup


def get_KEGG_info(term):
    try:
        response = requests.get(f'https://rest.kegg.jp/get/{term}')

        name_string = re.search(r'^\s*NAME\s+(.+)$', response.text, re.MULTILINE).group()
        modified_string = re.sub(r'NAME\s+', '', name_string)
        name = re.sub(r'\s*\[EC:[^\]]*\]', '', modified_string)

    except:
        return "NA"

    return name


def get_EC_info(term):
    
    try:
        r = requests.get(f"https://enzyme.expasy.org/EC/{term}")
        html = r.text

        soup = BeautifulSoup(html, 'html.parser')

        enzyme_name_element = soup.find('strong', {'property': 'schema:name'})
        name = enzyme_name_element.get_text(strip=True)
    except:
        name = 'NA'
    
    return name


if __name__ == '__main__':

    complete_annotation_df = pd.read_csv('~/git/TGNE-2022/TGNE/eggnog/complete_eggnog_annotation.csv')

    kegg_terms = []

    for terms in complete_annotation_df['KEGG_ko'].values:
        kegg_terms += terms.split(',')

    kegg_terms = np.unique(kegg_terms)

    ec_terms = []

    for terms in complete_annotation_df['EC'].values:
        ec_terms += terms.split(',')

    ec_terms = np.unique(ec_terms)

    print(len(kegg_terms), len(ec_terms))

kegg_descs = []
ec_descs = []

k_idx = 0
e_idx = 0

none = '-'

while k_idx < len(kegg_terms) or e_idx < len(ec_terms):

    if k_idx < len(kegg_terms):
        if kegg_terms[k_idx] == none:
            kegg_descs.append(none)
        else:
            kegg_desc = get_KEGG_info(kegg_terms[k_idx])
            print(kegg_desc)
            kegg_descs.append(kegg_desc)
        k_idx += 1

    print('KEGG:', len(kegg_descs), '/', len(kegg_terms))

    if e_idx < len(ec_terms):
        if ec_terms[e_idx] == none:
            ec_descs.append(none)
        else:
            ec_desc = get_EC_info(ec_terms[e_idx])
            print(ec_desc)
            ec_descs.append(ec_desc)
        e_idx += 1

    print('EC:  ', len(ec_descs), '/', len(ec_terms))
    
    
kegg_df = pd.DataFrame({
    'KEGG_ko' : kegg_terms,
    'KEGG_ko_description' : kegg_descs,
})

ec_df = pd.DataFrame({
    'EC' : ec_terms,
    'EC_description' : ec_descs,
})

print(kegg_df.shape)

kegg_df.to_csv('~/git/TGNE-2022/TGNE/enrichment/new_kegg_annotations.csv', index=False)

print(ec_df.shape)

ec_df.to_csv('~/git/TGNE-2022/TGNE/enrichment/new_ec_annotations.csv', index=False)
