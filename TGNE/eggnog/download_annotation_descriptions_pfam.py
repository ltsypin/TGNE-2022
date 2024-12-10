import requests
import json
import pandas as pd
from multiprocessing import Pool
import numpy as np
import os


file_dir = os.path.dirname(os.path.abspath(__file__))


def get_PFAM_info(pfam_short_name):

    try:
        r = requests.get(f'https://www.ebi.ac.uk/ebisearch/ws/rest/proteinFamilies/?query={pfam_short_name}&size=15&requestFrom=searchBox&format=JSON&fieldurl=true&viewurl=true&fields=creation_date%2Cdescription%2Cgpcr-family%2Cname&hlfields=creation_date%2Cdescription%2Cgpcr-family%2Cname&entryattrs=score')    
        json_data = r.text

        data = json.loads(json_data)

        info = None

        for entry in data['entries']:
            info = {
                'id': entry['id'],
                'source': entry['source'],
                'description': entry['fields']['description'][0] if entry['fields']['description'] else '',
                'value': entry['fieldURLs'][0]['value'] if entry['fieldURLs'] else ''
            }

            if info['id'] == pfam_short_name:
                return info['description']
    
    except:
        return 'NA'
    
    return 'NA'
        

def process_pfam_chunk(annotation_df):

    pfam_terms = annotation_df['PFAMs'].values
    pfam_descriptions = []

    # for t in tqdm(pfam_terms, 'pfam_terms'):
    for t in pfam_terms:
        none = '-'
        if t == none:
            pfam_descriptions.append(none)
            continue
        desc = get_PFAM_info(t)
        pfam_descriptions.append(desc)

    annotation_df['PFAMs_description'] = pfam_descriptions

    return annotation_df


if __name__ == '__main__':

    complete_annotation_df = pd.read_csv(os.path.join(file_dir, '../../active_files/eggnog_annotations.csv'))

    pfam_terms = []

    for terms in complete_annotation_df['PFAMs'].values:
        pfam_terms += terms.split(',')

    pfam_df = pd.DataFrame({
        'PFAMs': np.unique(pfam_terms)
    })

    print(pfam_df.shape)

    pfam_num_chunks = 30

    pfam_chunks = np.array_split(pfam_df, pfam_num_chunks)

    with Pool(processes=pfam_num_chunks) as pool:
        pfam_processed_chunks = pool.map(process_pfam_chunk, pfam_chunks)

    pfam_result_df = pd.concat(pfam_processed_chunks, ignore_index=True)

    print(pfam_result_df.shape)

    pfam_result_df.to_csv(os.path.join(file_dir, '../../active_files/pfam_annotations.csv'), index=False)
