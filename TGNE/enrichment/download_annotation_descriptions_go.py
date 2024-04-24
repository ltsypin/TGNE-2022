import requests
import json
import pandas as pd
from multiprocessing import Pool
import numpy as np
import os

file_dir = os.path.dirname(os.path.abspath(__file__))

def get_GO_info(go_term):

    try:
    
        r = requests.get(f'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_term}/complete', 'html5lib')
        
        go_info = json.loads(r.text)
        
        name = go_info['results'][0]['name']
        
        definition = go_info['results'][0]['definition']['text']
        
        obsolete = go_info['results'][0]['isObsolete']

    except:
        name = 'NA'
        definition = 'NA'
        obsolete = 'NA'
    
    return name, definition, obsolete


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
                break

        return info['description']
    
    except:
        return 'NA'
        


def process_go_chunk(annotation_df):

    go_terms = annotation_df['GOs'].values
    go_descriptions = []
    go_definitions = []
    go_obsoletes = []

    # for t in tqdm(go_terms, 'go_terms'):
    for t in go_terms:
        none = '-'
        if t == none:
            go_descriptions.append(none)
            go_definitions.append(none)
            go_obsoletes.append(none)
            continue
        desc, defi, obso = get_GO_info(t)
        go_descriptions.append(desc)
        go_definitions.append(defi)
        go_obsoletes.append(obso)

    annotation_df['GOs_description'] = go_descriptions
    annotation_df['GOs_definition'] = go_definitions
    annotation_df['GOs_obsolete'] = go_obsoletes

    return annotation_df


if __name__ == '__main__':

    complete_annotation_df = pd.read_csv(os.path.join(file_dir, '../../active_fastas/annotations.csv'))

    go_terms = []

    for terms in complete_annotation_df['GOs'].values:
        go_terms += terms.split(',')

    go_df = pd.DataFrame({
        'GOs': np.unique(go_terms)
    })

    go_num_chunks = 10

    go_chunks = np.array_split(go_df, go_num_chunks)

    with Pool(processes=go_num_chunks) as pool:
        go_processed_chunks = pool.map(process_go_chunk, go_chunks)

    go_result_df = pd.concat(go_processed_chunks, ignore_index=True)

    print(go_result_df.shape)

    go_result_df.to_csv(os.path.join(file_dir, './go_annotations.csv'), index=False)
