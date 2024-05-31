import pandas as pd
import sys
import os
from multiprocessing import Pool
import tqdm

file_dir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(file_dir, '../../'))

from utils import clustering_utils

def process_community(community: tuple):

    rna_seq_mod = []
    rna_seq_overlap_fraction = []

    microarray_mod = []
    microarray_overlap_fraction = []

    r_mod = community[0]
    r_set = community[1]

    for m_mod, m_set in microarray_communities_sets.items():
        rna_seq_mod.append(r_mod)
        rna_seq_overlap_fraction.append(len(set.intersection(r_set, m_set)) / len(r_set))

        microarray_mod.append(m_mod)
        microarray_overlap_fraction.append(len(set.intersection(r_set, m_set)) / len(m_set))

    comparison_df = pd.DataFrame(
    {
    'rna_seq_mod': rna_seq_mod,
    'rna_seq_overlap_fraction': rna_seq_overlap_fraction,
    'microarray_mod': microarray_mod,
    'microarray_overlap_fraction': microarray_overlap_fraction,
    }
    )

    return comparison_df

def init_pool(data):
    global microarray_communities_sets
    microarray_communities_sets = data

if __name__ == '__main__':

    rna_seq_clusters_df = pd.read_csv(os.path.join(file_dir, './rna_seq_label_df_round_1.csv'))

    rna_seq_labels = list(rna_seq_clusters_df['TTHERM_ID'].values)

    rna_seq_parition = list(rna_seq_clusters_df['label'].values)

    rna_seq_communities = clustering_utils.compute_communities(rna_seq_parition, rna_seq_labels)

    rna_seq_communities_sets = {k : set(rna_seq_communities[k]) for k in rna_seq_communities.keys()}

    # print(rna_seq_clusters_df.head())
    # print(rna_seq_labels[0])
    # print(rna_seq_parition[0])
    # print(rna_seq_communities[max(rna_seq_communities.keys())])

    microarray_clusters_df = pd.read_csv(os.path.join(file_dir, './test_nn3_leiden_label_df_round_1.csv'))

    microarray_labels = list(microarray_clusters_df['TTHERM_ID'].values)

    microarray_parition = list(microarray_clusters_df['label'].values)

    microarray_communities = clustering_utils.compute_communities(microarray_parition, microarray_labels)

    microarray_communities_sets = {k : set(microarray_communities[k]) for k in microarray_communities.keys()}

    # print(len(set.intersection(rna_seq_communities_sets[403], microarray_communities_sets[616])) / len(microarray_communities_sets[616]))

    input_data = tuple(rna_seq_communities_sets.items())

    # Process data in parallel
    with Pool(initializer=init_pool, initargs=(microarray_communities_sets,)) as pool:
        module_dfs = list(tqdm.tqdm(pool.imap(process_community, input_data), total=len(input_data)))

    # Concatenate results
    all_results = pd.concat(module_dfs)

    all_results.to_csv(os.path.join(file_dir, './module_comparison.csv'), index=False)
