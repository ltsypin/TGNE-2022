from sklearn.manifold import MDS

import pandas as pd
import numpy as np

import sys
import os

file_dir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(file_dir, '../../'))

from utils import clustering_utils, rna_seq_utils, microarray_utils

def mds_embedding(full_filtered_norm_df: pd.DataFrame, file_name: str):
    data = full_filtered_norm_df[list(full_filtered_norm_df.columns)[1:]].values

    data_dists = clustering_utils.compute_pairwise_distance_matrix(data, 'manhattan')

    mds_mapper = MDS(n_components=2, normalized_stress='auto', dissimilarity='precomputed', n_jobs=-1)
    embedding = mds_mapper.fit_transform(data_dists)

    mds_df = pd.DataFrame(np.array(embedding), columns=('x', 'y'))

    mds_df.to_csv(os.path.join(file_dir, f'./{file_name}'), index=False)

def nmds_embedding(full_filtered_norm_df: pd.DataFrame, file_name: str):
    data = full_filtered_norm_df[list(full_filtered_norm_df.columns)[1:]].values

    data_dists = clustering_utils.compute_pairwise_distance_matrix(data, 'manhattan')

    mds_mapper = MDS(n_components=2, normalized_stress='auto', dissimilarity='precomputed', n_jobs=-1, metric=False)
    embedding = mds_mapper.fit_transform(data_dists)

    mds_df = pd.DataFrame(np.array(embedding), columns=('x', 'y'))

    mds_df.to_csv(os.path.join(file_dir, f'./{file_name}'), index=False)


full_filtered_df = pd.read_csv(os.path.join(file_dir, '../../active_fastas/rna_seq.csv'))
full_filtered_norm_df = rna_seq_utils.normalize_expression_per_gene(full_filtered_df)

full_filtered_norm_df = full_filtered_norm_df.sample(500)

mds_embedding(full_filtered_norm_df, 'rna_seq_mds.csv')

nmds_embedding(full_filtered_norm_df, 'rna_seq_nmds.csv')


full_filtered_df = pd.read_csv(os.path.join(file_dir, '../microarray_probe_alignment_and_filtering/allgood_filt_agg_tidy_2021aligned_qc_rma_expression_full.csv'))
full_filtered_df = full_filtered_df.rename(columns={'Unnamed: 0': 'TTHERM_ID'})
full_filtered_norm_df = microarray_utils.normalize_expression_per_gene(full_filtered_df, z=True)

full_filtered_norm_df = full_filtered_norm_df.sample(500)

mds_embedding(full_filtered_norm_df, 'microarray_mds.csv')

nmds_embedding(full_filtered_norm_df, 'microarray_nmds.csv')