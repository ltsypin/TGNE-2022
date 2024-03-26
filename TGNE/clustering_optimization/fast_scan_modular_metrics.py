import tqdm
import numpy as np
import pandas as pd
from datetime import datetime
import sys
import os

file_dir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(file_dir, '../../'))
from utils import file_utils, microarray_utils, clustering_utils


def min_max_scale_2d_arr(arr: np.array):

    flat_arr = arr.flatten()

    non_inf_mask = flat_arr != float('inf')

    max_val = max(flat_arr[non_inf_mask])
    min_val = min(flat_arr)

    scaled_arr = (arr - min_val) / (max_val - min_val)

    return scaled_arr


def get_clr_dist_arr(nn: str):
    clr_df = pd.read_csv(os.path.join(file_dir, f'./clr_network_for_distances_{nn}.csv.gz'), compression='gzip')
    clr_df.rename(columns={'Unnamed: 0':'TTHERM_ID'}, inplace=True)

    zscore_arr = clr_df.loc[:,clr_df.columns[1:]].to_numpy()

    info = np.finfo(np.float64)
    smallest_float = info.eps

    scaled_zscore_arr = min_max_scale_2d_arr(zscore_arr)
    clr_dist_arr = np.sqrt(2 * (1 - scaled_zscore_arr)) + smallest_float

    np.fill_diagonal(clr_dist_arr, 0)

    return clr_dist_arr


partition_type = 'EXP'
full_filtered_df = pd.read_csv(os.path.join(file_dir, '../microarray_probe_alignment_and_filtering/allgood_filt_agg_tidy_2021aligned_qc_rma_expression_full.csv'))
full_filtered_df = full_filtered_df.rename(columns={'Unnamed: 0': 'TTHERM_ID'})

# full_filtered_df = shuffle_rows(full_filtered_df)
# partition_type = 'NC'

full_filtered_norm_df = microarray_utils.normalize_expression_per_gene(full_filtered_df)
raw_data = full_filtered_norm_df[list(full_filtered_norm_df.columns)[1:]].values

# SCAN START
curr_datetime = str(datetime.now())
idx_labels = list(range(raw_data.shape[0]))

p_minkowski = None

metrics = sys.argv[1:]

n_jobs = -1
random_state = 42


# scan_nns = np.arange(2, 13, 1)
scan_nns = [5]

# scan_rps = np.arange(0.005, 1.1, 0.005)
# scan_rps = np.arange(0.1, 1.1, 0.1)
scan_rps = [0.030, 0.035]


scan_dict = {}

for metric_p in metrics:
    metric_p_split = metric_p.split('_')

    metric = metric_p

    if metric_p_split[0] == 'minkowski':
        metric = metric_p_split[0]
        p_minkowski = float(metric_p_split[1])

    if metric != 'clr':
        distance_matrix = clustering_utils.compute_pairwise_distance_matrix(raw_data, metric, n_jobs, p_minkowski)
        nn_idxs, nn_dists = clustering_utils.compute_nns(raw_data, max(scan_nns), metric, random_state, n_jobs, p_minkowski, distance_matrix)

    for idx, nn in enumerate(scan_nns):     
        print(idx+1,'of',len(scan_nns))     
        print('NNs: ', nn)

        if metric == 'clr':
            distance_matrix = get_clr_dist_arr(int(nn))
            nn_idxs, nn_dists = clustering_utils.compute_nns(raw_data, max(scan_nns), metric, random_state, n_jobs, p_minkowski, distance_matrix)

        scan_dict[nn] = {}

        scan_dict[nn]['nn_idxs'] = nn_idxs
        scan_dict[nn]['nn_dists'] = nn_dists

        nn_graph = clustering_utils.compute_umap_graph(raw_data, nn, metric, nn_idxs, nn_dists)
        scan_dict[nn]['nn_graph'] = nn_graph

        for rp in tqdm.tqdm(scan_rps, 'RESOLUTION PARAMETERS COMPUTED'):

            scan_dict[nn][rp] = {}
            
            partition = clustering_utils.compute_leiden_partition(nn_graph, rp, random_state)
            scan_dict[nn][rp]['partition'] = partition

            communities = clustering_utils.compute_communities(partition, idx_labels)
            scan_dict[nn][rp]['communities'] = communities

            sil_score = clustering_utils.compute_silhouette_score(distance_matrix, partition)
            scan_dict[nn][rp]['sil_score'] = sil_score

            modularity = clustering_utils.compute_modularity(nn_graph, communities.values())
            scan_dict[nn][rp]['modularity'] = modularity

            enrichment_df = clustering_utils.compute_enrichment(full_filtered_norm_df, partition)
            scan_dict[nn][rp]['enrichment_df'] = enrichment_df

            num_clusters = clustering_utils.compute_num_clusters(partition, communities.values())
            scan_dict[nn][rp]['num_clusters'] = num_clusters

            num_enriched_clusters = clustering_utils.compute_num_enriched_clusters(enrichment_df)
            scan_dict[nn][rp]['num_enriched_clusters'] = num_enriched_clusters

            num_enriched_cluster_genes = clustering_utils.compute_num_enriched_cluster_genes(enrichment_df, partition)
            scan_dict[nn][rp]['num_enriched_cluster_genes'] = num_enriched_cluster_genes

            cluster_sizes = clustering_utils.compute_cluster_sizes(communities)
            scan_dict[nn][rp]['cluster_sizes'] = cluster_sizes

            enriched_cluster_sizes = clustering_utils.compute_enriched_cluster_sizes(communities, enrichment_df)
            scan_dict[nn][rp]['enriched_cluster_sizes'] = enriched_cluster_sizes

            cluster_stats = {
            'partition_type': partition_type,

            'dimensionality': 'baseline',

            'metric': metric_p,
            'graph': 'umap_fuzzy_simplicial_set',
            'nns': nn,

            'clustering': 'leiden_cpm',
            'parameter': rp,

            'silhouette_score': sil_score,
            'modularity': modularity,

            'nclusters': num_clusters,
            'mean_cluster_size': clustering_utils.compute_cluster_size_mean(cluster_sizes),
            'median_cluster_size': clustering_utils.compute_cluster_size_median(cluster_sizes),
            'sd_cluster_size': clustering_utils.compute_cluster_size_sd(cluster_sizes),

            'nenriched_clusters': num_enriched_clusters,
            'mean_enriched_cluster_size': clustering_utils.compute_cluster_size_mean(enriched_cluster_sizes),
            'median_enriched_cluster_size': clustering_utils.compute_cluster_size_median(enriched_cluster_sizes),
            'sd_enriched_cluster_size': clustering_utils.compute_cluster_size_sd(enriched_cluster_sizes),
            'nenriched_cluster_genes': num_enriched_cluster_genes,

            'datetime': curr_datetime
            }

            file_utils.write_to_csv(os.path.join(file_dir, (f'./{"_".join([m for m in metrics])}_{curr_datetime.replace(" ", "_").replace(":", "-")}_scan_stats.csv')), cluster_stats, list(cluster_stats.keys()))
