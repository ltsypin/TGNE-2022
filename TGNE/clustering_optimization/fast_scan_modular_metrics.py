import tqdm
import numpy as np
import pandas as pd
from datetime import datetime
import sys
import os
import warnings
from glob import glob

file_dir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(file_dir, '../../'))
from utils import file_utils, microarray_utils, clustering_utils, dataframe_utils

# SCAN START
curr_datetime = str(datetime.now())

p_minkowski = None

# PARAMETERS
################################################################################

expression_dataset = sys.argv[1]

if expression_dataset == 'microarray':
    expression_data_path = os.path.join(file_dir, '../microarray_probe_alignment_and_filtering/allgood_filt_agg_tidy_2021aligned_qc_rma_expression_full.csv')
elif expression_dataset == 'rna_seq':
    expression_data_path = os.path.join(file_dir, '../../active_fastas/rna_seq.csv')
else:
    raise(ValueError(f'{expression_dataset} is an invalid choice for expression_dataset.'))


# # manually curated metrics + metrics refered to in the documentation
# all_doc_metrics = ['angular', 'clr'] + ['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan'] + ['nan_euclidean'] + ['braycurtis', 'canberra', 'chebyshev', 'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']
# # manually curated metrics + metrics refered to in metric parameter ValueError (sklearn documentation is likely not updated)
# all_metrics = ['angular', 'clr'] + ['euclidean', 'l2', 'l1', 'manhattan', 'cityblock', 'braycurtis', 'canberra', 'chebyshev', 'correlation', 'cosine', 'dice', 'hamming', 'jaccard', 'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule', 'wminkowski', 'nan_euclidean', 'haversine']
# # metrics that covert data to boolean (essentially destroying all information for our data) (20 clusters are always produced with the same exact size)
# boolean_metrics = [
# 'dice',
# 'jaccard',
# 'rogerstanimoto',
# 'russellrao',
# 'sokalmichener',
# 'sokalsneath',
# 'yule',
# ]
# metrics = ['clr', 'manhattan', 'euclidean', 'cosine'] + [f'minkowski_{str(p)}' for p in np.array([0.5, 1, 2, 3, 4, 5])]
# metrics = [m for m in all_metrics if m not in metrics + boolean_metrics and m[: len('minkowski')] != 'minkowski']
# metrics = ['manhattan']
# metrics = ['clr_lev']
metrics = [sys.argv[2]]




# scan_nns = np.arange(3, 6, 1)
# scan_nns = [3]

scan_nns = [int(sys.argv[3])]

# scan_nns = np.arange(2, 13, 1)


# scan_rps = np.arange(0.05, 0.5, 0.05)
# scan_rps = [0.030, 0.035]
# scan_rps = [0.030]
# scan_rps = np.arange(0.005, 0.07, 0.005)

scan_rps = np.arange(0.005, 1.1, 0.005)


partition_type = 'EXP'
# partition_type = 'NC'
# partition_type = 'TNC'

num_iterations = 1
# num_iterations = 10

gene_lists = {}

module_subset_files = glob(os.path.join(file_dir, './modules_of_interest/*'))

for msf in module_subset_files:

    # gene_ids = []

    file_name = os.path.basename(msf)

    msf_df = pd.read_csv(msf, comment='#', sep='\t')

    if file_name in gene_lists:
        name_num = 0
        new_file_name = file_name
        while new_file_name in gene_lists:
            name_num += 1
            new_file_name = f'{file_name}_{name_num}'
        print(f'WARNING: DUPLICATE GENE LIST NAME: \'{file_name}\'. RENAMING TO \'{new_file_name}\'.')
        file_name = new_file_name

    # for id in msf_df.loc[:, 0].values:
    #     if id != 'TTHERM_ID' and (id[:len('YF_')] == 'YF_' or id[:len('TTHERM_')] =='TTHERM_'):
    #         gene_ids.append(id)
    gene_ids = list(msf_df.iloc[:, 0].values)
    gene_lists[file_name] = gene_ids


################################################################################
output_file = None

n_jobs = -1
random_state = 42

scan_dict = {}

if num_iterations > 1 and partition_type != 'NC':
    raise(ValueError(f'PARTITION TYPE IS SET TO {partition_type}. {num_iterations} IDENTICAL PARTITIONS WILL BE COMPUTED. PLEASE SET NUM ITERATIONS TO 1.'))

for iteration in tqdm.tqdm(range(num_iterations), 'ITERATIONS COMPUTED'):

    full_filtered_df = pd.read_csv(expression_data_path)
    full_filtered_df = full_filtered_df.rename(columns={'Unnamed: 0': 'TTHERM_ID'})

    if partition_type == 'NC':
        full_filtered_df = dataframe_utils.shuffle_rows(full_filtered_df)

    if expression_dataset == 'microarray':
        full_filtered_norm_df = microarray_utils.normalize_expression_per_gene(full_filtered_df, z=True)
    elif expression_dataset == 'rna_seq':
        full_filtered_norm_df = full_filtered_df
    else:
        raise(ValueError(f'INVALID EXPRESSION DATASET: {expression_dataset}.'))

    if partition_type == 'TNC':
        raw_data = dataframe_utils.get_hypercube_sample(full_filtered_df.shape[1], full_filtered_df.shape[0])
    
    raw_data = full_filtered_norm_df[list(full_filtered_norm_df.columns)[1:]].values

    idx_labels = list(range(raw_data.shape[0]))

    for metric_p in metrics:

        metric_p_split = metric_p.split('_')

        metric = metric_p

        if metric_p_split[0] == 'minkowski':
            metric = metric_p_split[0]
            p_minkowski = float(metric_p_split[1])

        if metric not in ['clr', 'clr_lev']:
            try:
                distance_matrix = clustering_utils.compute_pairwise_distance_matrix(raw_data, metric, n_jobs, p_minkowski)
            except ValueError as e:
                warnings.warn(f'The distance metric {metric_p} resulted in the following error:\n{e}')
                continue

            nn_idxs, nn_dists = clustering_utils.compute_nns(raw_data, max(scan_nns), metric, random_state, n_jobs, p_minkowski, distance_matrix)
        
        for idx, nn in enumerate(scan_nns):     
            print('COMPUTING', idx+1,'of',len(scan_nns), 'NEAREST NEIGHBORS')     
            print('NNs: ', nn)

            if expression_dataset == 'rna_seq':
                clr_networks_folder = 'rna_seq_clr_networks'

            elif expression_dataset == 'microarray':
                clr_networks_folder = 'microarray_clr_networks'

            if metric == 'clr':
                distance_matrix = clustering_utils.get_clr_dist_arr(int(nn), clr_networks_folder)
                nn_idxs, nn_dists = clustering_utils.compute_nns(raw_data, max(scan_nns), metric, random_state, n_jobs, p_minkowski, distance_matrix)

            if metric == 'clr_lev':
                distance_matrix = clustering_utils.get_clr_dist_arr_lev(int(nn), clr_networks_folder)
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
                'max_cluster_size': np.max(cluster_sizes),
                'min_cluster_size': np.min(cluster_sizes),
                'ngenes': len(partition),

                'nenriched_clusters': num_enriched_clusters,
                'mean_enriched_cluster_size': clustering_utils.compute_cluster_size_mean(enriched_cluster_sizes),
                'median_enriched_cluster_size': clustering_utils.compute_cluster_size_median(enriched_cluster_sizes),
                'sd_enriched_cluster_size': clustering_utils.compute_cluster_size_sd(enriched_cluster_sizes),
                'max_enriched_cluster_size': np.max(enriched_cluster_sizes),
                'min_enriched_cluster_size': np.min(enriched_cluster_sizes),
                'nenriched_cluster_genes': num_enriched_cluster_genes,

                'datetime': curr_datetime
                }

                for file_name, id_list in gene_lists.items():
                    print(file_name) # FIXME

                    cluster_stats[f'max_fraction_same_cluster_{file_name}'] = clustering_utils.fraction_max_same_cluster_genes(
                    id_list, clustering_utils.format_partition_for_enrichment(
                        full_filtered_norm_df, partition), 
                    print_mode=False)

                    print('SCAN: ', len(id_list)) #FIXME

                    print()
                    print()

                try:
                    output_file = os.path.join(file_dir, (f'./{expression_dataset}_{partition_type}_{"_".join([m for m in metrics])}_{"_".join([str(n) for n in scan_nns])}_{curr_datetime.replace(" ", "_").replace(":", "-")}_scan_stats.csv'))
                    file_utils.write_to_csv(output_file, cluster_stats, list(cluster_stats.keys()))
                except Exception as e:
                    output_file = os.path.join(file_dir, (f'./{expression_dataset}_{partition_type}_{curr_datetime.replace(" ", "_").replace(":", "-")}_scan_stats.csv'))
                    file_utils.write_to_csv(output_file, cluster_stats, list(cluster_stats.keys()))

if os.path.exists(output_file):
    print(f'SCAN RESULTS WRITTEN TO: {output_file}')
