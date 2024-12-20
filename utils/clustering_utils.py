import os
import sys
import pandas as pd
import numpy as np
from sklearn.metrics import silhouette_score, pairwise_distances
from sklearn.neighbors import NearestNeighbors
# from pynndescent import NNDescent
import umap
import igraph as ig
import leidenalg as la
import networkx as nx
from datetime import datetime
import pickle
import subprocess
from collections import OrderedDict

# from .file_utils import remove_file
from .dataframe_utils import get_hypercube_sample, shuffle_rows

file_dir = os.path.dirname(os.path.abspath(__file__))
main_dir = '../TGNE/'

# clustering_utils


def min_max_scale_2d_arr(arr: np.array):

    flat_arr = arr.flatten()

    non_inf_mask = flat_arr != float('inf')

    max_val = max(flat_arr[non_inf_mask])
    min_val = min(flat_arr)

    scaled_arr = (arr - min_val) / (max_val - min_val)

    return scaled_arr

def get_clr_dist_arr(nn: str, folder: str):
    path = os.path.join(file_dir, main_dir, 'clustering_optimization', folder, f'./clr_network_for_distances_{nn}.csv.gz')
    clr_df = pd.read_csv(path, compression='gzip')
    clr_df.rename(columns={'Unnamed: 0':'TTHERM_ID'}, inplace=True)

    zscore_arr = clr_df.loc[:,clr_df.columns[1:]].to_numpy()

    print(np.sum(np.isnan(zscore_arr)))

    info = np.finfo(np.float64)
    smallest_float = info.eps

    scaled_zscore_arr = min_max_scale_2d_arr(zscore_arr)
    print(np.sum(np.isnan(scaled_zscore_arr)))
    clr_dist_arr = np.sqrt(2 * (1 - scaled_zscore_arr)) + smallest_float

    print(np.sum(np.isnan(clr_dist_arr)))

    np.fill_diagonal(clr_dist_arr, 0)

    print(np.sum(np.isnan(clr_dist_arr)))

    return clr_dist_arr

def get_clr_dist_arr_lev(nn: str, folder: str):
    path = os.path.join(file_dir, main_dir, 'clustering_optimization', folder, f'./clr_network_for_distances_{nn}.csv.gz')
    print(path)
    clr_df = pd.read_csv(path, compression='gzip')
    clr_df.rename(columns={'Unnamed: 0':'TTHERM_ID'}, inplace=True)

    zscore_arr = clr_df.loc[:,clr_df.columns[1:]].to_numpy()

    info = np.finfo(np.float64)
    smallest_float = info.eps

    # 1 / zscore for all zscores != zero
    # scale values linearly 0 to 1
    # assign 1s to idxs where original zscores == zero
    zero_zscore_idxs = np.where(zscore_arr == 0)
    inverse_zscore_arr = 1 / zscore_arr

    scaled_inverse_zscore_arr = min_max_scale_2d_arr(inverse_zscore_arr)

    for idx in range(len(zero_zscore_idxs[0])):
        scaled_inverse_zscore_arr[zero_zscore_idxs[0][idx]][zero_zscore_idxs[1][idx]] = 1
    
    scaled_inverse_zscore_arr = scaled_inverse_zscore_arr + smallest_float

    np.fill_diagonal(scaled_inverse_zscore_arr, 0) # diagonal must be zeros for dist matrix

    return scaled_inverse_zscore_arr


def compute_pairwise_distance_matrix(data_arr, metric, n_jobs=-1, p_minkowski=1):

    if metric == 'minkowski':
        pair_dists = pairwise_distances(data_arr, metric=metric, n_jobs=n_jobs, p=p_minkowski)
        return pair_dists
    
    if metric == 'angular':
        pair_dists = pairwise_distances(data_arr, metric='cosine', n_jobs=n_jobs)
        angular_pair_dists = (2 * np.arccos(1 - pair_dists)) / np.pi
        np.fill_diagonal(angular_pair_dists, 0)
        return angular_pair_dists

    pair_dists = pairwise_distances(data_arr, metric=metric, n_jobs=n_jobs)

    np.fill_diagonal(pair_dists, 0)
    
    return pair_dists


def compute_nns(data_df, nn, metric, random_state=42, n_jobs=-1, p_minkowski=1, distance_matrix=None):
    
    # if metric == 'clr':
    num_neighbors = NearestNeighbors(n_neighbors=nn-1, metric='precomputed', n_jobs=-1).fit(distance_matrix)
    nn_dists, nn_idxs = num_neighbors.kneighbors(return_distance=True)

    nn_dists_list = []
    nn_idxs_list = []

    # add the node itself to the nearest neighbors data 
    for idx in range(len(nn_dists)):
        nn_dists_list.append(np.flip(np.append(np.flip(nn_dists[idx]), 0)))
        nn_idxs_list.append(np.flip(np.append(np.flip(nn_idxs[idx]), idx)))

    return np.array(nn_idxs_list), np.array(nn_dists_list)

    # the below inplementation is natively integrated into umap but does not support precomputed distance matrices

    # n_trees = min(64, 5 + int(round((data_df.shape[0]) ** 0.5 / 20.0)))
    # n_iters = max(5, int(round(np.log2(data_df.shape[0]))))

    # if metric == 'minkowski':
    #     knn_search_index = NNDescent(
    #             data_df,
    #             n_neighbors=nn,
    #             metric=metric,
    #             metric_kwds={'p': p_minkowski},
    #             random_state=random_state,
    #             n_trees=n_trees,
    #             n_iters=n_iters,
    #             max_candidates=60,
    #             # low_memory=low_memory,
    #             n_jobs=n_jobs,
    #             verbose=False,
    #             compressed=False,
    #         )
    # else:
    #     knn_search_index = NNDescent(
    #                 data_df,
    #                 n_neighbors=nn,
    #                 metric=metric,
    #                 # metric_kwds=metric_kwds,
    #                 random_state=random_state,
    #                 n_trees=n_trees,
    #                 n_iters=n_iters,
    #                 max_candidates=60,
    #                 # low_memory=low_memory,
    #                 n_jobs=n_jobs,
    #                 verbose=False,
    #                 compressed=False,
    #             )
    # nn_idxs, nn_dists = knn_search_index.neighbor_graph

    # return nn_idxs, nn_dists


def compute_umap_graph(data_df, nn, metric, nn_idxs, nn_dists, return_dists=False, random_state=42):
    
    result, sigmas, rhos, dists = umap.umap_.fuzzy_simplicial_set(data_df, nn, random_state, metric, knn_indices=nn_idxs, knn_dists=nn_dists, return_dists=True)

    sources, targets = result.nonzero()
    edge_list = zip(sources, targets)
    weights = result.data

    g = ig.Graph(edges=edge_list, edge_attrs={'weight': weights})

    if return_dists:
        return g, dists
    
    return g


def compute_leiden_partition(graph, resolution_parameter, random_state=42):
        
        partition = la.find_partition(graph, la.CPMVertexPartition, resolution_parameter = resolution_parameter, seed=random_state, weights='weight')
        # partition = la.find_partition(g, la.ModularityVertexPartition, seed=42, weights='weight')

        leiden_modules = np.array(partition.membership)

        return leiden_modules


def compute_communities(partition, labels):
    communities = {}

    for idx, membership in enumerate(partition):
        if membership not in communities:
            communities[membership] = []
        communities[membership].append(labels[idx])

    return communities


def compute_silhouette_score(distance_matrix, partition, random_state=42):
    ss = None
    try:
        ss = silhouette_score(distance_matrix, partition, metric='precomputed', random_state=random_state)
    except:
        ss = float('NaN')
    return ss


def compute_modularity(graph, communities):
    nx_g = nx.Graph(graph.get_edgelist())
    return nx.community.quality.modularity(nx_g, communities, weight='weight')


def format_partition_for_enrichment(df, partition):
    edf = pd.DataFrame.from_dict({'TTHERM_ID': []})
    edf['TTHERM_ID'] = df['TTHERM_ID'].values
    edf['label'] = partition
    return edf


def compute_enrichment(df, partition=None):
    edf = df

    if partition is not None:
        edf = format_partition_for_enrichment(df, partition)

    data_bytes = pickle.dumps(edf)

    process = subprocess.Popen([sys.executable, 'fast_enrichment_analysis.py'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, cwd=file_dir)
    stdout, _ = process.communicate(data_bytes)

    cedf = pickle.loads(stdout)

    return cedf


def compute_num_clusters(partition, communities=None):
    if communities is None:
        return len(set(partition))
    
    if len(set(partition)) != len(communities):
        raise ValueError(f'The number of clusters/modules ({len(set(partition))}) in the partition != the number of communities ({len(communities)}).')
    
    return len(set(partition))


def compute_cluster_sizes(communities):
    return [len(community) for community in communities.values()]


def compute_enriched_cluster_sizes(communities, cedf):
    enriched_cluster_mods = set(cedf['module'].values)
    return [len(community) for mod, community in communities.items() if mod in enriched_cluster_mods]


def compute_cluster_size_mean(cluster_sizes):
    return np.mean(cluster_sizes)


def compute_cluster_size_median(cluster_sizes):
    return np.median(cluster_sizes)


def compute_cluster_size_sd(cluster_sizes):
    return np.std(cluster_sizes)


def compute_cluster_size_sd(cluster_sizes):
    return np.std(cluster_sizes)


def compute_num_enriched_clusters(cedf):
    return len(set(cedf['module'].values))


def compute_num_enriched_cluster_genes(edf, partition):
    total_num_genes = 0

    for m in set(edf['module'].values):
        num_genes = np.count_nonzero(partition == int(m))
        total_num_genes += num_genes
    
    return total_num_genes


# def run_leiden(df, n_components=2, n_neighbors=3, random_state=42, metric='manhattan', leiden_type='modularity', la_res_param=0.1, return_dists=True, p_minkowski=1):
#     """
#     Function to compute the simplicial sets for coexpression using UMAP and to then apply
#     the Leiden algorithm to cluster the resulting graph.
    
#     Parameters:
#     -----------
#     df : pandas dataframe
#         the expression data
#     n_components : int (default 2)
#         the number of dimensions onto which the data should be projected
#     n_neighbors : int (default 15)
#         a parameter for the UMAP algorithm. I think it has to do with balancing
#         local vs. global topology in the data
#     random_state : float (default 42)
#         Constraining this parameter makes the output reproducible
#     metric : str (default "euclidean")
#         The distance function
#     return_dists : Bool (default True)
#         Whether the function should return the computed distances
        
#     Returns:
#     --------
#     leiden_modules : np array
#         An array of ints, each corresponding to the module (or cluster) to which a gene belongs,
#         listed in ortder of the input dataframe
#     """ # FIXME add all params and return objects to docstring
    
#     data = df[list(df.columns)[1:]].values

#     labels_idxs = list(range(data.shape[0]))  

# #     mapper = umap.UMAP(random_state=random_state, n_components=n_components, n_neighbors=n_neighbors).fit(data)
    
#     distance_matrix = compute_pairwise_distance_matrix(data, metric=metric, n_jobs=-1, p_minkowski=p_minkowski)

#     nn_idxs, nn_dists = compute_nns(data, n_neighbors, metric, random_state, distance_matrix=distance_matrix)

#     g, dists = compute_umap_graph(data, n_neighbors, metric, nn_idxs, nn_dists, return_dists=True)

#     if leiden_type == 'modularity':
#         partition = la.find_partition(g, la.ModularityVertexPartition, seed=random_state, weights='weight')
#         leiden_modules = np.array(partition.membership)
#     elif leiden_type == 'cpm':
#         leiden_modules =  compute_leiden_partition(g, la_res_param, random_state)
#     else:
#         raise ValueError('Invalid value for leiden_type parameter')
    
#     sscore = compute_silhouette_score(distance_matrix, leiden_modules)
#     communities = compute_communities(leiden_modules, labels_idxs)
#     modularity = compute_modularity(g, communities.values())
    
#     return leiden_modules, dists, sscore, modularity


# def build_label_df(data_df, phases, random_state=42, n_neighbors=3, metric='manhattan', leiden_type='modularity', la_res_param=1.0, lldf=None):
#     """
#     Function to build a dataframe of genes labeled according to their UMAP/Leiden modules
    
#     Parameters:
#     -----------
#     data_df : pandas DataFrame
#         The expression data
#     phases : str ('full', 'veg', or 'sex')
#         The physiological phases for which expression data is being provided
#     lldf : pandas DataFrame (default None)
#         Another leiden label df (lldf) to which to add a column
        
#     Returns:
#     --------
#     lldf : pandas DataFrame
#         Leiden Label DataFrame. Gene IDs and their corresponding UMAP/Leiden module
#         computed for a specific physiological regime (full set (full), vegetative only
#         (veg), or sexual only (sex))
#     """ # FIXME add all params and return objects to docstring
    
#     if type(lldf) == type(None):
#         lldf = pd.DataFrame.from_dict({'TTHERM_ID': []})
    
#     leiden_modules, dists, sscore, modularity = run_leiden(data_df, random_state=random_state, n_neighbors=n_neighbors, metric=metric, leiden_type=leiden_type, la_res_param=la_res_param)
#     lldf['TTHERM_ID'] = data_df['TTHERM_ID'].values
    
#     lldf['label'] = leiden_modules
    
#     return lldf, dists, sscore, modularity

def compute_parition_statistics(partition_type, metric_p, nn, rp, sil_score, modularity, num_clusters, cluster_sizes, partition, num_enriched_clusters, enriched_cluster_sizes, num_enriched_cluster_genes, curr_datetime):
    return OrderedDict([
            ('partition_type', partition_type),

            ('dimensionality', 'baseline'),

            ('metric', metric_p),
            ('graph', 'umap_fuzzy_simplicial_set'),
            ('nns', nn),

            ('clustering', 'leiden_cpm'),
            ('parameter', rp),

            ('silhouette_score', sil_score),
            ('modularity', modularity),

            ('nclusters', num_clusters),
            ('mean_cluster_size', compute_cluster_size_mean(cluster_sizes)),
            ('median_cluster_size', compute_cluster_size_median(cluster_sizes)),
            ('sd_cluster_size', compute_cluster_size_sd(cluster_sizes)),
            ('q1_cluster_size', np.percentile(cluster_sizes, 25)),
            ('q3_cluster_size', np.percentile(cluster_sizes, 75)),
            ('max_cluster_size', np.max(cluster_sizes)),
            ('min_cluster_size', np.min(cluster_sizes)),
            ('ngenes', len(partition)),

            ('nenriched_clusters', num_enriched_clusters),
            ('mean_enriched_cluster_size', compute_cluster_size_mean(enriched_cluster_sizes)),
            ('median_enriched_cluster_size', compute_cluster_size_median(enriched_cluster_sizes)),
            ('sd_enriched_cluster_size', compute_cluster_size_sd(enriched_cluster_sizes)),
            ('q1_enriched_cluster_size', float('NaN') if num_enriched_clusters == 0 else np.percentile(enriched_cluster_sizes, 25)),
            ('q3_enriched_cluster_size', float('NaN') if num_enriched_clusters == 0 else np.percentile(enriched_cluster_sizes, 75)),
            ('max_enriched_cluster_size', float('NaN') if num_enriched_clusters == 0 else np.max(enriched_cluster_sizes)),
            ('min_enriched_cluster_size', float('NaN') if num_enriched_clusters == 0 else np.min(enriched_cluster_sizes)),
            ('nenriched_cluster_genes', num_enriched_cluster_genes),

            ('datetime', curr_datetime),
            ])


def build_label_df(data_df, dataset, metric='manhattan', n_neighbors=3, resolution_param=0.5, partition_type = 'EXP', n_jobs = -1, random_state=42):

    curr_datetime = str(datetime.now())

    metric_p = metric

    nn = n_neighbors

    rp = resolution_param

    partition_type = partition_type

    p_minkowski = None

    if partition_type not in ['EXP', 'NC', 'TNC']:
        raise(ValueError('partition_type must be \'EXP\' (experimental), \'NC\' (negative control), or \'TNC\' (true negative control).'))

    if partition_type == 'NC':
        data_df = shuffle_rows(data_df)

    if partition_type == 'TNC':
        raw_data = get_hypercube_sample(data_df.shape[1], data_df.shape[0])
    
    raw_data = data_df[list(data_df.columns)[1:]].values

    idx_labels = list(range(raw_data.shape[0]))

    metric_p_split = metric_p.split('_')

    metric = metric_p

    if metric_p_split[0] == 'minkowski':
        metric = metric_p_split[0]
        p_minkowski = float(metric_p_split[1])

    if metric != 'clr':
        distance_matrix = compute_pairwise_distance_matrix(raw_data, metric, n_jobs, p_minkowski)
        nn_idxs, nn_dists = compute_nns(raw_data, nn, metric, random_state, n_jobs, p_minkowski, distance_matrix)

    if metric == 'clr':
        distance_matrix = get_clr_dist_arr(int(nn))
        nn_idxs, nn_dists = compute_nns(raw_data, nn, metric, random_state, n_jobs, p_minkowski, distance_matrix)

    nn_graph = compute_umap_graph(raw_data, nn, metric, nn_idxs, nn_dists, random_state=random_state)
    
    partition = compute_leiden_partition(nn_graph, rp, random_state=random_state)

    communities = compute_communities(partition, idx_labels)

    sil_score = compute_silhouette_score(distance_matrix, partition, random_state=random_state)

    modularity = compute_modularity(nn_graph, communities.values())

    partition_df = format_partition_for_enrichment(data_df, partition)

    enrichment_df = compute_enrichment(partition_df)

    num_clusters = compute_num_clusters(partition, communities.values())

    num_enriched_clusters = compute_num_enriched_clusters(enrichment_df)

    num_enriched_cluster_genes = compute_num_enriched_cluster_genes(enrichment_df, partition)

    cluster_sizes = compute_cluster_sizes(communities)

    enriched_cluster_sizes = compute_enriched_cluster_sizes(communities, enrichment_df)

    cluster_stats = compute_parition_statistics(partition_type, metric_p, nn, rp, sil_score, modularity, num_clusters, cluster_sizes, partition, num_enriched_clusters, enriched_cluster_sizes, num_enriched_cluster_genes, curr_datetime)

    return partition_df, cluster_stats, cluster_sizes, enriched_cluster_sizes


def get_gene_module_assignments(all_gene_labels: list, gene_list: list, parition: list):
    gene_module_assignments = {}

    for gene in gene_list:
        if gene not in all_gene_labels:
            raise ValueError(f'The gene {gene} is not in the list of all gene labels.')
        gene_idx = all_gene_labels.index(gene)
        module_num = parition[gene_idx]
        if module_num not in gene_module_assignments:
            gene_module_assignments[module_num] = []
        gene_module_assignments[module_num].append(gene)

    return gene_module_assignments


# def multi_fraction_max_same_cluster_genes(gene_list_dict: dict, label_df: pd.DataFrame, fraction_threshold=0, print_mode=False):
#     for name, gene_list in gene_list_dict.items():
#         print('GENE_LIST:', name)
#         fraction_max_same_cluster_genes(gene_list, label_df, fraction_threshold=fraction_threshold, print_mode=print_mode)
#         print()
#         print()
#         print()

# def fraction_max_same_cluster_genes(gene_list: list, label_df: pd.DataFrame, fraction_threshold=0, print_mode=True):

#     df_y_to_ttherm = pd.read_csv(os.path.join(file_dir, '../new_raw_data/tgd2024/yf_ttherm_mapping_feb2024.csv'))
#     dict_ttherm_to_y = {ttherm: yf for yf, ttherm in zip(df_y_to_ttherm['yf2024'].values, df_y_to_ttherm['ttherm2021'].values)}
    
#     target_ids = []

#     for gene in gene_list:
#         if gene in dict_ttherm_to_y:
#             target_ids.append(f'{dict_ttherm_to_y[gene]}.t1')

#         # ACCOUNT FOR TRANSLATE YF TO TTHERM BEING ENABLED
#         if gene in label_df['TTHERM_ID'].values:
#             target_ids.append(gene)

#     gene_cluster_assignments = label_df.loc[label_df['TTHERM_ID'].isin(target_ids)]

#     fraction = 0
#     gene_cluster_assignments_mode_only_len = 0
#     gene_cluster_assignments_len = 0

#     if gene_cluster_assignments_len > 0:
#         cluster_assignment_mode = gene_cluster_assignments['label'].mode()[0]
#         gene_cluster_assignments_mode_only = gene_cluster_assignments.loc[gene_cluster_assignments['label'] == cluster_assignment_mode]

#         gene_cluster_assignments_mode_only_len = gene_cluster_assignments_mode_only.shape[0]
#         gene_cluster_assignments_len = gene_cluster_assignments.shape[0]

#         fraction = gene_cluster_assignments_mode_only_len/gene_cluster_assignments_len

#     print(fraction_threshold, '>', fraction)

#     if fraction_threshold > fraction:
#         return None

#     if print_mode:
#         print(gene_cluster_assignments_mode_only_len, '/', gene_cluster_assignments_len, '=',
#                 fraction
#                 )
#         print(gene_cluster_assignments_mode_only)
#         print(','.join(list(gene_cluster_assignments_mode_only['TTHERM_ID'].values)))

#     return fraction


def multi_fraction_max_same_cluster_genes(gene_list_dict: dict, label_df: pd.DataFrame, fraction_threshold=0, print_mode=False):
    for name, gene_list in gene_list_dict.items():
        if print_mode:
            print('GENE_LIST:', name)
        fraction_max_same_cluster_genes(gene_list, label_df, print_mode=print_mode)
        print()
        print()
        print()

def fraction_max_same_cluster_genes(gene_list: list, label_df: pd.DataFrame, fraction_threshold=0, print_mode=True):

    df_y_to_ttherm = pd.read_csv(os.path.join(file_dir, '../new_raw_data/tgd2024/yf_ttherm_mapping_feb2024.csv'))
    dict_ttherm_to_y = {ttherm: yf for yf, ttherm in zip(df_y_to_ttherm['yf2024'].values, df_y_to_ttherm['ttherm2021'].values)}
    
    target_ids = []

    for gene in gene_list:
        if gene in dict_ttherm_to_y:
            target_ids.append(f'{dict_ttherm_to_y[gene]}.t1')

        # ACCOUNT FOR TRANSLATE YF TO TTHERM BEING ENABLED
        elif gene in label_df['TTHERM_ID'].values:
            target_ids.append(gene)

        else:
            pass
            # print(gene)

    gene_cluster_assignments = label_df.loc[label_df['TTHERM_ID'].isin(target_ids)]

    fraction = 0

    if gene_cluster_assignments.shape[0] > 0:
        cluster_assignment_mode = gene_cluster_assignments['label'].mode()[0]
        gene_cluster_assignments_mode_only = gene_cluster_assignments.loc[gene_cluster_assignments['label'] == cluster_assignment_mode]

        fraction = gene_cluster_assignments_mode_only.shape[0]/gene_cluster_assignments.shape[0]

        if print_mode:
            print(gene_cluster_assignments_mode_only.shape[0], '/', gene_cluster_assignments.shape[0], '=',
                    fraction
                    )
            print(gene_cluster_assignments_mode_only)
            print(','.join(list(gene_cluster_assignments_mode_only['TTHERM_ID'].values)))

    else:
        if print_mode:
            print(gene_cluster_assignments)

    return fraction


def ari_mean_nexpr_per_mod(full_filtered_norm_df: pd.DataFrame, leiden_label_df_round_1_arranged_sorted: pd.DataFrame):    
    avg_df = None

    for m in leiden_label_df_round_1_arranged_sorted['label'].unique():

        curr_df = (full_filtered_norm_df.loc[full_filtered_norm_df['TTHERM_ID'].isin(
                        (leiden_label_df_round_1_arranged_sorted.loc[leiden_label_df_round_1_arranged_sorted['label'] == m]['TTHERM_ID'].values)
                    )].iloc[:, 1:].mean()).to_frame().T
        curr_df['label'] = m

        if avg_df is None:
            avg_df = curr_df
            continue

        avg_df = pd.concat((avg_df, curr_df), ignore_index=True)

    avg_df = avg_df.loc[: , list(avg_df.columns)[avg_df.shape[1] - 1:] + list(avg_df.columns)[0: avg_df.shape[1] - 1]]

    return avg_df

def compute_fraction_clusters_enriched(row):
    return row['nenriched_clusters'] / (row['nclusters'])