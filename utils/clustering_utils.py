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
from .file_utils import remove_file
import subprocess

file_dir = os.path.dirname(os.path.abspath(__file__))

# clustering_utils


def compute_pairwise_distance_matrix(data_arr, metric, n_jobs=-1, p_minkowski=1):

    if metric == 'minkowski':
        pair_dists = pairwise_distances(data_arr, metric=metric, n_jobs=n_jobs, p=p_minkowski)
    else:
        pair_dists = pairwise_distances(data_arr, metric=metric, n_jobs=n_jobs)
    
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


def compute_umap_graph(data_df, nn, metric, nn_idxs, nn_dists, return_dists=False):
    
    result, sigmas, rhos, dists = umap.umap_.fuzzy_simplicial_set(data_df, nn, 42, metric, knn_indices=nn_idxs, knn_dists=nn_dists, return_dists=True)

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


def compute_communities(partition, idx_labels):
    communities = {}

    for idx, membership in enumerate(partition):
        if membership not in communities:
            communities[membership] = []
        communities[membership].append(idx_labels[idx])

    return communities


def compute_silhouette_score(distance_matrix, partition):
    return silhouette_score(distance_matrix, partition, metric='precomputed')


def compute_modularity(graph, communities):
    nx_g = nx.Graph(graph.get_edgelist())
    return nx.community.quality.modularity(nx_g, communities, weight='weight')


def format_partition_for_enrichment(df, partition):
    edf = pd.DataFrame.from_dict({'TTHERM_ID': []})
    edf['TTHERM_ID'] = df['TTHERM_ID'].values
    edf['label'] = partition
    return edf

import pickle

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


def run_leiden(df, n_components=2, n_neighbors=3, random_state=42, metric='manhattan', leiden_type='modularity', la_res_param=0.1, return_dists=True, p_minkowski=1):
    """
    Function to compute the simplicial sets for coexpression using UMAP and to then apply
    the Leiden algorithm to cluster the resulting graph.
    
    Parameters:
    -----------
    df : pandas dataframe
        the expression data
    n_components : int (default 2)
        the number of dimensions onto which the data should be projected
    n_neighbors : int (default 15)
        a parameter for the UMAP algorithm. I think it has to do with balancing
        local vs. global topology in the data
    random_state : float (default 42)
        Constraining this parameter makes the output reproducible
    metric : str (default "euclidean")
        The distance function
    return_dists : Bool (default True)
        Whether the function should return the computed distances
        
    Returns:
    --------
    leiden_modules : np array
        An array of ints, each corresponding to the module (or cluster) to which a gene belongs,
        listed in ortder of the input dataframe
    """ # FIXME add all params and return objects to docstring
    
    data = df[list(df.columns)[1:]].values

    labels_idxs = list(range(data.shape[0]))  

#     mapper = umap.UMAP(random_state=random_state, n_components=n_components, n_neighbors=n_neighbors).fit(data)
    
    distance_matrix = compute_pairwise_distance_matrix(data, metric=metric, n_jobs=-1, p_minkowski=p_minkowski)

    nn_idxs, nn_dists = compute_nns(data, n_neighbors, metric, random_state, distance_matrix=distance_matrix)

    g, dists = compute_umap_graph(data, n_neighbors, metric, nn_idxs, nn_dists, return_dists=True)

    if leiden_type == 'modularity':
        partition = la.find_partition(g, la.ModularityVertexPartition, seed=random_state, weights='weight')
        leiden_modules = np.array(partition.membership)
    elif leiden_type == 'cpm':
        leiden_modules =  compute_leiden_partition(g, la_res_param, random_state)
    else:
        raise ValueError('Invalid value for leiden_type parameter')
    
    sscore = compute_silhouette_score(distance_matrix, leiden_modules)
    communities = compute_communities(leiden_modules, labels_idxs)
    modularity = compute_modularity(g, communities.values())
    
    return leiden_modules, dists, sscore, modularity


def build_label_df(data_df, phases, random_state=42, n_neighbors=3, metric='manhattan', leiden_type='modularity', la_res_param=1.0, lldf=None):
    """
    Function to build a dataframe of genes labeled according to their UMAP/Leiden modules
    
    Parameters:
    -----------
    data_df : pandas DataFrame
        The expression data
    phases : str ('full', 'veg', or 'sex')
        The physiological phases for which expression data is being provided
    lldf : pandas DataFrame (default None)
        Another leiden label df (lldf) to which to add a column
        
    Returns:
    --------
    lldf : pandas DataFrame
        Leiden Label DataFrame. Gene IDs and their corresponding UMAP/Leiden module
        computed for a specific physiological regime (full set (full), vegetative only
        (veg), or sexual only (sex))
    """ # FIXME add all params and return objects to docstring
    
    if type(lldf) == type(None):
        lldf = pd.DataFrame.from_dict({'TTHERM_ID': []})
    
    leiden_modules, dists, sscore, modularity = run_leiden(data_df, random_state=random_state, n_neighbors=n_neighbors, metric=metric, leiden_type=leiden_type, la_res_param=la_res_param)
    lldf['TTHERM_ID'] = data_df['TTHERM_ID'].values
    
    lldf['label'] = leiden_modules
    
    return lldf, dists, sscore, modularity
