import tqdm
import numpy as np
import pandas as pd
from datetime import datetime
import sys
import os
import warnings
from glob import glob
import pickle
import argparse

file_dir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(file_dir, '../../'))

from utils import file_utils, clustering_utils, dataframe_utils, expr_data_utils, print_utils

# SCAN START
curr_datetime = str(datetime.now())

p_minkowski = None

parser = argparse.ArgumentParser(description='')

# PARAMETERS
################################################################################

parser.add_argument('--expression_dataset', 
                    required=True, 
                    type=str, 
                    help='\'microarray\' or \'rna_seq\'')

parser.add_argument('--norm_type', 
                    required=True, 
                    type=str, 
                    help='\'z_score\' or \'min_max\'')

parser.add_argument('--metric', 
                    required=True, 
                    type=str, 
                    help='distance metric (e.g., manhattan)')

parser.add_argument('--scan_nn_start', 
                    required=True, 
                    type=int, 
                    help='number of nearest neighbors to use for building the graph (range start inclusive)')

parser.add_argument('--scan_nn_end', 
                    required=True, 
                    type=int, 
                    help='number of nearest neighbors to use for building the graph (range end exclusive)')

parser.add_argument('--scan_nn_step', 
                    required=True, 
                    type=int, 
                    help='number of nearest neighbors to use for building the graph (step size)')

parser.add_argument('--scan_rps_start', 
                    required=True, 
                    type=float, 
                    help='resolution parameter to use for Leiden network detection (range start inclusive)')

parser.add_argument('--scan_rps_end', 
                    required=True, 
                    type=float, 
                    help='resolution parameter to use for Leiden network detection (range end exclusive)')

parser.add_argument('--scan_rps_step', 
                    required=True, 
                    type=float, 
                    help='resolution parameter to use for Leiden network detection (step size)')

parser.add_argument('--partition_type', 
                    required=True, 
                    type=str, 
                    help='experimental (\'EXP\'), scrambled negtive control (\'NC\'), or simulated negative control (\'TNC\')')

parser.add_argument('--print_stats', 
                    default='n', 
                    type=str, 
                    help='yes (\'y\') or no (\'n\')')

parser.add_argument('--num_simulations', 
                    default=1, 
                    type=int, 
                    help='number of simulations to run (relevant for scrambled negtive control (\'NC\') and simulated negative control (\'TNC\'))')

args = parser.parse_args()

print_stats = args.print_stats.lower()

expression_dataset = args.expression_dataset

norm_type = args.norm_type

metrics = [args.metric]

scan_nns = np.arange(args.scan_nn_start, args.scan_nn_end, args.scan_nn_step)

scan_rps = np.arange(args.scan_rps_start, args.scan_rps_end, args.scan_rps_step)

partition_type = args.partition_type

num_iterations = args.num_simulations

################################################################################

if expression_dataset == 'microarray':
    expression_data_path = os.path.join(file_dir, '../../active_files/allgood_filt_agg_tidy_2021aligned_qc_rma_expression_full.csv')
elif expression_dataset == 'rna_seq':
    expression_data_path = os.path.join(file_dir, '../../active_files/rna_seq.csv')
else:
    raise(ValueError(f'INVALID EXPRESSION DATASET: {expression_dataset}.'))

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

    gene_ids = list(msf_df.iloc[:, 0].values)
    gene_lists[file_name] = gene_ids


################################################################################

output_file = ''

n_jobs = -1
random_state = 42

if num_iterations > 1 and partition_type not in ['NC', 'TNC']:
    raise(ValueError(f'PARTITION TYPE IS SET TO {partition_type}. {num_iterations} IDENTICAL PARTITIONS WILL BE COMPUTED. PLEASE SET NUM ITERATIONS TO 1.'))

if num_iterations > 1 and partition_type in ['NC', 'TNC'] and len(scan_rps) > 1:
    raise(ValueError(f'PARTITION TYPE IS SET TO {partition_type}. {num_iterations} PARTITIONS WILL BE COMPUTED FOR {str(len(scan_rps))} RESOLUTION PARAMETERS. PLEASE SELECT A SINGLE RESOLUTION PARAMETER.'))

if num_iterations == 1 and partition_type != 'EXP':
    raise(ValueError(f'PARTITION TYPE IS SET TO {partition_type}. ONLY {num_iterations} PARTITIONS WILL BE COMPUTED FOR {str(len(scan_rps))} RESOLUTION PARAMETERS. PLEASE SET PARTITION TYPE TO \'EXP\'.'))

if partition_type not in ['EXP', 'NC', 'TNC']:
    raise(ValueError(f"PARTITION TYPE IS SET TO {partition_type}. PLEASE SELECT A VALID PARTITION TYPE ('EXP', 'NC', or 'TNC')."))

for idx, iteration in enumerate(range(num_iterations)):
    print('COMPUTING', idx+1,'of', num_iterations, 'ITERATIONS')     

    full_filtered_df = pd.read_csv(expression_data_path)
    

    data_min = np.min([full_filtered_df[c].min() for c in (list(full_filtered_df.columns)[1:])])
    data_max = np.max([full_filtered_df[c].max() for c in (list(full_filtered_df.columns)[1:])])

    # if partition_type == 'EXP':
    #     print(data_min)
    #     print(data_max)


    if partition_type == 'NC':
        full_filtered_df = dataframe_utils.shuffle_rows(full_filtered_df)
    
    if partition_type == 'TNC':
        full_filtered_df_shape = full_filtered_df.shape

        full_filtered_df_columns = list(full_filtered_df.columns)[1:]

        full_filtered_df_ttherm_ids = full_filtered_df['TTHERM_ID'].values

        full_filtered_df = None

        raw_data = dataframe_utils.get_hypercube_sample(full_filtered_df_shape[1], full_filtered_df_shape[0])
        # print(np.min([raw_data[c].min() for c in raw_data.columns]))
        # print(np.min([raw_data[c].max() for c in raw_data.columns]))

        raw_data = dataframe_utils.scale_df_values(raw_data, data_min, data_max)
        # print(np.min([raw_data[c].min() for c in raw_data.columns]))
        # print(np.min([raw_data[c].max() for c in raw_data.columns]))


        cols_to_add = full_filtered_df_columns
        full_filtered_df = pd.DataFrame({'TTHERM_ID': full_filtered_df_ttherm_ids})

        for idx, col in enumerate(cols_to_add):
            full_filtered_df[col] = raw_data[idx].values
    
    if expression_dataset == 'microarray':
        full_filtered_norm_df = expr_data_utils.normalize_expression_per_gene(full_filtered_df, norm_type=norm_type, add_scalar=0)
        full_filtered_norm_df = expr_data_utils.mean_df_of_duplicates(full_filtered_norm_df, mean_type='arithmetic')
    elif expression_dataset == 'rna_seq':
        full_filtered_norm_df = expr_data_utils.normalize_expression_per_gene(full_filtered_df, norm_type=norm_type, add_scalar=1)
        full_filtered_norm_df = expr_data_utils.mean_df_of_duplicates(full_filtered_norm_df, mean_type='arithmetic')
    else:
        raise(ValueError(f'INVALID EXPRESSION DATASET: {expression_dataset}.'))
    
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

            nn_graph = clustering_utils.compute_umap_graph(raw_data, nn, metric, nn_idxs, nn_dists, random_state=random_state)

            for rp in tqdm.tqdm(scan_rps, 'RESOLUTION PARAMETERS COMPUTED'):
                
                partition = clustering_utils.compute_leiden_partition(nn_graph, rp, random_state=random_state)

                enrichment_df = clustering_utils.compute_enrichment(full_filtered_norm_df, partition)

                try:
                    output_dir = os.path.join(file_dir, (f'./{expression_dataset}_{partition_type}_{"_".join([m for m in metrics])}_{"_".join([str(n) for n in scan_nns])}_{curr_datetime.replace(" ", "_").replace(":", "-")}_scan_stats'))
                    file_utils.create_directories(output_dir)
                    with open(os.path.join(output_dir, f'./{expression_dataset}_{partition_type}_{metric_p}_{str(nn)}_{curr_datetime.replace(" ", "_").replace(":", "-")}.pkl'), 'wb') as f:
                        pickle.dump({
                                    'partition': partition,
                                    # 'nn_graph': nn_graph,
                                    # 'distance_matrix': distance_matrix,
                                    'enrichment_df': enrichment_df,
                                    },
                                    f)
                        
                except Exception as e:
                    output_dir = os.path.join(file_dir, (f'./{expression_dataset}_{partition_type}_{curr_datetime.replace(" ", "_").replace(":", "-")}_scan_stats'))
                    file_utils.create_directories(output_dir)
                    with open(os.path.join(output_dir, f'./{expression_dataset}_{partition_type}_{metric_p}_{str(nn)}_{curr_datetime.replace(" ", "_").replace(":", "-")}.pkl'), 'wb') as f:
                        pickle.dump({
                                    'partition': partition,
                                    # 'nn_graph': nn_graph,
                                    # 'distance_matrix': distance_matrix,
                                    'enrichment_df': enrichment_df,
                                    },
                                    f)
                        
                communities = clustering_utils.compute_communities(partition, idx_labels)

                sil_score = clustering_utils.compute_silhouette_score(distance_matrix, partition, random_state=random_state)

                modularity = clustering_utils.compute_modularity(nn_graph, communities.values())

                num_clusters = clustering_utils.compute_num_clusters(partition, communities.values())

                num_enriched_clusters = clustering_utils.compute_num_enriched_clusters(enrichment_df)

                num_enriched_cluster_genes = clustering_utils.compute_num_enriched_cluster_genes(enrichment_df, partition)

                cluster_sizes = clustering_utils.compute_cluster_sizes(communities)

                enriched_cluster_sizes = clustering_utils.compute_enriched_cluster_sizes(communities, enrichment_df)

                cluster_stats = clustering_utils.compute_parition_statistics(partition_type, metric_p, nn, rp, sil_score, modularity, num_clusters, cluster_sizes, partition, num_enriched_clusters, enriched_cluster_sizes, num_enriched_cluster_genes, curr_datetime)

                if print_stats == 'y':
                    print_utils.print_aligned_dict(cluster_stats)

                for file_name, id_list in gene_lists.items():
                    cluster_stats[f'max_fraction_same_cluster_{file_name}'] = clustering_utils.fraction_max_same_cluster_genes(
                    id_list, clustering_utils.format_partition_for_enrichment(
                        full_filtered_norm_df, partition), 
                    print_mode=False)
                    
                try:
                    output_file = os.path.join(file_dir, (f'./{expression_dataset}_{partition_type}_{"_".join([m for m in metrics])}_{"_".join([str(n) for n in scan_nns])}_{curr_datetime.replace(" ", "_").replace(":", "-")}_scan_stats.csv'))
                    file_utils.write_to_csv(output_file, cluster_stats)
                except Exception as e:
                    output_file = os.path.join(file_dir, (f'./{expression_dataset}_{partition_type}_{curr_datetime.replace(" ", "_").replace(":", "-")}_scan_stats.csv'))
                    file_utils.write_to_csv(output_file, cluster_stats)

if os.path.exists(output_file):
    print(f'SCAN RESULTS WRITTEN TO: {output_file}')
