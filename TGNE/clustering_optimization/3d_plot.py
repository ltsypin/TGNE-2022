import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import distinctipy
import glob
import os

file_dir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(file_dir, '../../'))

from utils import dataframe_utils, clustering_utils
from mpl_toolkits.mplot3d import Axes3D

##########################

# norm_type = 'min_max'
norm_type = 'z_score'

dataset = 'microarray'
# dataset = 'rna_seq'

norm_types = ['min_max', 'z_score']
datasets = ['microarray', 'rna_seq']
distance_metrics = ['manhattan', 'euclidean', 'clr', 'angular', 'correlation']

##########################

def smallest_unit(number):
    # Convert the number to a string to handle decimal places
    num_str = str(number)
    
    # Find the position of the decimal point
    decimal_pos = num_str.find('.')
    
    # If there's no decimal point, return 1 (for integers)
    if decimal_pos == -1:
        return 1
    
    # Calculate the length of the fractional part
    fractional_length = len(num_str) - decimal_pos - 1
    
    # Calculate the smallest unit
    smallest_unit = 10 ** (-fractional_length)
    
    return smallest_unit

def rgb_to_hex(rgb):
    """
    Convert RGB tuple to hexadecimal color code.
    """
    return '#{:02x}{:02x}{:02x}'.format(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))

# https://en.wikipedia.org/wiki/Pareto_front
# https://en.wikipedia.org/wiki/Maxima_of_a_point_set

def maxima_of_a_point_set(points):
    sorted_points = sorted(points, key=lambda p: p[0])

    max_y = float('-inf')
    maximal_points = []

    for point in reversed(sorted_points):
        x, y, idx = point
        if y > max_y:
            maximal_points.append(point)
            max_y = y

    return maximal_points


def pareto_optimal_points(points):
    # Sort points by the first dimension
    sorted_points = sorted(points, key=lambda p: p[0])
    
    maximal_points = []
    max_y = float('-inf')
    max_z = float('-inf')
    
    for point in reversed(sorted_points):
        x, y, z, idx = point
        if y > max_y or z > max_z:
            maximal_points.append(point)
            max_y = max(max_y, y)
            max_z = max(max_z, z)
    
    return maximal_points

def plotter(dataset, norm_type, distance_metric):

    microarray_data_pattern = os.path.join(file_dir, f'./microarray_scan/{norm_type}/exp/*.csv')
    
    rna_seq_data_pattern = os.path.join(file_dir, f'./rna_seq_scan/{norm_type}/exp/*.csv')
    
    microarray_stats_files = glob.glob(microarray_data_pattern)
    
    rna_seq_stats_files = glob.glob(rna_seq_data_pattern)
    
    microarray_stats_df = dataframe_utils.csv_files_to_df(microarray_stats_files)
    microarray_stats_df.drop_duplicates(inplace=True, ignore_index=True)
    
    rna_seq_stats_df = dataframe_utils.csv_files_to_df(rna_seq_stats_files)
    rna_seq_stats_df.drop_duplicates(inplace=True, ignore_index=True)
    
    microarray_stats_df['fraction_clusters_enriched'] = microarray_stats_df.apply(clustering_utils.compute_fraction_clusters_enriched, axis=1)
    
    rna_seq_stats_df['fraction_clusters_enriched'] = rna_seq_stats_df.apply(clustering_utils.compute_fraction_clusters_enriched, axis=1)
    
    rna_seq_stats_df['iqr_cluster_size'] = rna_seq_stats_df['q3_cluster_size'] - rna_seq_stats_df['q1_cluster_size']
    rna_seq_stats_df['iqr_enriched_cluster_size'] = rna_seq_stats_df['q3_enriched_cluster_size'] - rna_seq_stats_df['q1_enriched_cluster_size']
    rna_seq_stats_df['cv_cluster_size'] = rna_seq_stats_df['mean_cluster_size'] / rna_seq_stats_df['sd_cluster_size']
    rna_seq_stats_df['cv_enriched_cluster_size'] = rna_seq_stats_df['mean_enriched_cluster_size'] / rna_seq_stats_df['sd_enriched_cluster_size']
    microarray_stats_df['iqr_cluster_size'] = microarray_stats_df['q3_cluster_size'] - microarray_stats_df['q1_cluster_size']
    microarray_stats_df['iqr_enriched_cluster_size'] = microarray_stats_df['q3_enriched_cluster_size'] - microarray_stats_df['q1_enriched_cluster_size']
    microarray_stats_df['cv_cluster_size'] = microarray_stats_df['mean_cluster_size'] / microarray_stats_df['sd_cluster_size']
    microarray_stats_df['cv_enriched_cluster_size'] = microarray_stats_df['mean_enriched_cluster_size'] / microarray_stats_df['sd_enriched_cluster_size']
    
    
    
    # Enable interactive mode
    plt.ion()
    
    x_all_stat = 'modularity'
    
    y_all_stat = 'fraction_clusters_enriched'
    
    # y_all_stat = 'nenriched_clusters'
    
    z_all_stat = 'iqr_cluster_size'
    
    ###############
    
    if dataset == 'microarray':
        df = microarray_stats_df.loc[microarray_stats_df['metric'] == distance_metric]
    
    if dataset == 'rna_seq':
        df = rna_seq_stats_df.loc[rna_seq_stats_df['metric'] == distance_metric]
    
    df = df.reset_index()
    
    ###############
    
    # dfs = [microarray_stats_df.loc[(microarray_stats_df['parameter'] > 0) & (microarray_stats_df['parameter'] < 1)].reset_index(), 
    #        rna_seq_stats_df.loc[(rna_seq_stats_df['parameter'] > 0) & (rna_seq_stats_df['parameter'] < 1)].reset_index()]
    
    
    x_all = df[x_all_stat].values
    
    y_all = df[y_all_stat].values
    
    z_all = df[z_all_stat].values
    
    # df_idxs = df.index.tolist()
    
    # points = zip(x_all, y_all, z_all, df_idxs)
    
    # max_points = np.array(pareto_optimal_points(points))
    # max_point_idxs = [int(pt[3]) for pt in max_points]
    
    # other_idxs = list(set(df_idxs) - set(max_point_idxs))
    
    global_s = 30
    
    
    fig = plt.figure()
    
    ax = fig.add_subplot(111, projection='3d')
    
    for nn in df['nns'].unique():
        # print(nn)
        mini_df = df.loc[df['nns'] == nn]
        mini_x_all = mini_df[x_all_stat].values
        
        mini_y_all = mini_df[y_all_stat].values
        
        mini_z_all = mini_df[z_all_stat].values
    
        mini_df = mini_df.reset_index()
        
        mini_df_idxs = mini_df.index.tolist()
        
        mini_points = zip(mini_x_all, mini_y_all, mini_z_all, mini_df_idxs)
        
        max_points = np.array(pareto_optimal_points(mini_points))
        max_point_idxs = [int(pt[3]) for pt in max_points]
        
        other_idxs = list(set(mini_df_idxs) - set(max_point_idxs))
        
        ax.plot(
            mini_x_all[mini_df_idxs[1:]], 
            mini_y_all[mini_df_idxs[1:]], 
            mini_z_all[mini_df_idxs[1:]], 
            color='black', 
            label='filtered out', 
            zorder=1
        )
    
        ax.scatter(
            mini_x_all[max_point_idxs], 
            mini_y_all[max_point_idxs], 
            mini_z_all[max_point_idxs], 
            color='red', 
            marker='*', 
            s=global_s, 
            label='Pareto efficient', 
            zorder=1
        )
    
    
    df_idxs = df.index.tolist()
    
    points = zip(x_all, y_all, z_all, df_idxs)
    
    max_points = np.array(pareto_optimal_points(points))
    max_point_idxs = [int(pt[3]) for pt in max_points]
    
    max_modularity_idx = np.argmax(x_all[max_point_idxs])
    
    optimal_idx = None
    
    for idx, val in enumerate(z_all[max_point_idxs]):
        if val > 10:
            optimal_idx = idx
            break
    
    center_x = x_all[max_point_idxs][optimal_idx]
    center_y = y_all[max_point_idxs][optimal_idx]
    center_z = z_all[max_point_idxs][optimal_idx]
    ax.scatter([center_x], [center_y], [center_z], facecolors='none', edgecolor='#47EA00', linewidth=2, s=200, zorder=1,
                    # label=f'({round(center_x, 2)}, {round(center_y, 2)}, {round(center_z, 2)})'
                )
    
    label = f'({round(center_x, 2)}, {round(center_y, 2)}, {round(center_z, 2)})'
    ax.text(center_x + 0.04, center_y + 0.04, center_z + 4, label, color='black', fontsize=8, ha='left')
    
    
    # plt.annotate('', xy=(center_x, center_y), xytext=(center_x, center_y),
    #              arrowprops=dict(facecolor='#47EA00', 
    #                              arrowstyle="->",
    #                             #  shrink=0.05
    #                              ),
    #             #  rotation=45
    #              )
    
    x_extreme = x_all.max() if x_all.max() > abs(x_all.min()) else x_all.min()
    y_exterme = y_all.max() if y_all.max() > abs(y_all.min()) else y_all.min()
    
    #   ax.annotate('optimal', xy=(center_x + x_extreme*0.02, center_y + y_exterme/100), xytext=(center_x + x_extreme*0.1, center_y + y_exterme/20), zorder=1,
    #               arrowprops=dict(
    #                  #  color='#47EA00', 
    #                  arrowstyle="->", 
    #                  relpos=(1, 1), 
    #                  lw=1.5,
    #                  ))
    
    #   ax.xlabel('Modularity')
    #   ax.ylabel('# of Enriched Clusters' if y_all_stat == 'nenriched_clusters' else 'Fraction of All Clusters Enriched')
    #   ax.legend()
    
    #   for nns_value, group in df.groupby('nns'):
    #      group_df = group.copy(deep=True).sort_values(by='modularity')
    #      ax.plot(group_df['modularity'], group_df[y_all_stat], color='grey', linestyle='-', linewidth=1, zorder=-1)
    
    ax.set_xlabel('x='+x_all_stat)
    ax.set_ylabel('y='+y_all_stat)
    ax.set_zlabel('z='+z_all_stat)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 80)
    
    plt.savefig(f'./plots/3d_{dataset}_{norm_type}_{distance_metric}.svg', format='svg')
    
    # plt.show(block=True)
    
    
    # fig = plt.figure()
    
    # ax = fig.add_subplot(111, projection='3d')
    
    # ax.scatter(x_all[other_idxs], y_all[other_idxs], z_all[other_idxs], color='black', marker='x', s=global_s, label='filtered out', zorder=1)
    
    # ax.scatter(x_all[max_point_idxs], y_all[max_point_idxs], z_all[max_point_idxs], color='red', marker='*', s=global_s, label='Pareto efficient', zorder=1)
    
    # max_modularity_idx = np.argmax(x_all[max_point_idxs])
    
    # optimal_idx = None
    
    # for idx, val in enumerate(z_all[max_point_idxs]):
    #     if val > 10:
    #         optimal_idx = idx
    #         break
    
    # center_x = x_all[max_point_idxs][optimal_idx]
    # center_y = y_all[max_point_idxs][optimal_idx]
    # center_z = z_all[max_point_idxs][optimal_idx]
    # ax.scatter([center_x], [center_y], [center_z], facecolors='none', edgecolor='#47EA00', linewidth=2, s=200, zorder=1,
    #                 # label=f'({round(center_x, 2)}, {round(center_y, 2)}, {round(center_z, 2)})'
    #             )
    
    # label = f'({round(center_x, 2)}, {round(center_y, 2)}, {round(center_z, 2)})'
    # ax.text(center_x + 0.04, center_y + 0.04, center_z + 4, label, color='black', fontsize=8, ha='left')
    
    
    # # plt.annotate('', xy=(center_x, center_y), xytext=(center_x, center_y),
    # #              arrowprops=dict(facecolor='#47EA00', 
    # #                              arrowstyle="->",
    # #                             #  shrink=0.05
    # #                              ),
    # #             #  rotation=45
    # #              )
    
    # x_extreme = x_all.max() if x_all.max() > abs(x_all.min()) else x_all.min()
    # y_exterme = y_all.max() if y_all.max() > abs(y_all.min()) else y_all.min()
    
    # #   ax.annotate('optimal', xy=(center_x + x_extreme*0.02, center_y + y_exterme/100), xytext=(center_x + x_extreme*0.1, center_y + y_exterme/20), zorder=1,
    # #               arrowprops=dict(
    # #                  #  color='#47EA00', 
    # #                  arrowstyle="->", 
    # #                  relpos=(1, 1), 
    # #                  lw=1.5,
    # #                  ))
    
    # #   ax.xlabel('Modularity')
    # #   ax.ylabel('# of Enriched Clusters' if y_all_stat == 'nenriched_clusters' else 'Fraction of All Clusters Enriched')
    # #   ax.legend()
    
    # #   for nns_value, group in df.groupby('nns'):
    # #      group_df = group.copy(deep=True).sort_values(by='modularity')
    # #      ax.plot(group_df['modularity'], group_df[y_all_stat], color='grey', linestyle='-', linewidth=1, zorder=-1)
    
    # ax.set_xlabel('x='+x_all_stat)
    # ax.set_ylabel('y='+y_all_stat)
    # ax.set_zlabel('z='+z_all_stat)
    
    # plt.savefig(f'3d_{dataset}_{norm_type}.svg', format='svg')
    
    # plt.show(block=True)

for d in datasets:
    for n in norm_types:
        for m in distance_metrics:
            plotter(d, n, m)