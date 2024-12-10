import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import distinctipy
import glob
import os
from collections import OrderedDict

file_dir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(file_dir, '../../'))

from utils import dataframe_utils, clustering_utils

##########################

# norm_type = 'min_max'
norm_type = 'z_score'

dataset = 'microarray'
# dataset = 'rna_seq'

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

microarray_data_pattern = os.path.join(file_dir, f'./2024-08-08_microarray/{norm_type}/exp/*.csv')

rna_seq_data_pattern = os.path.join(file_dir, f'./2024-08-08_rna_seq/{norm_type}/exp/*.csv')

microarray_stats_files = glob.glob(microarray_data_pattern)

rna_seq_stats_files = glob.glob(rna_seq_data_pattern)

microarray_stats_df = dataframe_utils.csv_files_to_df(microarray_stats_files)
microarray_stats_df.drop_duplicates(inplace=True, ignore_index=True)

rna_seq_stats_df = dataframe_utils.csv_files_to_df(rna_seq_stats_files)
rna_seq_stats_df.drop_duplicates(inplace=True, ignore_index=True)

microarray_stats_df['fraction_clusters_enriched'] = microarray_stats_df.apply(clustering_utils.compute_fraction_clusters_enriched, axis=1)

rna_seq_stats_df['fraction_clusters_enriched'] = rna_seq_stats_df.apply(clustering_utils.compute_fraction_clusters_enriched, axis=1)

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

rna_seq_stats_df['iqr_cluster_size'] = rna_seq_stats_df['q3_cluster_size'] - rna_seq_stats_df['q1_cluster_size']
rna_seq_stats_df['iqr_enriched_cluster_size'] = rna_seq_stats_df['q3_enriched_cluster_size'] - rna_seq_stats_df['q1_enriched_cluster_size']
rna_seq_stats_df['cv_cluster_size'] = rna_seq_stats_df['mean_cluster_size'] / rna_seq_stats_df['sd_cluster_size']
rna_seq_stats_df['cv_enriched_cluster_size'] = rna_seq_stats_df['mean_enriched_cluster_size'] / rna_seq_stats_df['sd_enriched_cluster_size']
microarray_stats_df['iqr_cluster_size'] = microarray_stats_df['q3_cluster_size'] - microarray_stats_df['q1_cluster_size']
microarray_stats_df['iqr_enriched_cluster_size'] = microarray_stats_df['q3_enriched_cluster_size'] - microarray_stats_df['q1_enriched_cluster_size']
microarray_stats_df['cv_cluster_size'] = microarray_stats_df['mean_cluster_size'] / microarray_stats_df['sd_cluster_size']
microarray_stats_df['cv_enriched_cluster_size'] = microarray_stats_df['mean_enriched_cluster_size'] / microarray_stats_df['sd_enriched_cluster_size']

from mpl_toolkits.mplot3d import Axes3D

# Enable interactive mode
plt.ion()

x_all_stat = 'modularity'

y_all_stat = 'fraction_clusters_enriched'

# y_all_stat = 'nenriched_clusters'

z_all_stat = 'iqr_cluster_size'

###############

if dataset == 'microarray':
    df = microarray_stats_df

if dataset == 'rna_seq':
    df = rna_seq_stats_df

df = df.reset_index()

data_dict = OrderedDict()

colors = ['#FFB000', '#785EF0', '#648FFF', '#DC267F', '#FE6100']

for idx, m in enumerate(df['metric'].unique()):

    filtered_df = df.loc[df['metric'] == m]

    x_all = filtered_df[x_all_stat].values

    y_all = filtered_df[y_all_stat].values

    z_all = filtered_df[z_all_stat].values

    data_dict[m] = {'x': x_all, 
                    'y': y_all, 
                    'z': z_all, 
                    'color': colors[idx]}

###############

global_s = 10

fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')

for plot_metric, plot_dict in data_dict.items():

    ax.scatter(plot_dict['x'], plot_dict['y'], plot_dict['z'], color=plot_dict['color'], s=global_s, label=plot_metric, zorder=1)

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
ax.legend()

#   for nns_value, group in df.groupby('nns'):
#      group_df = group.copy(deep=True).sort_values(by='modularity')
#      ax.plot(group_df['modularity'], group_df[y_all_stat], color='grey', linestyle='-', linewidth=1, zorder=-1)

ax.set_xlabel('x='+x_all_stat)
ax.set_ylabel('y='+y_all_stat)
ax.set_zlabel('z='+z_all_stat)

plt.savefig(f'3d_{dataset}_{norm_type}.svg', format='svg')

plt.show(block=True)