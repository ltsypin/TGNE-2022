#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import pandas as pd
from statistics import mean, stdev
import sys


# In[2]:


files_consensus_partitions = glob.glob('./TAG/*.labels')

files_consensus_partitions.sort(key=lambda name: int((name.split('res')[1]).split('.')[0]))

files_consensus_partitions


# In[3]:


list_consensus_partitions = {}

for file in files_consensus_partitions:
    consensus_partition = []
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            consensus_partition.append(line.split('\t'))
    list_consensus_partitions[file.split('.')[2]] = consensus_partition


# In[4]:


for name, partition in list_consensus_partitions.items():
    print(name)
    print('                  NUM CLUSTERS:', len(partition))

    max_cluster_size = -1

    all_cluster_lens = []

    len_ge_100 = 0
    len_lt_100 = 0
    len_lt_50 = 0
    len_lt_10 = 0
    len_eq_1 = 0

    for cluster in partition:
        if len(cluster) > max_cluster_size:
            max_cluster_size = len(cluster)

        all_cluster_lens.append(len(cluster))
        if len(cluster) >= 100:
            len_ge_100 += 1
        if len(cluster) < 100:
            len_lt_100 += 1
        if len(cluster) < 50:
            len_lt_50 += 1
        if len(cluster) < 10:
            len_lt_10 += 1
        if len(cluster) == 1:
            len_eq_1 += 1

    mean_len = mean(all_cluster_lens)
    stdev_len = stdev(all_cluster_lens)

    cv_len = (stdev_len / mean_len) * 100

    print('               CLUSTER SIZE CV:', cv_len)
    print('              MAX CLUSTER SIZE:', max_cluster_size)
    print('NUM CLUSTERS WITH >= 100 GENES:', len_ge_100)
    print(' NUM CLUSTERS WITH < 100 GENES:', len_lt_100)
    print('  NUM CLUSTERS WITH < 50 GENES:', len_lt_50)
    print('  NUM CLUSTERS WITH < 10 GENES:', len_lt_10)
    print('  NUM CLUSTERS WITH == 1 GENES:', len_eq_1)
    print()


# In[5]:


df = pd.DataFrame(columns=['TTHERM_ID', 'leiden_label_full'])

cluster_num = 0

chosen_resolution = ('res' + str(sys.argv[1]))

for cluster in list_consensus_partitions[chosen_resolution]:
    for gene in cluster:
        new_row = {'TTHERM_ID': gene.strip(), 'leiden_label_full': cluster_num}
        df = df.append(new_row, ignore_index=True)
    cluster_num += 1

df.reset_index(inplace=True, drop=True)
df


# In[6]:


df.to_csv('./mcl_rcl_label_df_round_1.csv', index=False)


# # In[16]:


# gene_list_test = []
# for cluster in list_consensus_partitions['res400']:
#     if 'TTHERM_00321680' in cluster:
#         print(len(cluster))
#         for gene in cluster:
#             gene_list_test.append(gene.strip())
#             print('\'',gene.strip(),'\'', end=', ', sep='')


# # In[8]:


# for cluster in list_consensus_partitions['res400']:
#     if 'TTHERM_00527180' in cluster:
#         print(len(cluster))
#         for gene in cluster:
#             print(gene.strip(), end=', ')


# # In[13]:


# for gene in gene_list_test:
#     print((df.loc[df['TTHERM_ID'] == gene])['leiden_label_full'].values[0])

