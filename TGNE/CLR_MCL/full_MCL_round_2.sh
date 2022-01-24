# First, take a look at what different z-score cutoffs will do to the network
mcx query -imx full_clr_graph_round_2.mci -vary-threshold 7/20/130

# It looks like 14.4 is as high as we can go without really losing genes from the network.

mcx alter -imx full_clr_graph_round_2.mci -tf 'gq(14.4), add(-14.4)' -o full_clr_graph_round_2_14.4.mci

# To assess the quality of the results, first perform coarse clustering
clm close -imx full_clr_graph_round_2_14.4.mci --write-cc -o full_clr_graph_round_2_14.4.base

# Cluster and check stats. We want to maximimize efficiency while keeping things reasonable
mcl full_clr_graph_round_2_14.4.mci -I 1.2 -scheme 7

mcl full_clr_graph_round_2_14.4.mci -I 1.4 -scheme 7

mcl full_clr_graph_round_2_14.4.mci -I 1.6 -scheme 7

mcl full_clr_graph_round_2_14.4.mci -I 1.8 -scheme 7

mcl full_clr_graph_round_2_14.4.mci -I 2 -scheme 7

mcl full_clr_graph_round_2_14.4.mci -I 2.2 -scheme 7

mcl full_clr_graph_round_2_14.4.mci -I 2.4 -scheme 7

mcl full_clr_graph_round_2_14.4.mci -I 2.6 -scheme 7

mcl full_clr_graph_round_2_14.4.mci -I 2.8 -scheme 7

mcl full_clr_graph_round_2_14.4.mci -I 3 -scheme 7

clm info full_clr_graph_round_2_14.4.mci out.full_clr_graph_round_2_14.4.mci.I{12,14,16,18,20,22,24,26,28,30}
clm dist --chain out.full_clr_graph_round_2_14.4.mci.I{12,14,16,18,20,22,24,26,28,30}

# Based on the summary stats for the different approaches, it looks like -I 1.8 is the way to go here.
# We get efficiency=0.41976 massfrac=0.65759 areafrac=0.00965 The MCL algorithm also gives "superb"
# pruning jury stats (98,96,98 out of 100). This gives 443 clusters, most of which are very small.
# Using a higher inflation does give better efficiency, but the quality really starts to drop off
# below 1.8. I think this gives a healthy balance between wanting modules of reasonable size and
# having the most performant clustering.

mcxdump -icl out.full_clr_graph_round_2_14.4.mci.I18 -o full_labeled_MCL_round_2_clusters.txt -tabr full_clr_graph_round_2.dict


# repeat process for vegetative and sexual cases!
mcx query -imx veg_clr_graph_round_2.mci -vary-threshold 7/20/130

# Looks like a cutoff of 10.9 is the way to go here.
mcx alter -imx veg_clr_graph_round_2.mci -tf 'gq(10.9), add(-10.9)' -o veg_clr_graph_round_2_10.9.mci

mcl veg_clr_graph_round_2_10.9.mci -I 1.2 -scheme 7

mcl veg_clr_graph_round_2_10.9.mci -I 1.4 -scheme 7

mcl veg_clr_graph_round_2_10.9.mci -I 1.6 -scheme 7

mcl veg_clr_graph_round_2_10.9.mci -I 1.8 -scheme 7

mcl veg_clr_graph_round_2_10.9.mci -I 2 -scheme 7

mcl veg_clr_graph_round_2_10.9.mci -I 2.2 -scheme 7

mcl veg_clr_graph_round_2_10.9.mci -I 2.4 -scheme 7

mcl veg_clr_graph_round_2_10.9.mci -I 2.6 -scheme 7

mcl veg_clr_graph_round_2_10.9.mci -I 2.8 -scheme 7

mcl veg_clr_graph_round_2_10.9.mci -I 3 -scheme 7

clm info veg_clr_graph_round_2_10.9.mci out.veg_clr_graph_round_2_10.9.mci.I{12,14,16,18,20,22,24,26,28,30}
clm dist --chain out.veg_clr_graph_round_2_10.9.mci.I{12,14,16,18,20,22,24,26,28,30}

# Here, I also see that -I 1.8 is the way to go: efficiency=0.38906 massfrac=0.60813 areafrac=0.01712
mcxdump -icl out.veg_clr_graph_round_2_10.9.mci.I18 -o veg_labeled_MCL_round_2_clusters.txt -tabr veg_clr_graph_round_2.dict 

# Now for the sexual case
mcx query -imx sex_clr_graph_round_2.mci -vary-threshold 7/20/130

# Use a cutoff of 14.9
mcx alter -imx sex_clr_graph_round_2.mci -tf 'gq(14.9), add(-14.9)' -o sex_clr_graph_round_2_14.9.mci

mcl sex_clr_graph_round_2_14.9.mci -I 1.2 -scheme 7

mcl sex_clr_graph_round_2_14.9.mci -I 1.4 -scheme 7

mcl sex_clr_graph_round_2_14.9.mci -I 1.6 -scheme 7

mcl sex_clr_graph_round_2_14.9.mci -I 1.8 -scheme 7

mcl sex_clr_graph_round_2_14.9.mci -I 2 -scheme 7

mcl sex_clr_graph_round_2_14.9.mci -I 2.2 -scheme 7

mcl sex_clr_graph_round_2_14.9.mci -I 2.4 -scheme 7

mcl sex_clr_graph_round_2_14.9.mci -I 2.6 -scheme 7

mcl sex_clr_graph_round_2_14.9.mci -I 2.8 -scheme 7

mcl sex_clr_graph_round_2_14.9.mci -I 3 -scheme 7

clm info sex_clr_graph_round_2_14.9.mci out.sex_clr_graph_round_2_14.9.mci.I{12,14,16,18,20,22,24,26,28,30}
clm dist --chain out.sex_clr_graph_round_2_14.9.mci.I{12,14,16,18,20,22,24,26,28,30}

# Sticking with -I 1.8 again efficiency=0.41547 massfrac=0.62523 areafrac=0.01015
mcxdump -icl out.sex_clr_graph_round_2_14.9.mci.I18 -o sex_labeled_MCL_round_2_clusters.txt -tabr sex_clr_graph_round_2.dict 