# TGNE-2022
Analysis for and development of the Tetrahymena Gene Network Explorer

Brief guide to starting:

1) Build the conda environments for python and R (provided spec files called cdh2 and cdh2_r). 
Note: for the R environment, will still need to install some packages using BiocManager; there 
might be some pip packages for the python environment...

You will also need to install hisat2 separately and make sure the PATH works.

2) Download the raw data by running the notebook in the top level directory 
(Download_raw_data.ipynb).

3) Begin by running 
TGNE/microarray_probe_alignment_and_filtering/microarray_probe_alignment.ipynb until you get 
to the microarray QC step

4) Run TGNE/microarray_QC/microarray_QC.Rmd in the appropriate conda environment

5) Return to the mmnicroarray_probe_alignment notebook and go through the rest of it

6) Run TGNE/microarray_probe_alignment_and_filtering/gene_filtering.Rmd

7) Run TGNE/embedding/embedding_and_app.ipynb and follow directions for when to run the 
TGNE/enrichment/enrichment.ipynb notebook and when to return back to the embedding file

There are a bunch of frayed ends in nearly all of these notebooks, and we can talk about where 
they're supposed to go.
