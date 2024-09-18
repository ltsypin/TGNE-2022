source ~/.bash_profile

conda clean --all -y

CONDA_SUBDIR=osx-64 conda create -n tgne.env python=3.8.16 -c conda-forge --override-channels -y

conda activate tgne.env

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda config --env --set subdir osx-64

conda install pandas matplotlib matplotlib-venn numpy scikit-learn umap-learn python-igraph leidenalg networkx scipy bokeh openpyxl distinctipy jupyter gdown tqdm biopython=1.76 psutil=5.7.0 xlsxwriter=1.4.3 eggnog-mapper r-base r-ggplot2 r-bigmemory r-snow r-rmarkdown r-tidyr r-rcolorbrewer r-actuar r-fitdistrplus r-hmisc r-wgcna r-biocmanager bioconductor-geneplotter bioconductor-oligo bioconductor-pdinfobuilder bioconductor-masigpro bioconductor-limma r-dplyr r-ggrepel r-kernsmooth r-parmigene -y

Rscript -e 'install.packages("admisc", dependencies=TRUE, repos="https://cran.rstudio.com/")'
