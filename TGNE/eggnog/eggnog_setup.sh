# First, make sure that everything is install and paths are set as described
# here: https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.6

# This is essentially the download_eggnog_data.py scripts commands, but I had to make
# adjustments for this to work correctly on my machine.

# Make sure to set $EGGNOG_DATA_DIR and add to $PATH

# make sure to be in correct conda environment

# To download the necessary data:
# Don't download Diamond DB

# eggnog db
cd $EGGNOG_DATA_DIR && wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.db.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog.db.gz && echo Decompressing... && gunzip eggnog.db.gz -f

# eggnog taxa db
cd $EGGNOG_DATA_DIR && wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.taxa.tar.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz && echo Decompressing... && tar -zxf eggnog.taxa.tar.gz && rm eggnog.taxa.tar.gz

# pfam db
cd $EGGNOG_DATA_DIR && wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O pfam.tar.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/pfam.tar.gz && echo Decompressing... && tar -zxf pfam.tar.gz && rm pfam.tar.gz

# hmmer db for eukaryotes (note: must pre-install coreutils on mac to use gmv and gzcat)
mkdir $EGGNOG_DATA_DIR/hmmer/2759
cd $EGGNOG_DATA_DIR/hmmer/2759; echo Downloading HMMs... && wget -N -nH --user-agent=Mozilla/5.0 --relative -r --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2759/2759_hmms.tar.gz && echo Decompressing HMMs... && tar zxf 2759_hmms.tar.gz && echo 2759/* | xargs gmv -t ./ && rm -r 2759 && rm 2759_hmms.tar.gz; numf=$(find ./ | grep -c ".hmm$"); curr=0; cat /dev/null > 2759.hmm_tmp; for file in $(find ./ | grep ".hmm$"); do curr=$((curr+1)); echo "merging HMMs... ${file} (${curr}/${numf})"; cat "${file}" | sed -e "s/.faa.final_tree.fa//" -e "s/.faa.final_tree//" >> 2759.hmm_tmp; rm "${file}"; done; mv 2759.hmm_tmp 2759.hmm; (if [ -f 2759.hmm.h3i ]; then rm 2759.hmm.h3*; fi) && echo "hmmpress-ing HMMs... " && hmmpress 2759.hmm && echo "generating idmap file... " && cat 2759.hmm | grep "^NAME" | sed -e "s/^NAME *//" | awk '{print NR"    "$0}' > 2759.hmm.idmap && echo "removing single OG hmm files... " && echo ./*hmm | xargs rm;

# fastas
cd $EGGNOG_DATA_DIR/hmmer/2759; echo Downloading FASTAs... && wget -N -nH --user-agent=Mozilla/5.0 --relative -r --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2759/2759_raw_algs.tar && echo Decompressing FASTAs... && tar xf 2759_raw_algs.tar && echo 2759/* | xargs gmv -t ./ && rm -r 2759 && rm 2759_raw_algs.tar; numf=$(find ./ | grep -c ".faa.gz$"); curr=0; for file in $(find ./ | grep ".faa.gz$"); do curr=$((curr+1)); echo "processing FASTAs...  ${file} (${curr}/${numf})"; outf=$(echo "$file" | sed "s/\.raw_alg\.faa\.gz/\.fa/"); gzcat "$file" | awk '/^>/{print; next}{gsub("-", ""); print}' > "$outf" && rm "$file"; done





