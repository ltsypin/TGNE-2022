# First, make sure that everything is install and paths are set as described
# here: https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.6

# To download the necessary data:
# Don't download Diamond DB
# Keep default as the DB (for some reason specifying Eukaryota didn't work)
download_eggnog_data.py -H -d 2759 # fails after downloading hmmer data because of differences between GNU and Mac OS X mv and zcat commands

# Install coreutils using brew. Rerun failed commands with gmv and gzcat
cd /Users/eukarya/eggnog-mapper-data/hmmer/2759; echo Downloading HMMs... && wget -N -nH --user-agent=Mozilla/5.0 --relative -r --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2759//2759_hmms.tar.gz && echo Decompressing HMMs... && tar zxf 2759_hmms.tar.gz && echo 2759/* | xargs gmv -t ./ && rm -r 2759 && rm 2759_hmms.tar.gz; numf=$(find ./ | grep -c ".hmm$"); curr=0; cat /dev/null > 2759.hmm_tmp; for file in $(find ./ | grep ".hmm$"); do curr=$((curr+1)); echo "merging HMMs... ${file} (${curr}/${numf})"; cat "${file}" | sed -e "s/.faa.final_tree.fa//" -e "s/.faa.final_tree//" >> 2759.hmm_tmp; rm "${file}"; done; mv 2759.hmm_tmp 2759.hmm; (if [ -f 2759.hmm.h3i ]; then rm 2759.hmm.h3*; fi) && echo "hmmpress-ing HMMs... " && /anaconda3/envs/eggnog/bin/hmmpress 2759.hmm && echo "generating idmap file... " && cat 2759.hmm | grep "^NAME" | sed -e "s/^NAME *//" | awk '{print NR"    "$0}' > 2759.hmm.idmap && echo "removing single OG hmm files... " && echo ./*hmm | xargs rm;

cd /Users/eukarya/eggnog-mapper-data/hmmer/2759; echo Downloading FASTAs... && wget -N -nH --user-agent=Mozilla/5.0 --relative -r --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2759//2759_raw_algs.tar && echo Decompressing FASTAs... && tar xf 2759_raw_algs.tar && echo 2759/* | xargs gmv -t ./ && rm -r 2759 && rm 2759_raw_algs.tar; numf=$(find ./ | grep -c ".faa.gz$"); curr=0; for file in $(find ./ | grep ".faa.gz$"); do curr=$((curr+1)); echo "processing FASTAs...  ${file} (${curr}/${numf})"; outf=$(echo "$file" | sed "s/\.raw_alg\.faa\.gz/\.fa/"); gzcat "$file" | awk '/^>/{print; next}{gsub("-", ""); print}' > "$outf" && rm "$file"; done

# Get the pfams
download_eggnog_data.py -P -f

# Run eggnog-mapper
emapper.py  -m hmmer -d 2759 -i ../../raw_data/Tthermophila_MAC_protein_2021.fasta -o Ttherm2021MAC_eggnog_hmmer_default_evalue --cpu 0 --tax_scope Eukaryota --target_taxa 2759 --report_orthologs --report_no_hits --go_evidence non-electronic --pfam_realign realign --dbtype hmmdb --dbmem



