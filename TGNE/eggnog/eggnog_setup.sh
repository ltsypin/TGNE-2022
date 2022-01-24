# First, make sure that everything is install and paths are set as described
# here: https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.6

# To download the necessary data:
# Don't download Diamond DB
# specify 2759 as the DB (for some reason Eukaryota didn't work)
download_eggnog_data.py -H -d 2759

download_eggnog_data.py -P -f

# Run eggnog-mapper
emapper.py  -m hmmer -d 2759 -i ../../raw_data/Tthermophila_MAC_protein_2021.fasta -o Ttherm2021MAC_eggnog_hmmer_default_evalue --cpu 0 --tax_scope Eukaryota --target_taxa 2759 --report_orthologs --report_no_hits --go_evidence non-electronic --pfam_realign realign --dbtype hmmdb --dbmem



