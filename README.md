# TGNE-2022
Analysis for and development of the Tetrahymena Gene Network Explorer

## User Guide:

### Building the conda environment for python and R: 

Navigate to the top level directory and run one of the following:

Mac – Apple Silicon ([Rosetta](https://support.apple.com/en-us/102527) required):
```
bash build_tgne_env_silicon.sh
```

Linux or Mac – Intel:
```
bash build_tgne_env.sh
```

### Building the TGNE from precomputed data:

Download the [precomputed data](LINK # FIXME).

Create necessary directories within the top level directory:

```
mkdir -p input_data/eggnog_annotations_file/
mkdir -p input_data/interproscan_tsv_file/
mkdir -p active_files/
```

Place ```eggnog_compiled_2024_jun28.annotations``` in the ```input_data/eggnog_annotations_file/``` directory.

Place ```interproscan_compiled.tsv``` in the ```input_data/interproscan_tsv_file/``` directory.

Place ```allgood_filt_agg_tidy_2021aligned_qc_rma_expression_full.csv``` in the ```active_files/``` directory.

Place ```rna_seq.csv``` in the ```active_files/``` directory.

```
bash s-execute-pipeline.sh pipeline_precomputed.txt
```

### Building the TGNE from raw data:

1) Kallisto:

Create necessary directories within the top level directory:

```
mkdir -p input_data/kallisto_data_folders
```

Build the conda environment.

Activate the conda environment:
```
conda activate tgne.env
```

Install the [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit). Ensure that the toolkit is in your system path.

Download the [RNA-seq data from SRA](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA861835&o=acc_s%3Aa) (this will take a few hours):

```
prefetch PRJNA861835
```

The ```declare -A``` flag requires bash version > 4.
```
declare -A time_points=(
    [SRR20576712]="240min_B"
    [SRR20576713]="240min_A"
    [SRR20576714]="210min_B"
    [SRR20576715]="210min_A"
    [SRR20576716]="180min_B"
    [SRR20576717]="180min_A"
    [SRR20576718]="150min_B"
    [SRR20576719]="150min_A"
    [SRR20576720]="120min_B"
    [SRR20576721]="120min_A"
    [SRR20576722]="090min_B"
    [SRR20576723]="090min_A"
    [SRR20576724]="060min_B"
    [SRR20576725]="060min_A"
    [SRR20576726]="030min_B"
    [SRR20576727]="030min_A"
    [SRR20576728]="000min_B"
    [SRR20576729]="000min_A"
)

for i in {12..29}; do
    srr_id="SRR205767${i}"
    
    time_sample="${time_points[$srr_id]}"
    
    fasterq-dump -e ${num_cores} -o ${time_sample}_# -S $srr_id
    
    gzip ./${time_sample}*.fastq
done
```

Or (if you have [pigz](https://zlib.net/pigz/) installed):

```
num_cores=<NUMBER OF CORES TO USE>

declare -A time_points=(
    [SRR20576712]="240min_B"
    [SRR20576713]="240min_A"
    [SRR20576714]="210min_B"
    [SRR20576715]="210min_A"
    [SRR20576716]="180min_B"
    [SRR20576717]="180min_A"
    [SRR20576718]="150min_B"
    [SRR20576719]="150min_A"
    [SRR20576720]="120min_B"
    [SRR20576721]="120min_A"
    [SRR20576722]="090min_B"
    [SRR20576723]="090min_A"
    [SRR20576724]="060min_B"
    [SRR20576725]="060min_A"
    [SRR20576726]="030min_B"
    [SRR20576727]="030min_A"
    [SRR20576728]="000min_B"
    [SRR20576729]="000min_A"
)

for i in {12..29}; do
    srr_id="SRR205767${i}"
    
    time_sample="${time_points[$srr_id]}"
    
    fasterq-dump -e ${num_cores} -o ${time_sample}_# -S $srr_id
    
    pigz -p ${num_cores} ./${time_sample}*.fastq
done
```

Move the samples to a new folder:
```
mkdir -p raw/
mv *.fastq.gz raw/
```

Run FastQC on each of the fastq.gz files:
```
cd raw
mkdir -p qc/
for filename in *.fastq.gz; do fastqc --threads 16 --outdir ./qc $filename; done
```

Compile the FastQC results with MultiQC:
```
multiqc ./qc/ -o ./qc/
```

Download the proper Trimmomatic adapter sequence:
```
curl -O https://raw.githubusercontent.com/usadellab/Trimmomatic/refs/heads/main/adapters/TruSeq3-PE.fa
```

Run Trimmomatic on each pair of fastq.gz files: 
```
output_folder="trimmed/"

mkdir -p ${output_folder}

file_name_1=<FILE_NAME_1>
file_name_2=<FILE_NAME_2>

trimmomatic PE ${FASTQ_DIR%/}/${file_name_1} ${FASTQ_DIR%/}/${file_name_2} ${output_folder%/}/p_trimmed_${file_name_1} ${output_folder%/}/up_trimmed_${file_name_1} ${output_folder%/}/p_trimmed_${file_name_2} ${output_folder%/}/up_trimmed_${file_name_2} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 -threads <NUM_THREADS>
```

Move the unpaired reads to a new folder:
```
cd trimmed/

mkdir -p up/

mv up_trimmed_* up/
```

Run FastQC on each of the trimmed fastq.gz files:
```
mkdir -p qc/
fastqc --threads <NUM_THREADS> --outdir ./qc/ ./<FILE_NAME>
```

Compile the FastQC results with MultiQC:
```
multiqc ./qc/ -o ./qc/
```

Download the latest CDS file:
```
# Uncomment the two lines below if on a cluster
# mkdir -p cds_index/
# curl -o ./cds_index/cds.fasta https://github.com/yefei521/Tetrahymena_Genome_annotation_V2024/releases/download/V2024.2/Tetrahymena_Genome_annotation_V2024_CDS.fasta

# Comment two three lines below if on a cluster
wget https://github.com/yefei521/Tetrahymena_Genome_annotation_V2024/releases/download/V2024.2/Tetrahymena_Genome_annotation_V2024_CDS.fasta
mv ./Tetrahymena_Genome_annotation_V2024_CDS.fasta ./cds_index/cds.fasta
```

Index the latest CDS file with Kallisto:
```
kallisto index --index ./cds_index/kallisto_cds_index ./cds_index/cds.fasta
```

Run Kallisto on each pair of trimmed fastq.gz files: 
```
output_folder="counts/"

mkdir -p ${output_folder}

file_name_1=<FILE_NAME_1>
file_name_2=<FILE_NAME_2>

file_basename=basename="${file_name_1%_*}"

kallisto quant -t {num_cores} -o ${output_folder%/}/kallisto_quant_${file_basename} -i ./cds_index/kallisto_cds_index ${file_name_1} ${file_name_2}
```

Place each of the resulting Kallisto output folders in the ```counts/``` directory in the ```input_data/kallisto_data_folders``` directory.

2) Eggnog:

Create necessary directories within the top level directory:

```
mkdir -p input_data/eggnog_annotations_file
```

Build the conda environment.

Activate the conda environment:
```
conda activate tgne.env
```

Download the latest protein sequences file:
```
curl -o ./pep.fasta https://github.com/yefei521/Tetrahymena_Genome_annotation_V2024/releases/download/V2024.2/Tetrahymena_Genome_annotation_V2024_Protein_addAnno.fasta
```

Prepare the latest protein sequences file:

Place ```scripts/eggnog_pep_prep.py``` in your current directory.

```
python eggnog_pep_prep.py
```

Download the Eggnog database:
```
mkdir -p eggnog_data/

download_eggnog_data.py --data_dir eggnog_data/ -P -H -d 2759
```

Respond to the prompts:
Download main annotation database? [y,n] ```y```
Download taxa database? [y,n] ```y```
Download diamond database (~4GB after decompression)? [y,n] ```n```
Download pfam database (~3GB after decompression)? [y,n] ```y```
Download HMMER database of tax ID 2759? [y,n] ```y```
Please, specify a non-empty name for the database (e.g. Bacteria) [default:2759]: ```2759```

Use Eggnog to gather ortholog hits:
```
emapper.py --data_dir /eggnog_data -m hmmer -d 2759 --no_annot -i ./pep.fasta -o pep_a --cpu 12 --tax_scope Eukaryota --target_taxa 2759 --report_orthologs --report_no_hits --go_evidence non-electronic --pfam_realign none --dbtype hmmdb --override
```

Use Eggnog to gather orthologous annotations:
```
emapper.py --data_dir /project2/apturkew/michael/eggnog_data -m no_search --annotate_hits_table pep_a.emapper.seed_orthologs -o pep_b --dbmem --excel --override
```

Place the resulting Eggnog output file, ```pep_b.emapper.annotations```, in the ```input_data/eggnog_annotations_file/``` directory.

3) InterProScan:

Create necessary directories within the top level directory:

```
mkdir -p input_data/eggnog_annotations_file
```

Install [InterProScan](https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html) (Linux only).

Build the conda environment.

Activate the conda environment:
```
conda activate tgne.env
```

Download the latest protein sequences file:
```
curl -o ./pep.fasta https://github.com/yefei521/Tetrahymena_Genome_annotation_V2024/releases/download/V2024.2/Tetrahymena_Genome_annotation_V2024_Protein_addAnno.fasta
```

Prepare the latest protein sequences file:

Place ```scripts/interproscan_pep_prep.py``` in your current directory.

```
python interproscan_pep_prep.py
```

```
path_to_interproscan=<PATH TO AND INLCUDING interproscan-5.68-100.0/>

${path_to_interproscan%/}/interproscan.sh -i pep_cleaned.fasta -f tsv -d ./ -cpu <NUM_CPUS>
```

Place the resulting InterProScan output file, ```pep_cleaned.fasta.tsv```, in the ```input_data/eggnog_annotations_file/``` directory.

4. Build the TGNE:

Create necessary directories within the top level directory:

```
mkdir -p active_files/
```

```
bash s-execute-pipeline.sh pipeline_precomputed.txt
```



### Mucocyst regranulation analysis:

Build the TGNE following "Building the TGNE from precomputed data" or "Building the TGNE from raw data."

Navigate to the top level directory and run the following:

```
bash s-execute-pipeline.sh ./regranulation/regranulation_pipeline.txt
```
