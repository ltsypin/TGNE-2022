# TGNE-2022
Analysis for and development of the Tetrahymena Gene Network Explorer

## User Guide:

### Building the conda environment for python and R. Navigate to the top level directory and run on of the following:

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

Place each of the following files in the specificed directory:

```
eggnog_compiled_2024_jun28.annotations -> input_data/eggnog_annotations_file/
```

```
interproscan_compiled.tsv -> input_data/interproscan_tsv_file/
```

```
allgood_filt_agg_tidy_2021aligned_qc_rma_expression_full.csv -> active_files/
```

```
rna_seq.csv -> active_files/
```

```
bash s-execute-pipeline.sh -> pipeline_precomputed.txt
```

### Building the TGNE from raw data:

1) Kallisto:

Install the [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit).

Download the [RNA-seq data from SRA](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA861835&o=acc_s%3Aa):

```
prefetch PRJNA861835
```

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

Install FastQC.

Run FastQC on each of the fastq.gz files:
```
cd raw
mkdir -p qc/
fastqc --threads <NUM_THREADS> --outdir ./qc/ ./<FILE_NAME>
```

Install MultiQC.

Compile the FastQC results with MultiQC:
```
multiqc ./qc/ -o ./qc/
```

Install Trimmomatic.

Run Trimmomatic on each pair of fastq.gz files: 
```
output_folder="trimmed/"

mkdir -p ${output_folder}

trimmomatic_path=<PATH TO Trimmomatic-0.39>

file_name_1=<FILE_NAME_1>
file_name_2=<FILE_NAME_2>

java -jar ${trimmomatic_path%/}/trimmomatic-0.39.jar PE ${FASTQ_DIR%/}/${file_name_1} ${FASTQ_DIR%/}/${file_name_2} ${output_folder%/}/p_trimmed_${file_name_1} ${output_folder%/}/up_trimmed_${file_name_1} ${output_folder%/}/p_trimmed_${file_name_2} ${output_folder%/}/up_trimmed_${file_name_2} ILLUMINACLIP:${trimmomatic_path%/}/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 -threads <NUM_THREADS>
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

Install MultiQC.

Compile the FastQC results with MultiQC:
```
multiqc ./qc/ -o ./qc/
```

Install Kallisto.

Download the latest CDS file:
```
mkdir -p cds_index/
curl -o ./cds_index/cds.fasta https://github.com/yefei521/Tetrahymena_Genome_annotation_V2024/releases/download/V2024.2/Tetrahymena_Genome_annotation_V2024_CDS.fasta
```

Index the latest CDS file with Kallisto:
```
kallisto index --index ./cds_index/kallisto_cds_index ./cds_index/cds.fasta
```

Run Kallisto on each pair of trimmed fastq.gz files: 
```
output_folder="counts/"

mkdir -p ${output_folder}

trimmomatic_path=<PATH TO Trimmomatic-0.39>

file_name_1=<FILE_NAME_1>
file_name_2=<FILE_NAME_2>

file_basename=basename="${file_name_1%_*}"

kallisto quant -t {num_cores} -o ${output_folder%/}/kallisto_quant_${file_basename} -i ./cds_index/kallisto_cds_index ${file_name_1} ${file_name_2}
```

Place each of the resulting Kallisto output folders in the ```counts/``` directory in the ```input_data/kallisto_data_folders``` directory.

2) Eggnog:


3) InterProScan:



