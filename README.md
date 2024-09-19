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
mkdir raw
mv *.fastq.gz raw/
```



Install FastQC.

Install MultiQC.

Install Trimmomatic.

Install Kallisto.

Set paths and SLURM details in scripts 
FastQC MultiQC
Trim the reads
FastQC MultiQC
Index the CDS
Kallisto

2) Eggnog:


3) InterProScan:



