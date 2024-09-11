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

### Building the TGNE from [precomputed data](LINK):

Download the precomputed data

Deposit the file in the folders

```
bash s-execute-pipeline.sh pipeline_precomputed.txt
```

### Building the TGNE from raw data:


1) Download the [RNA-seq data from SRA](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA861835&o=acc_s%3Aa):

```
prefetch PRJNA861835
```

```
num_cores=<NUMBER OF CORES TO USE>

fasterq-dump -e ${num_cores} --outdir ./ SRR20576712
gzip ./SRR20576712*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576713
gzip ./SRR20576713*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576714
gzip ./SRR20576714*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576715
gzip ./SRR20576715*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576716
gzip ./SRR20576716*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576717
gzip ./SRR20576717*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576718
gzip ./SRR20576718*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576719
gzip ./SRR20576719*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576720
gzip ./SRR20576720*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576721
gzip ./SRR20576721*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576722
gzip ./SRR20576722*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576723
gzip ./SRR20576723*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576724
gzip ./SRR20576724*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576725
gzip ./SRR20576725*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576726
gzip ./SRR20576726*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576727
gzip ./SRR20576727*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576728
gzip ./SRR20576728*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576729
gzip ./SRR20576729*.fastq
```

OR 

```
num_cores=<NUMBER OF CORES TO USE>

fasterq-dump -e ${num_cores} --outdir ./ SRR20576712
pigz -p ${num_cores} ./SRR20576712*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576713
pigz -p ${num_cores} ./SRR20576713*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576714
pigz -p ${num_cores} ./SRR20576714*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576715
pigz -p ${num_cores} ./SRR20576715*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576716
pigz -p ${num_cores} ./SRR20576716*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576717
pigz -p ${num_cores} ./SRR20576717*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576718
pigz -p ${num_cores} ./SRR20576718*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576719
pigz -p ${num_cores} ./SRR20576719*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576720
pigz -p ${num_cores} ./SRR20576720*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576721
pigz -p ${num_cores} ./SRR20576721*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576722
pigz -p ${num_cores} ./SRR20576722*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576723
pigz -p ${num_cores} ./SRR20576723*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576724
pigz -p ${num_cores} ./SRR20576724*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576725
pigz -p ${num_cores} ./SRR20576725*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576726
pigz -p ${num_cores} ./SRR20576726*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576727
pigz -p ${num_cores} ./SRR20576727*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576728
pigz -p ${num_cores} ./SRR20576728*.fastq

fasterq-dump -e ${num_cores} --outdir ./ SRR20576729
pigz -p ${num_cores} ./SRR20576729*.fastq
```

2) Kallisto:
Install FastQC, MultiQC, Trimmomatic, and Kallisto
Set paths and SLURM details in scripts 
FastQC MultiQC
Trim the reads
FastQC MultiQC
Index the CDS
Kallisto

3) Eggnog:


4) InterProScan:

