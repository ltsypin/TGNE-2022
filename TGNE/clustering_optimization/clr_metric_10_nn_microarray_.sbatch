#!/bin/bash
#SBATCH --job-name=run_scan_clr_10_microarray
#SBATCH --account=pi-apturkew
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2000
##SBATCH --output=myjob_%j.out      # Standard output log
##SBATCH --error=myjob_%j.err       # Standard error log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bertagna@rcc.uchicago.edu

source activate cdh2_bokeh.env

cmd="python3 fast_scan_modular_metrics.py microarray clr 10"

echo $cmd

eval $cmd
