#!/bin/sh
#$ -cwd     # Use current working directory
#$ -V     # Verbose
#$ -j y     # Maximum output, inc errors
#$ -r y     # Condense error files
#$ -pe smp 1     # Request CPU cores
#$ -l h_rt=20:00:00     # Request runtime (up to 240 hours)
#$ -l h_vmem=20G     # Request RAM per core
#$ -m bea     # Status emails
#$ -t 91-98



module load R

#INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" /data/home/btw863/network_OTU_comparisons/numbers.csv)

Rscript /data/scratch/btw863/bat-diet/scripts/r/array_lulu.R ${SGE_TASK_ID}


