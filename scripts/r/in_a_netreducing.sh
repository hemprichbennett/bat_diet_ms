#!/bin/sh
#$ -cwd     # Use current working directory
#$ -V     # Verbose
#$ -j y     # Maximum output, inc errors
#$ -r y     # Condense error files
#$ -pe smp 1     # Request CPU cores
#$ -l h_rt=150:0:0     # Request runtime (up to 240 hours)
#$ -l h_vmem=10G     # Request RAM per core
#$ -m bea     # Status emails
#$ -t 1-14


module load R/3.4.3



Rscript /data/scratch/btw863/bat-diet/scripts/r/netreducing.R ${SGE_TASK_ID}
