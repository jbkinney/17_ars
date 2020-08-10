#!/bin/bash
# specify BASH shell
#$ -S /bin/bash
# run job in the current working directory where qsub is executed from
#$ -cwd
#  specify that the job requires 8GB of memory
#$ -l m_mem_free=8G

# run commands and application
pwd
date
./map_reads.sh reads/Pd7-G_S7_R1_001.fastq  reads/Pd7-G_S7_R2_001.fastq  genome/s288c.fa tmp/ edu_exp_2_WT.txt stats_WT.txt
date
