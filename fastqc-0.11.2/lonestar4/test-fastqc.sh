#!/bin/bash
#$ -V
#$ -cwd # Start job in submission directory
#$ -N test-fastqc # Job Name
#$ -o $JOB_NAME.o$JOB_ID        # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way 12
#$ -q development       # Queue name normal
#$ -A iPlant-Collabs
#$ -l h_rt=01:00:00     # Run time (hh:mm:ss) 
#$ -M jcarson@tacc.utexas.edu   # Address for email notification
#$ -m be        # Email at Begin and End of job

# for DNA Subway w/ Agave

# inputs (for testing)
SEQ1="/work/02570/jcarson/iPlant/testdata/WT_rep1_1.fastq.gz"

tar zxf FastQC.tgz

export PATH=$PWD/FastQC:$PATH

mkdir fastqc_out
FastQC/fastqc --quiet -f fastq -o fastqc_out "$SEQ1"

rm -rf FastQC

# don't delete input during testing
#rm "$SEQ1"


