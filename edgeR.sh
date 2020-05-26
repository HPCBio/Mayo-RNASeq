#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 8G
#SBATCH --job-name=edgeR
#SBATCH -p classroom

#Change directory to the main one containing data/, results/ and src/
cd ~/mouse-rnaseq-2020/

#Load the R module; name may change in future so check each time with 'module avail' on command line
module load R/3.6.0-IGB-gcc-8.2.0


#Create output directory for the results.
# NOTE: the -p option will not only create parental directories if needed, but also will not
# throw an error or overwrite the directory if it already exists
mkdir -p results/edgeR

#run the R script
echo "start edgeR"

Rscript src/stats_edgeR.R

echo "end edgeR; check results/edgeR for output files"
