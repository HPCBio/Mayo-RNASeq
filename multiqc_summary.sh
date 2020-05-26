#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 8G
#SBATCH --job-name=multiqc
#SBATCH -p classroom

#Change directory to the main results 
cd ~/mouse-rnaseq-2020/results

# Load the multiqc module 

module load MultiQC/1.7-IGB-gcc-8.2.0-Python-3.7.2

# Run multiqc on the star and featureCounts outputs

echo "start multiqc"

multiqc ./

echo "end multiqc"

module purge


#Load the R module; name may change in future so check each time with 'module avail' on command line
module load R/3.6.0-IGB-gcc-8.2.0


#run the R script
echo "start makeTargetsFinal.R"

Rscript ../src/makeTargetsFinal.R

echo "end makeTargetsFinal.R"

echo "download and inspect the multiqc_report.html, ReadFatePlot.jpeg and Targets_Final.txt"

