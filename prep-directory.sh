#!/bin/bash
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --mem 2G
#SBATCH -J prep-directory
#SBATCH -p classroom

# This script sets up the directory structure and makes copies or links to the necessary files
# for the RNA-Seq module..

#Create main project directory
cd ~/
mkdir ~/mouse-rnaseq-2020/
cd ~/mouse-rnaseq-2020/

#Set up directory structure. data/ folder will contain raw data, results/ will contain results from analyses, and src/ will contain scripts.
mkdir data results src
mkdir data/rawseq data/genome

# Create soft links (a.k.a shortcuts) to the raw data
cd data/rawseq/
ln -s /home/classroom/hpcbio/mayo-rnaseq/mouse-rnaseq-2020/data/rawseq/*.fastq .
cd ../genome/
ln -s /home/classroom/hpcbio/mayo-rnaseq/mouse-rnaseq-2020/data/genome/mouse_chr12.* . 

# Create copies of scripts
cd ../../src/
cp /home/classroom/hpcbio/mayo-rnaseq/mouse-rnaseq-2020/src/Mayo-RNASeq/* ./

# Create copy of Targets0.txt
cd ../results
cp /home/classroom/hpcbio/mayo-rnaseq/mouse-rnaseq-2020/results/Targets0.txt ./

echo "Finished prep-directory.sh"

