#!/bin/bash
#SBATCH -N 1 
#SBATCH -n 2
#SBATCH --mem 32G
#SBATCH -J make.index
#SBATCH -p classroom


#Change directory to the main one containing data/, results/ and src/.
cd ~/mouse-rnaseq-2020/

#Load the STAR module; name may change in future so check each time with 'module avail STAR' on command line
module load STAR/2.7.3a-IGB-gcc-8.2.0

#Create output directory for the STAR index
# NOTE: the -p option will not only create parental directories if needed, but also will not
# throw an error or overwrite the directory if it already exists
mkdir -p data/genome/STAR-2.7.3a_mouse-chr12_Index/

#Make the mouse STAR index
date
STAR --runThreadN $SLURM_NTASKS \
     --runMode genomeGenerate \
     --genomeDir data/genome/STAR-2.7.3a_mouse-chr12_Index \
     --genomeFastaFiles data/genome/mouse_chr12.fna \
     --limitGenomeGenerateRAM 32000000000 \
     --genomeSAindexNbases 12 \
     --outTmpDir /scratch/$SLURM_JOB_ID
date

