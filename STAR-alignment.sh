#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem 16G
#SBATCH --job-name=align_star
#SBATCH -p classroom
#SBATCH --array=1-4%2

#Change directory to the main one containing data/, results/ and src/
cd ~/mouse-rnaseq-2020/

#Load the STAR module; name may change in future so check each time with 'module avail' on command line
module load STAR/2.7.3a-IGB-gcc-8.2.0

#Set up variable to pull out base file names
line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" src/basenames.txt)

#Create output directory for the STAR results
# NOTE: the -p option will not only create parental directories if needed, but also will not
# throw an error or overwrite the directory if it already exists 
mkdir -p results/star

#run star
echo "Start STAR alignment"

STAR --runThreadN $SLURM_NTASKS \
     --genomeDir data/genome/STAR-2.7.3a_mouse-chr12_Index \
     --readFilesIn data/rawseq/${line}.fastq  \
     --sjdbGTFfile data/genome/mouse_chr12.gtf \
     --outFileNamePrefix results/star/${line}_ \
     --limitGenomeGenerateRAM 32000000000 \
     --outSAMtype BAM SortedByCoordinate \
     --outTmpDir /scratch/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}

echo "end STAR alignment"

echo "Start creating index of BAM file"

module load SAMtools/1.10-IGB-gcc-8.2.0

samtools index results/star/${line}_Aligned.sortedByCoord.out.bam

echo "Finished creating index of BAM file"
