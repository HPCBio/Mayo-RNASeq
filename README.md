# Mayo-RNASeq
Scripts for RNA-Seq module of Computational Genomics Course

Run the following commands on the Biocluster classroom partition to complete the lab. 

$ cd ~/

$ cp /home/classroom/hpcbio/mayo-rnaseq/mouse-rnaseq-2020/src/Mayo-RNASeq/prep-directory.sh ./

$ sbatch prep-directory.sh

$ sbatch STAR-index-mouse-genome.sh

$ sbatch STAR-alignment.sh

$ sbatch featureCounts.sh

#ADD EDGER AFTER THIS
