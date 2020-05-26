# Mayo-Illinois RNASeq Lab
Run the following commands on the Biocluster classroom partition to complete the lab. 

## Set up directory. This takes less than a minute.
$ cd ~/

$ cp /home/classroom/hpcbio/mayo-rnaseq/mouse-rnaseq-2020/src/Mayo-RNASeq/prep-directory.sh ./

$ sbatch prep-directory.sh 

## Create a STAR index of the mouse genome (chromosome 12 only). This takes ~3-5 minutes to run.
$ cd mouse-rnaseq-2020/src/

$ sbatch STAR-index-mouse-genome.sh

## Run STAR alignments on 4 mouse samples. This takes ~3 minutes per sample, but only runs 2 samples at once, so could take 6-12 minutes total depending on how long they're waiting for resources to open up.
$ sbatch STAR-alignment.sh

## Run featureCounts on 4 alignment files to get gene counts. This takes ~X minutes to run.
$ sbatch featureCounts.sh

## ADD EDGER HERE
