# Mayo-Illinois RNASeq Lab
Run the following commands on the Biocluster classroom partition to complete the lab. 

## Set up directory. This takes less than a minute.
$ cd ~/

$ cp /home/classroom/hpcbio/mayo-rnaseq/mouse-rnaseq-2020/src/Mayo-RNASeq/prep-directory.sh ./

$ sbatch prep-directory.sh 

## Create a STAR index of the mouse genome (chromosome 12 only). This takes ~3-5 minutes to run.
$ cd mouse-rnaseq-2020/src/

$ sbatch STAR-index-mouse-genome.sh

## Run STAR alignments on 4 mouse samples. 
### This takes ~1 minute per sample, but only runs 2 samples at once, so could take 2-4 minutes total depending on how long they're waiting for resources to open up.
$ sbatch STAR-alignment.sh

## Run featureCounts on 4 alignment files to get gene counts. 
### This takes ~1 minute per sample to run. Runs all 4 samples at once, so could take 1-4 minutes total depending on how long they're waiting for resources to open up.
$ sbatch featureCounts.sh

## Run multiqc to collate STAR and featureCounts numbers. This takes ~1 min to run.
$ sbatch multiqc_summary.sh

*OUTPUT* in results/ directory: 
1. multiqc_report.html 
2. ReadFatePlot.jpeg 
3. Targets_Final.txt

Download these to VMs to inspect

## Run edgeR statistical analysis. This takes ~30 seconds to run
$ sbatch edgeR.sh  

*OUTPUT* in results/edgeR/ directory:
1. MDSclustering.jpeg  
2. NumSigGenes_FDR0.05.csv  
3. RawCounts.txt  
4. t8_vs_t0_AllResults.txt  
5. t8_vs_t0_MeanDifferencePlot.jpeg

Download these to VMs to inspect
