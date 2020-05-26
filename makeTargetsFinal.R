# R script to take input sample names + number of reads from sequencing center
# and collate reads lost to other fates from multiQC reports

# Needs:
# 1. File with samples names and number of reads from sequencing report. Should be 
#    called "Targets0.txt"; samples names must be in column called "Sample", number of reads in column called "NumReads".
#    Optionally, if the sequencing center filtered out any reads, these can be in a column called "ReadsLost". 
#    Sample names are assumed to be the SHORT version - not including the 
#    index or anything after. Any number of additional columns are allowed.
#    However, column names "files", "QCfiltered", "unmapped", "multimapped",
#    "not.in.gene", "ambiguous", "in.a.gene" and "Total" are reserved.
# 2. MultiQC run on STAR alignments and featureCounts; these can be in any sub-directory
# 3. Folder of count files from featureCounts (i.e.,"*.txt" and matching
#    "*txt.summary" files); these can be in any subdirectory

# How to run (simple):
# Login to worker node. Swtich to directory containing the Targets0.txt file and multiQC/featureCounts/STAR results
# load R module:
# $ module load R/3.4.1-IGB-gcc-4.9.4 
# $ Rscript pathToFolder/makeTargetsFinal.R

# Outputs:
# 1. Targets_ReadFates.txt file with additional categories of lost reads
#     and file names of the featureCounts .txt files
# 2. JPEG of read fate plot: ReadFatePlot.jpeg 

#Read in targets file 

targets <- read.delim(dir(pattern = "Targets0.txt", recursive = TRUE, ignore.case = TRUE), 
as.is = TRUE)


#Read in featureCount file names; remove directory appendation and .summary at end

fc.names <- dir(pattern = ".txt.summary", recursive = TRUE)
temp <- strsplit(fc.names, "/", fc.names)
fc.names <-  gsub(".summary", "", sapply(temp, function(x) x[length(x)]),fixed = TRUE)



#Read in multQC output for STAR and featureCounts

star <- read.delim(dir(pattern = "multiqc_star.txt", recursive = TRUE), as.is = TRUE)
fc <- read.delim(dir(pattern = "multiqc_featureCounts.txt", recursive = TRUE), as.is = TRUE)

#add in sample names

names(fc.names) <- targets$Sample 

rownames(star) <- star$Sample

rownames(fc) <- fc$Sample


#if no ReadsLost in targets, add all 0

if(!"ReadsLost" %in% names(targets))
  targets$ReadsLost <- 0

#add in featureCounts file names

if(sum(targets$Sample %in% names(fc.names)) == nrow(targets))
  targets$files <- fc.names[targets$Sample]


 #allow for other columns
 n.targCols <- ncol(targets)

#check and adjust orders of names

if(sum(targets$Sample %in% rownames(star)) == nrow(targets))
  star <- star[targets$Sample,]
if(sum(targets$Sample %in% rownames(fc)) == nrow(targets)) {
  fc <- fc[targets$Sample,]
#  targets$QCfiltered <- targets$ReadsLost + targets$NumReads - star$total_reads
  targets$unmapped <- star$unmapped_other + star$unmapped_tooshort
  targets$multimapped <- star$multimapped_toomany + star$multimapped
  targets$not.in.gene <- fc$Unassigned_NoFeatures
  targets$ambiguous <- fc$Unassigned_Ambiguity
  targets$in.a.gene <- fc$Assigned
  targets$Total <- rowSums(targets[,(n.targCols+1):(n.targCols+5)])

  write.table(targets, file = "Targets_Final.txt", row.names = FALSE, sep = "\t")

  read.fate <- targets[,(n.targCols+1):(n.targCols+5)] / targets$Total * 100
# make read fate plot on biocluster2
    jpeg("ReadFatePlot.jpeg", width = 10, height = 6, units = "in", res = 300, quality = 100, type = "cairo")
        barplot( t(read.fate), beside = T, col = 1:5, las = 2, ylim = c( 0,110), legend.text = F,
             cex.names = 0.8, names.arg = targets$Sample)
        abline(h = seq(0,100, by = 10), col = "gray")
        barplot( t(read.fate), beside = T, col = 1:5, las = 2, ylim = c( 0,110), legend.text = F,
             ylab = "Percent of total reads", main = "Read fates per sample", cex.names = 0.8,
             names.arg = targets$Sample, add = TRUE)
        legend ("topright", legend = colnames( read.fate), col = 1:5, fill = 1:5, ncol = 5 )
    dev.off()
} else {
    write.table("Sample IDs do not match between Targets.txt file and STAR and featureCount
    multiQC reports", file = "error_log_ReadFates.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
  }
