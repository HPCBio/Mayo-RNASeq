#Script to analysis 2020 Mayo CompGen RNA-Seq lab
#Adapted from https://f1000research.com/articles/5-1438

options(stringsAsFactors = FALSE)

library(edgeR)
library(org.Mm.eg.db)

#setwd("C:/Keck/Workshops_spring2020/CompGen")
setwd("~/mouse-rnaseq-2020/")

#Read in Targets_Final.txt file, make some additional columns and re-order

targets <- readTargets("results/Targets_Final.txt")

#Read in the counts

d <- readDGE(targets, path = "results/featureCounts/", columns = c(1,7),
             labels = targets$Sample, comment.char = "#", header = TRUE)

#add in annotations

d$genes <- data.frame(ENTREZID = gsub("GeneID:","", rownames(d$counts)))

temp <- select(org.Mm.eg.db, keys = d$genes$ENTREZID, ketype = "ENTREZID",
               columns = c("SYMBOL","GENENAME"))

d$genes <- dplyr::left_join(d$genes, temp)
rownames(d$genes) <- d$genes$ENTREZID
rownames(d$counts) <- d$genes$ENTREZID


#Filter out genes without 0.5 cpm in at least 2 samples

keep <- rowSums(cpm(d) > 0.5) >= 2

#Write out raw counts and keep information

write.table(cbind(d$genes, keep = keep, d$counts), file = "results/edgeR/RawCounts.txt",
            row.names = FALSE, sep = "\t")


d <- d[keep, , keep.lib.sizes=FALSE]

#Do TMM normalization

d <- calcNormFactors(d)

#Do MDS clustering

jpeg("results/edgeR/MDSclustering.jpeg", width = 5, height = 5, units = "in", res = 300, quality = 100, type = "cairo")
plotMDS(d, col=rep(1:2,2))
dev.off()


#Make design matrix

design <- model.matrix(~d$samples$Time)
colnames(design)[2] <- "t8_vs_t0"


#Estimate dispersion and fit edgeR-quasi model

d <- estimateDisp(d, design, robust=TRUE)
fit <- glmQLFit(d, design, robust=TRUE)

# Get results

res <- glmQLFTest(fit, coef = 2)

is.de <- decideTestsDGE(res)

write.csv(summary(is.de), file = "results/edgeR/NumSigGenes_FDR0.05.csv")

jpeg("results/edgeR/t8_vs_t0_MeanDifferencePlot.jpeg", width = 5, height = 5, units = "in",
    res = 300, quality = 100, type = "cairo")
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")
dev.off()

#write out all results

temp <- topTags(res, n = Inf)$table

#back-translate fold-change

temp$F <- 2^abs(temp$logFC) * sign(temp$logFC)
names(temp)[6] <- "FoldChange"

#add in normalized values

logCPM <- data.frame(cpm(d, log = TRUE))

temp <- dplyr::left_join(temp, tibble::rownames_to_column(logCPM, var = "ENTREZID"))

write.table(temp, file = "results/edgeR/t8_vs_t0_AllResults.txt",
            row.names = FALSE, sep = "\t")

print("Finished edgeR analysis")

q()

