#Exon counts created on compute1 as described here:
#/storage1/fs1/alberthkim/Active/data/hao_usp7_project/20220120_dexseq_analysis/exon_counts/not_aggregated
#/Volumes/alberthkim/Active/data/hao_usp7_project/20220120_dexseq_analysis/exon_counts/not_aggregated

#Set up parallel compute
#library("BiocParallel")
#BPPARAM = MulticoreParam(6)

#Read in the exon count data
inDirDiv6 = "/Users/mgriffit/Google Drive/My Drive/Manuscripts/USP7-AlbertKim/exon_counts/not_aggregated/div6"
countFilesDiv6 = list.files(inDirDiv6, pattern="exon_counts.tsv$", full.names=TRUE)
basename(countFilesDiv6)
inDirDiv7 = "/Users/mgriffit/Google Drive/My Drive/Manuscripts/USP7-AlbertKim/exon_counts/not_aggregated/div7"
countFilesDiv7 = list.files(inDirDiv7, pattern="exon_counts.tsv$", full.names=TRUE)
basename(countFilesDiv7)
inDirDiv8 = "/Users/mgriffit/Google Drive/My Drive/Manuscripts/USP7-AlbertKim/exon_counts/not_aggregated/div8"
countFilesDiv8 = list.files(inDirDiv8, pattern="exon_counts.tsv$", full.names=TRUE)
basename(countFilesDiv8)

#Master list of all files
countFilesDivN = c(countFilesDiv6, countFilesDiv7, countFilesDiv8)

#Read in exon annotations
inDir = "/Users/mgriffit/Google Drive/My Drive/Manuscripts/USP7-AlbertKim/exon_counts/not_aggregated"
flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
basename(flattenedFile)

#Set output dir
outDir = "/Users/mgriffit/Google Drive/My Drive/Manuscripts/USP7-AlbertKim/exon_counts/not_aggregated/results/"
setwd(outDir)

#Prepare a sample table for each batch
sampleTableDiv6 = data.frame(
  row.names = c("div6_cd_1.exon_counts.tsv", "div6_cd_2.exon_counts.tsv", "div6_cd_3.exon_counts.tsv", 
                "div6_wt_1.exon_counts.tsv", "div6_wt_2.exon_counts.tsv", "div6_wt_3.exon_counts.tsv"),
  condition = c("div6_cd", "div6_cd", "div6_cd",  
                "div6_wt", "div6_wt", "div6_wt"),
  libType = c( "paired-end", "paired-end", "paired-end", 
               "paired-end", "paired-end", "paired-end") )
sampleTableDiv6

sampleTableDiv7 = data.frame(
  row.names = c("div7_cd_1.exon_counts.tsv", "div7_cd_2.exon_counts.tsv", "div7_cd_3.exon_counts.tsv",
                "div7_wt_1.exon_counts.tsv", "div7_wt_2.exon_counts.tsv", "div7_wt_3.exon_counts.tsv"),
  condition = c("div7_cd", "div7_cd", "div7_cd",  
                "div7_wt", "div7_wt", "div7_wt"),
  libType = c( "paired-end", "paired-end", "paired-end",
               "paired-end", "paired-end", "paired-end") )
sampleTableDiv7

sampleTableDiv8 = data.frame(
  row.names = c("div8_cd_1.exon_counts.tsv", "div8_cd_2.exon_counts.tsv", "div8_cd_3.exon_counts.tsv",
                "div8_wt_1.exon_counts.tsv", "div8_wt_2.exon_counts.tsv", "div8_wt_3.exon_counts.tsv"),
  condition = c("div8_cd", "div8_cd", "div8_cd",  
                "div8_wt", "div8_wt", "div8_wt"),
  libType = c( "paired-end", "paired-end", "paired-end",
               "paired-end", "paired-end", "paired-end") )
sampleTableDiv8

#Prepare a sample table that combines all batches together
sampleTableDivN = data.frame(
  row.names = c("div6_cd_1.exon_counts.tsv", "div6_cd_2.exon_counts.tsv", "div6_cd_3.exon_counts.tsv", 
                "div6_wt_1.exon_counts.tsv", "div6_wt_2.exon_counts.tsv", "div6_wt_3.exon_counts.tsv",
                "div7_cd_1.exon_counts.tsv", "div7_cd_2.exon_counts.tsv", "div7_cd_3.exon_counts.tsv", 
                "div7_wt_1.exon_counts.tsv", "div7_wt_2.exon_counts.tsv", "div7_wt_3.exon_counts.tsv",
                "div8_cd_1.exon_counts.tsv", "div8_cd_2.exon_counts.tsv", "div8_cd_3.exon_counts.tsv", 
                "div8_wt_1.exon_counts.tsv", "div8_wt_2.exon_counts.tsv", "div8_wt_3.exon_counts.tsv"),
  condition = c("divN_cd", "divN_cd", "divN_cd",  
                "divN_wt", "divN_wt", "divN_wt",
                "divN_cd", "divN_cd", "divN_cd",  
                "divN_wt", "divN_wt", "divN_wt",
                "divN_cd", "divN_cd", "divN_cd",  
                "divN_wt", "divN_wt", "divN_wt"),
  libType = c( "paired-end", "paired-end", "paired-end", 
               "paired-end", "paired-end", "paired-end",
               "paired-end", "paired-end", "paired-end", 
               "paired-end", "paired-end", "paired-end",
               "paired-end", "paired-end", "paired-end", 
               "paired-end", "paired-end", "paired-end") )
sampleTableDivN


library( "DEXSeq" )
library( "openxlsx" )

#CREATE A FUNCTION THAT PROCESSES A PAIR THROUGH DEXSEQ
runDEXseq = function(countFiles, sampleTable, flattenedFile){
  dxd = DEXSeqDataSetFromHTSeq( countFiles, sampleData=sampleTable, 
                                design= ~ sample + exon + condition:exon,
                                flattenedfile=flattenedFile )

  #Inspect the example data
  colData(dxd)
  head( counts(dxd), 5 )

  #Note that the number of columns is 12, the first six (we have six samples) 
  # corresponding to the number of reads mapping to out exonic regions and the 
  # last six correspond to the sum of the counts mapping to the rest of the 
  # exons from the same gene on each sample.
  split( seq_len(ncol(dxd)), colData(dxd)$exon )
  head( featureCounts(dxd), 5 )
  head( rowRanges(dxd), 3 )
  sampleAnnotation( dxd )

  #Perform normalization
  dxd = estimateSizeFactors( dxd )

  #Perform dispersion estimation
  dxd = estimateDispersions( dxd )
  #dxd = estimateDispersions( dxd, BPPARAM=BPPARAM )
 
  #As a shrinkage diagnostic, the DEXSeqDataSet use the method plotDispEsts() 
  # that plots the per-exon dispersion estimates versus the mean normalised count,
  # the resulting fitted values and the a posteriori (shrinked) dispersion 
  # estimates (Figure 1). 
  plotDispEsts( dxd )

  #Test for differential expression 
  #dxd = testForDEU( dxd, BPPARAM=BPPARAM )
  dxd = testForDEU( dxd )

  #Get fold change values for differential exon usage
  dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
  #dxd = estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=BPPARAM)

  #Get the final results
  dxr = DEXSeqResults( dxd )
  
  return(dxr)
}

#dxr_div6 = runDEXseq(countFilesDiv6, sampleTableDiv6, flattenedFile)
#saveRDS(dxr_div6, file = "dxr_div6.rds")
dxr_div6 = readRDS(file = "dxr_div6.rds")

#dxr_div7 = runDEXseq(countFilesDiv7, sampleTableDiv7, flattenedFile)
#saveRDS(dxr_div7, file = "dxr_div7.rds")
dxr_div7 = readRDS(file = "dxr_div7.rds")

#dxr_div8 = runDEXseq(countFilesDiv8, sampleTableDiv8, flattenedFile)
#saveRDS(dxr_div8, file = "dxr_div8.rds")
dxr_div8 = readRDS(file = "dxr_div8.rds")

#dxr_divN = runDEXseq(countFilesDivN, sampleTableDivN, flattenedFile)
#saveRDS(dxr_divN, file = "dxr_divN.rds")
dxr_divN = readRDS(file = "dxr_divN.rds")

#Show column descriptions
mcols(dxr_div6)$description
mcols(dxr_div7)$description
mcols(dxr_div8)$description
mcols(dxr_divN)$description

#How many exonic regions are significant with a false discovery rate of 10%:
table ( dxr_div6$padj < 0.1 )  #7
table ( dxr_div7$padj < 0.1 )  #34
table ( dxr_div8$padj < 0.1 )  #53
table ( dxr_divN$padj < 0.1 )  #49

#We may also ask how many genes are affected:
table ( tapply( dxr_div6$padj < 0.1, dxr_div6$groupID, any ) ) #3
table ( tapply( dxr_div7$padj < 0.1, dxr_div7$groupID, any ) ) #30
table ( tapply( dxr_div8$padj < 0.1, dxr_div8$groupID, any ) ) #41
table ( tapply( dxr_divN$padj < 0.1, dxr_divN$groupID, any ) ) #34

#To see how the power to detect differential exon usage depends on the number 
# of reads that map to an exon, a so-called MA plot is useful, which plots the 
# logarithm of fold change versus average normalized count per exon and marks 
# by red color the exons which are considered significant; here, the exons with 
# an adjusted p values of less than 0.1 (There is of course nothing special 
# about the number 0.1, and you can specify other thresholds in the call to plotMA().
plotMA( dxr_div6, cex=0.8 )
plotMA( dxr_div7, cex=0.8 )
plotMA( dxr_div8, cex=0.8 )
plotMA( dxr_divN, cex=0.8 )

#Pull some examples with strong pvalues and fold changes
examples_div6 = dxr_div6[(which (dxr_div6$padj < 0.1 & abs(dxr_div6$log2fold_div6_wt_div6_cd) > 1.5)),]
examples_div7 = dxr_div7[(which (dxr_div7$padj < 0.1 & abs(dxr_div7$log2fold_div7_wt_div7_cd) > 1.5)),]
examples_div8 = dxr_div8[(which (dxr_div8$padj < 0.1 & abs(dxr_div8$log2fold_div8_wt_div8_cd) > 1.5)),]



#Usp7 positive control. Expected skipped junction: chr16:8706729-8708859
#Skipped exon coordinates: chr16:8707638-8707746
#This exon in DEXseq annotations looks like: ENSMUSG00000022710:E049, 16:8707638-8707746:-

#Try some visualizations starting with Usp7 itself
plotDEXSeq( dxr_div6, "ENSMUSG00000022710", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr_div7, "ENSMUSG00000022710", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr_div8, "ENSMUSG00000022710", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#Try some visualizations starting with Usp7 itself - show replicates
plotDEXSeq( dxr_div6, "ENSMUSG00000022710", expression=FALSE, norCounts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr_div7, "ENSMUSG00000022710", expression=FALSE, norCounts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr_div8, "ENSMUSG00000022710", expression=FALSE, norCounts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#Try some visualizations starting with Usp7 itself - show transcripts
plotDEXSeq( dxr_div6, "ENSMUSG00000022710", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr_div7, "ENSMUSG00000022710", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr_div8, "ENSMUSG00000022710", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#Try some visualizations starting with Usp7 itself - remove overall changes in expression from the plots
plotDEXSeq( dxr_div6, "ENSMUSG00000022710", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr_div7, "ENSMUSG00000022710", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr_div8, "ENSMUSG00000022710", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#Try some visualizations starting with Usp7 itself - Combine the options above
plotDEXSeq( dxr_div6, "ENSMUSG00000022710", displayTranscripts=TRUE, expression=FALSE, norCounts=TRUE, splicing=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr_div7, "ENSMUSG00000022710", displayTranscripts=TRUE, expression=FALSE, norCounts=TRUE, splicing=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr_div8, "ENSMUSG00000022710", displayTranscripts=TRUE, expression=FALSE, norCounts=TRUE, splicing=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


#Create visualizations for all genes passing a threshold
outDir2 = "/Users/mgriffit/Google Drive/My Drive/Manuscripts/USP7-AlbertKim/exon_counts/not_aggregated/results/div6"
setwd(outDir2)
x = dxr_div6[(which (dxr_div6$padj < 0.1)),]
c = length(x[,1])
o = order(x$padj)
y = x[o[1:c],]
z = unique(y$groupID)
c1 = length(z)
for (i in 1:c1){
  gene = z[i]
  name = paste(gene, ".pdf", sep="")
  pdf(file=name)
  plotDEXSeq( dxr_div6, gene, expression=FALSE, norCounts=TRUE, splicing=TRUE,
              legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
  dev.off()
}
#Write out the significant results
filtered_data = dxr_div6[(which (dxr_div6$padj < 0.1)),]
dim(filtered_data)
length(unique(filtered_data[,1]))
write.table(as.data.frame(filtered_data[,1:12]), file="DEXseq_Significant_Div6.tsv", quote=FALSE, sep="\t", row.names=FALSE)

outDir2 = "/Users/mgriffit/Google Drive/My Drive/Manuscripts/USP7-AlbertKim/exon_counts/not_aggregated/results/div7"
setwd(outDir2)
x = dxr_div7[(which (dxr_div7$padj < 0.1)),]
c = length(x[,1])
o = order(x$padj)
y = x[o[1:c],]
z = unique(y$groupID)
c1 = length(z)
for (i in 1:c1){
  gene = z[i]
  name = paste(gene, ".pdf", sep="")
  pdf(file=name)
  plotDEXSeq( dxr_div7, gene, expression=FALSE, norCounts=TRUE, splicing=TRUE,
              legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
  dev.off()
}
#Write out the significant results
filtered_data = dxr_div7[(which (dxr_div7$padj < 0.1)),]
dim(filtered_data)
length(unique(filtered_data[,1]))
write.table(as.data.frame(filtered_data[,1:12]), file="DEXseq_Significant_Div7.tsv", quote=FALSE, sep="\t", row.names=FALSE)

outDir2 = "/Users/mgriffit/Google Drive/My Drive/Manuscripts/USP7-AlbertKim/exon_counts/not_aggregated/results/div8"
setwd(outDir2)
x = dxr_div8[(which (dxr_div8$padj < 0.1)),]
c = length(x[,1])
o = order(x$padj)
y = x[o[1:c],]
z = unique(y$groupID)
c1 = length(z)
for (i in 1:c1){
  gene = z[i]
  name = paste(gene, ".pdf", sep="")
  pdf(file=name)
  plotDEXSeq( dxr_div8, gene, expression=FALSE, norCounts=TRUE, splicing=TRUE,
              legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
  dev.off()
}
#Write out the significant results
filtered_data = dxr_div8[(which (dxr_div8$padj < 0.1)),]
dim(filtered_data)
length(unique(filtered_data[,1]))
write.table(as.data.frame(filtered_data[,1:12]), file="DEXseq_Significant_Div8.tsv", quote=FALSE, sep="\t", row.names=FALSE)

#Create a web report
outDir3 = "/Users/mgriffit/Google Drive/My Drive/Manuscripts/USP7-AlbertKim/exon_counts/not_aggregated/results/div6/html/"
setwd(outDir3)
DEXSeqHTML( dxr_div6, FDR=0.01, color=c("#FF000080", "#0000FF80") )
setwd(outDir)
