library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

#preparing the samples

sampleinfo <- read.table("d:/LATEST_STAGES_count.txt")
colnames(sampleinfo)<-(sampleinfo[1,])
head(sampleinfo)
rownames(sampleinfo)<-(sampleinfo[,1])
colnames(sampleinfo)<-sub("/mnt/d/featurecounts/LATEST_STAGES/","",colnames(sampleinfo))
colnames(sampleinfo)<-sub("_sequence_aligned_sorted.bam","",colnames(sampleinfo))
sampleinfo<-(sampleinfo[-1,-1])
sampleinfo<-sampleinfo[,c(5:15)]
mycount<-as.data.frame(sapply(sampleinfo,as.numeric))
mycount2<-(mycount[,-1])
rownames(mycount2)<-rownames(sampleinfo)

#method1
mycount2<-calcNormFactors(mycount2)
myCPM <- cpm(mycount2)
head(myCPM)
thresh <- myCPM > 0.5
head(thresh)
table(rowSums(thresh))
keep <- rowSums(thresh) >= 2
counts.keep <- mycount2[keep,]
summary(keep)
dim(counts.keep)

plot(myCPM[,1],mycount2[,1])
plot(myCPM[,1],mycount2[,1],ylim=c(0,50),xlim=c(0,3))
abline(v=0.5)

library("compGenomRData")
counts_file <- system.file("/mnt/d/featurecounts/LATEST_STAGES/LATEST_STAGES_count.txt",
                           package = "compGenomRData")

#method2

cpm <- cpm(mycount2, log = FALSE, normalized.lib.sizes=TRUE)
           
Patient <- factor(c("C19","C19","S35","S35","S42","S42","S10","S10","S47","S47"))
Tissue <- factor(c("T","N","T","N","T","N","T","N","T","N"))

dgeObj <- DGEList(counts.keep)
dgeObj
names(dgeObj)
dgeObj$samples


dgeObj$samples$lib.size

barplot(dgeObj$samples$lib.size, names=colnames(dgeObj), las=2)
logcounts <- cpm(dgeObj,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

plotMDS(dgeObj)
par(mfrow=c(1,2))
levels(sampleinfo$CellType)
