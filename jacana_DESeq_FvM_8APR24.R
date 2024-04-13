library(DESeq2)
#read in gene count data
setwd("C:/Users/tmp15/Downloads/R/jacana")
counts_genes <- as.matrix(read.csv("jacana_POATnAcounts_14FEB24.csv",row.names="Gene.stable.ID"))
counts <- subset(counts_genes,select=-c(Gene.name)) #remove gene names column for DESeq steps
counts <- as.matrix(counts)
mode(counts) <- "integer"

treat <- read.csv("jacana_allsamples_sexStageInfo.csv",header=TRUE,row.names=1)
treat$group <- paste(treat$Brain.Region,treat$Sex,sep="_")
#treat <- subset(treat,select=-c(Sex,Stage))

all(rownames(treat) %in% colnames(counts))
#true - all sample IDs are in 
all(rownames(treat) == colnames(counts))
#true - all sample IDs are in the correct order

#make a DESeq object
dds <- DESeqDataSetFromMatrix(countData=counts, colData=treat, design= ~ group)
dds
yes = rowSums(counts > 10) >= 5
dds2 <- dds[yes,]
#make F the reference
#dds2$Sex <-relevel(dds$SexStage, ref="Female")
dds2<-DESeq(dds2)

#Make FvM comparisons in each brain region
resPOAFvM <- results(dds2, contrast=c("group", "POA_Female", "POA_Male"))
resTnAFvM <- results(dds2, contrast=c("group", "TnA_Female", "TnA_Male"))

summary(resPOAFvM)
summary(resTnAFvM)

allgenesPOAFvM <- merge(as.matrix(resPOAFvM), counts_genes, by="row.names") #add gene names and counts
allgenesTnAFvM <- merge(as.matrix(resTnAFvM), counts_genes, by="row.names")

#adjusted pval less than 0.05
siggenesPOAFvM <- as.matrix(subset(resPOAFvM,padj <0.05))
siggenesTnAFvM <- as.matrix(subset(resTnAFvM,padj <0.05))

#create file with counts so all the info for significant DEGs is together
sigcountsPOAFvM <- merge(siggenesPOAFvM, counts_genes, by="row.names", all=F) #merge the matrices by gene IDs, all=F so we don't have a bunch of NA lines for non-significant genes
sigcountsTnAFvM <- merge(siggenesTnAFvM, counts_genes, by="row.names", all=F)

#write this file out
write.csv(sigcountsPOAFvM, "jacanaPOA_FvM_DESeq2Out_8APR24.csv")
write.csv(sigcountsTnAFvM, "jacanaTnA_FvM_DESeq2Out_8APR24.csv")


####Volcano Plots####
#POA FvM
with(allgenesPOAFvM, plot(log2FoldChange, -log10(padj), pch=19, main="Females vs Males POA", cex=0.7, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value), xlim=c(-3,13), ylim=c(-1,300)))

#Add lines for cut-offs: logFC>0.5 and p-value cut-off at padj<0.05
abline(h=-log10(0.05), col="black", lty=3, lwd=1)
abline(v=-0.5, col="black", lty=3, lwd=1)
abline(v=0.5, col="black", lty=3, lwd=1)

#color significant genes based on whether they were up- or down-regulated
with(subset(allgenesPOAFvM, padj<0.05 & log2FoldChange< -0.5), points(log2FoldChange, -log10(padj), pch=19, col="blue", cex=0.7))
with(subset(allgenesPOAFvM, padj<0.05 & log2FoldChange>0.5), points(log2FoldChange, -log10(padj), pch=19, col="red", cex=0.7))

#Add genes names to the significant genes using the code below. Try adjusting how many gene names are shown.
cutoff=sort(allgenesPOAFvM$padj)[10] #selects the top 10 smallest p values
sign.genes=which(allgenesPOAFvM$padj <= cutoff)
text(x=allgenesPOAFvM$log2FoldChange[sign.genes] , y=-log10(allgenesPOAFvM$padj[sign.genes]), label=allgenesPOAFvM$Gene.name[sign.genes], cex=0.5)
