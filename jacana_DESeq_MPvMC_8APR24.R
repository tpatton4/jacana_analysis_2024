library(DESeq2)
#read in gene count data
setwd("C:/Users/tmp15/Downloads/R/jacana")
counts_genes <- as.matrix(read.csv("jacana_POATnAcounts_14FEB24.csv",row.names="Gene.stable.ID"))
counts <- subset(counts_genes,select=-c(Gene.name)) #remove gene names column for DESeq steps
nom_cols <- grep("NOM", colnames(counts), value = FALSE)
counts <- counts[, nom_cols]
counts <- as.matrix(counts)
mode(counts) <- "integer"

treat <- read.csv("jacana_allsamples_sexStageInfo.csv",header=TRUE,row.names=1)
nom_treat <- grep("NOM", rownames(treat), value= FALSE)
treat <- treat[nom_treat, ]
treat$group <- paste(treat$Brain.Region,treat$Stage,sep="_")
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

#Make MPvMC comparisons in each brain region
resPOAMPvMC <- results(dds2, contrast=c("group", "POA_Parenting", "POA_Courting"))
resTnAMPvMC <- results(dds2, contrast=c("group", "TnA_Parenting", "TnA_Courting"))

summary(resPOAMPvMC)
summary(resTnAMPvMC)

allgenesPOAMPvMC <- merge(as.matrix(resPOAMPvMC), counts_genes, by="row.names") #add gene names and counts
allgenesTnAMPvMC <- merge(as.matrix(resTnAMPvMC), counts_genes, by="row.names")

#adjusted pval less than 0.05
siggenesPOAMPvMC <- as.matrix(subset(resPOAMPvMC,padj <0.05))
siggenesTnAMPvMC <- as.matrix(subset(resTnAMPvMC,padj <0.05))

#create file with counts so all the info for significant DEGs is together
sigcountsPOAMPvMC <- merge(siggenesPOAMPvMC, counts_genes, by="row.names", all=F) #merge the matrices by gene IDs, all=F so we don't have a bunch of NA lines for non-significant genes
sigcountsTnAMPvMC <- merge(siggenesTnAMPvMC, counts_genes, by="row.names", all=F)

#write this file out
write.csv(sigcountsPOAMPvMC, "jacanaPOA_MPvMC_DESeq2Out_8APR24.csv")
write.csv(sigcountsTnAMPvMC, "jacanaTnA_MPvMC_DESeq2Out_8APR24.csv")

