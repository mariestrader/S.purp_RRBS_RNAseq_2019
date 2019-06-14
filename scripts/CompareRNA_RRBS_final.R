setwd("/Users/mariestrader/Dropbox/RRBS_Loaded")
############
library(methylKit)
library(genomation)
library(GenomicRanges)
library(reshape)
library(ggplot2)
library(vegan)
library('adegenet')
library('VennDiagram')
library(dplyr)
library(stringr)
library(pheatmap)

### load RNAseq data ###
d1=load("RNAseq/LoadedRNA_rld&res_gene_bm3_JUNE2019.Rdata")
# or if considering all genes: d1=load("RNAseq/LoadedRNA_rld&res_gene_all_JUNE2019.Rdata")
d2=load("RNAseq/rldpvals_all_bm3_JUNE2019.Rdata")

### load differential methylation analysis data ###
d3=load("bismarkOutput/Ensembl/methylKit_Spurp_ensembl_Maternal_CpGs.Rdata")
myDiff25pM=getData(myDiff25p)
d4=load("bismarkOutput/Ensembl/methylKit_Spurp_ensembl_Developmental_CpGs.Rdata")
myDiff25pD=getData(myDiff25p)

### load GS-CpG data ###
d5=load("bismarkOutput/Ensembl/GS-CpGs.Rdata")

### adjust gene names in rld file ###
row.names(rld.df)<-gsub("-tr","",row.names(rld.df))
row.names(rld.df)<-sub("transcript:","", row.names(rld.df))

### get means of gene counts across samples ###
rld.df=as.data.frame(rld.df)
rld.df$means=apply(rld.df[,1:12], 1, mean)
rld.df$gene=row.names(rld.df)

### rename gene annotation columns ###
dat_med$gene=dat_med$V41
dat5_med$gene=dat5_med$V41

dat_med_m=merge(dat_med, rld.df, by="gene") #6368 genes (BM3), 9204 genes (all genes)
dat5_med_m=merge(dat5_med, rld.df, by="gene") #3359 genes (BM3), 5024 (all genes)

# Figure S5 - plotted using all genes
par(mfrow=c(2,1))
plot(density(dat5_med_m$percent_methylated, bw=0.05), col='black', main='Percent Methylation in Genes')
plot(means ~ percent_methylated, dat5_med_m,pch=16,cex=0.5,col=rgb(0,0,0,alpha=0.2), xlab="Percent Methylated")
lm1 = lm(means ~ percent_methylated, dat5_med_m)
abline(lm1, col = 'turquoise4')
summary(lm1) #Adjusted R-squared:  0.07448,  p-value: < 2.2e-16 (BM3), Adjusted R-squared:  0.1836 (all genes)

#########
#### grab gene assignments for each DMCpG , use dataset with all genes ###
id_gene=as.data.frame(meth.annot.perc.meth[-c(2:13,15:18)])
colnames(id_gene)=c("id","gene")

myDiff25pM=as.data.frame(myDiff25pM)
myDiff25pM$id <- paste(myDiff25pM$chr,myDiff25pM$start, myDiff25pM$end, sep=".")
maternal=myDiff25pM$id
maternal_gene=merge(myDiff25pM, id_gene, by="id")

length(maternal_gene$gene) #215 sig maternal CpGs are GS-CpGs
length(unique(maternal_gene$gene)) #136 genes have sig maternal CpG within a gene

sig_maternal=merge(dat_med_m, maternal_gene, by="gene")

length(unique(sig_maternal$gene)) #118 genes have DEG data

####
myDiff25pD=as.data.frame(myDiff25pD)
myDiff25pD$id <- paste(myDiff25pD$chr,myDiff25pD$start, myDiff25pD$end, sep=".")
dev=myDiff25pD$id
dev_gene=merge(myDiff25pD, id_gene, by="id")

length(dev_gene$gene) #82 sig develop. CpGs are GS-CpGs. 
length(unique(dev_gene$gene)) #49 unique genes have sig developmental CpGs within a gene. 

sig_dev=merge(dat_med_m, dev_gene, by="gene")
length(unique(sig_dev$gene)) #40 genes have DEG data


# combine the two sig CpG datasets to plot all DM CpGs, there should be some overlap between them but will just plot the same point twice
maternal_gene$treat="M"
dev_gene$treat="D"
sig_CpG=rbind(maternal_gene, dev_gene)
dim(sig_CpG) #297 total GS-CpGs that are also DMCpGs (including those overlapping bewteen treatments)
length(unique(sig_CpG$gene)) #177 genes with DMCpGs
length(unique(sig_CpG$id)) #258 total DMCpGs that are also GS-CpGs

sigCpG.genes=unique(sig_CpG$gene)

annot=read.table("gene_info_table_header.txt",header=T,sep="\t",quote=NULL)
annot=as.data.frame(annot[-c(3:6, 11:19) ],)
genes = sig_CpG$gene
genes2annot = match(genes,annot$spu_id)
sig_CpG_annot=data.frame(GeneID=sig_CpG$gene, sig_CpG, annot[genes2annot,])
sig_CpG_annot = sig_CpG_annot[order(sig_CpG_annot$qvalue),]

save(sig_maternal,sig_dev,sig_CpG,sigCpG.genes, sig_CpG_annot, file="sigCpGdata.Rdata")
write.table(sig_CpG_annot,file="SigCpGs_annots.tab",sep="\t",quote=F, row.names=F)

######## get DEG information to compare with GS-CpGs and DMCpGs - use BM3 datatset
resM=as.data.frame(resM_NU)
resD=as.data.frame(resD_NU)

row.names(resM)<-gsub("-tr","",row.names(resM))
row.names(resM)<-sub("transcript:","", row.names(resM))

row.names(resD)<-gsub("-tr","",row.names(resD))
row.names(resD)<-sub("transcript:","", row.names(resD))

resM$gene=row.names(resM)
resD$gene=row.names(resD)

resM_meth=merge(resM, dat5_med_m, by="gene") #3359 genes have >5 CpGs/gene represented in the RRBS data, 9204 when considering all genes
resM_meth_sig=merge(resM_meth, maternal_gene, by="gene") #134 genes, 215 when considering all genes
resM_meth_sigG=resM_meth[resM_meth$padj<=0.05 & !is.na(resM_meth$padj),] #439 genes

resD_meth=merge(resD, dat5_med_m, by="gene") #3359 genes have >5 CpGs/gene represented in the RRBS data
resD_meth_sig=merge(resD_meth, dev_gene, by="gene") #50 genes, 82 when considering all genes
resD_meth_sigG=resD_meth[resD_meth$padj<=0.05 & !is.na(resD_meth$padj),] #941 genes, 1608 when considering all genes


##### plots for Figure 5
par(mfrow=c(2,2))
plot(density(resD_meth$percent_methylated, bw=0.05), col='black', ylim=c(0,7), xlim=c(0,1.1), main="DEGs developemental")
lines(density(resD_meth_sigG$percent_methylated, bw=0.05), col='turquoise4')
lines(density(resD_meth_sig$percent_methylated, bw=0.05), col='purple')

plot(density(resM_meth$percent_methylated, bw=0.05), col='black', ylim=c(0,7), xlim=c(0,1.1), main="DEGs maternal")
lines(density(resM_meth_sigG$percent_methylated, bw=0.05), col='coral2')
lines(density(resM_meth_sig$percent_methylated, bw=0.05), col='purple')

plot(log2FoldChange ~ percent_methylated, resD_meth, xlab="Percent Methylation", ylab="log2FoldChange Developmental",pch=16,cex=0.5,col="black",ylim=c(-6,6), main="DEGs developmental")
points(log2FoldChange ~ percent_methylated, resD_meth_sigG, col = 'turquoise4', pch=16,cex=0.5)
points(log2FoldChange ~ percent_methylated, resD_meth_sig, col = 'purple', pch=16,cex=0.5)
lm1 = lm(log2FoldChange ~ percent_methylated, resD_meth_sigG)
abline(lm1, col = 'turquoise4')
summary(lm1) #Adjusted R-squared:  0.1152 p-value: < 2.2e-16
lm2 = lm(log2FoldChange ~ percent_methylated, resD_meth)
abline(lm2, col = 'black')
summary(lm2) #Adjusted R-squared:  0.04867 p-value: < 2.2e-16
lm3 = lm(log2FoldChange ~ percent_methylated, resD_meth_sig)
abline(lm3, col = 'purple')
summary(lm3) #Adjusted R-squared:  0.1719 p-value: 0.001614

plot(log2FoldChange ~ percent_methylated, resM_meth, xlab="Percent Methylation", ylab="log2FoldChange Maternal",pch=16,cex=0.5,col="black",ylim=c(-6,6), main="DEGs maternal")
points(log2FoldChange ~ percent_methylated, resM_meth_sigG, col = 'coral2', pch=16,cex=0.5)
points(log2FoldChange ~ percent_methylated, resM_meth_sig, col = 'purple', pch=16,cex=0.5)
lm1 = lm(log2FoldChange ~ percent_methylated, resM_meth_sigG)
abline(lm1, col = 'coral2')
summary(lm1) #Adjusted R-squared:  -0.002091 p-value: 0.7696
lm2 = lm(log2FoldChange ~ percent_methylated, resM_meth)
abline(lm2, col = 'black') 
summary(lm2) #Adjusted R-squared:  0.0002878 p-value: 0.1609
lm3 = lm(log2FoldChange ~ percent_methylated, resM_meth_sig)
abline(lm3, col = 'purple')
summary(lm3) #Adjusted R-squared:  -0.007346 p-value: 0.8624

###### make heatmaps with DEGs and DMCpGs

setwd("/Users/mariestrader/Dropbox/RRBS_Loaded/RNAseq/WGCNA")

datOutput <- read.table("Loaded_annotatedNetworkAnalysisResult_bm10_merged0.2_signed_methGenes.txt", header=T, sep="\t")

rld.df$spu_id<-rld.df$gene

combo=merge(rld.df, datOutput, by="spu_id")
dim(combo) #124  59

resM_NU=as.data.frame(resM_NU)
row.names(resM_NU)<-gsub("-tr","",row.names(resM_NU))
row.names(resM_NU)<-sub("transcript:","", row.names(resM_NU))
resM_NU$spu_id<-row.names(resM_NU)

comboM=merge(combo, resM_NU, by="spu_id")
dim(comboM) #124  57

comboM=comboM[comboM$padj<= 0.05, ]
dim(comboM) #15  7

row.names(comboM) <- make.names(paste(comboM$spu_id,comboM$moduleColors, sep="."), unique=T)
#row.names(comboM)<-comboM$spu_id
comboM=comboM[2:13]
means=apply(comboM,1,mean) # means of rows
explc=comboM-means

heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=0.8)(100)
pdf("SigDifMeth_heatmap_maternalDEGs_spu.pdf",height=2.5,width=4.5)
pheatmap(explc,color=heat.colors,cluster_cols=F,border_color=NA,clustering_distance_rows="correlation")
dev.off()

##### for dev. 
resD_NU=as.data.frame(resD_NU)
row.names(resD_NU)<-gsub("-tr","",row.names(resD_NU))
row.names(resD_NU)<-sub("transcript:","", row.names(resD_NU))
resD_NU$spu_id<-row.names(resD_NU)

comboD=merge(combo, resD_NU, by="spu_id")
dim(comboD) #124  57

comboD=comboD[comboD$padj<= 0.05, ]
dim(comboD) #31  7

row.names(comboD) <- make.names(paste(comboD$common_name,comboD$moduleColors, sep="."), unique=T)
comboD=comboD[2:13]
means=apply(comboD,1,mean) # means of rows
explc=comboD-means

heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=0.8)(100)
pdf("SigDifMeth_heatmap_devDEGs_spu.pdf",height=5,width=4.5)
pheatmap(explc,color=heat.colors,cluster_cols=F,border_color=NA,clustering_distance_rows="correlation")
dev.off()

meth.data=load("sigCpGdata.Rdata")
length(unique(sig_CpG$gene)) #177

sigMethSum <- 
  # split data by gene: 
  group_by(sig_CpG, gene) %>%
  # make summary table. to get total cpg sites for each gene just take the length of the column: 
  summarize(total_sig_cpg = length(id))
names(sigMethSum)[1]<-"spu_id"

megaData=merge(combo, sigMethSum, by="spu_id")
megaData=as.data.frame(megaData[-c(17:19,26:34) ])
write.table(megaData,"sigCpG_ALLinfor.txt",quote=F,row.names=F, sep="\t")

##### output for fishers test
rld.ge=load("RNAseq/LoadedRNA_rld&res_gene_ALL_JUNE2019.Rdata")
rld.df=as.data.frame(rld.df)
row.names(rld.df)<-gsub("-tr","",row.names(rld.df))
row.names(rld.df)<-sub("transcript:","", row.names(rld.df))
rld.df$gene=row.names(rld.df)

methTF <- sapply(rld.df$gene, function(x) { # for each value in df1$x
  j <- which(sig_CpG$gene == x) # check if df2$y has a match
  ifelse(length(j) > 0, T, F) # if there is, give the location in the vector
}) 
methTF.df=as.data.frame(methTF, row.names=names(methTF))
methT=methTF.df[methTF.df$methTF==TRUE, ]
methTF.df$spu_id=row.names(methTF.df)

methTF.df$go[methTF.df$methTF==T]<-1
methTF.df$go[methTF.df$methTF==F]<-0 
methTF.df$go=as.character(methTF.df$go) 
summary(methTF.df$go)
#    0     1 
#30107   177

methTF.df=methTF.df[2:3]

write.csv(methTF.df,file="RNAseq/GO/sigCpG_fishers.csv",quote=F, row.names=F)

