#This script uses DESeq2 to identify DEGs associated with maternal and developmental conditioning in purple sea urchins. R V3.4.2

setwd("/Users/mariestrader/Dropbox/RRBS_Loaded/RNAseq")

library('DESeq2')
library('arrayQualityMetrics')
library('vegan')
library('rgl')
library('ape')
library('pheatmap')
library('dplyr')
library('adegenet')
library('VennDiagram')

gcountsALL=read.table("geneCounts_02122019.txt", header=T) 

gcounts=as.data.frame(gcountsALL[-c(2:6, 16) ]) #removing top rows with general count numbers
row.names(gcounts)<-gcounts$Geneid
gcounts$Geneid=NULL

length(gcounts[,1]) #30284
dim(gcounts) #30284    12

colnames(gcounts)=c("NN1","NN2","NN3","NU1","NU2","NU3","UN1","UN2","UN3","UU1","UU2","UU3")

summary(colSums(gcounts))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#6376636 7520221 7723652 7896547 8558910 8948115   

totalCounts=colSums(gcounts)
totalCounts
#    NN1     NN2     NN3     NU1     NU2     NU3     UN1     UN2     UN3     UU1     UU2     UU3 
#8555156 8577700 8948115 8545753 7170372 7825609 7570387 7621696 7541956 8570174 7455015 6376636 

min(totalCounts) #6376636
mean(totalCounts) #7896547
max(totalCounts)  #8948115

### REMOVE GENES WITH LOW MEAN COUNTS ###

mns = apply(gcounts, 1, mean)

gcounts=gcounts[mns>3,] #get rid of genes that show little or no expression
#gcounts=gcounts[mns>10,] #only highly expressed genes for WGCNA

table(mns > 3)
#FALSE  TRUE 
#12269 18015 
#### for means >10 for WGCNA
#FALSE  TRUE 
#15142 15142 

dim(gcounts) #18015    12

### BUILD A DATAFRAME ASSOCIATING SAMPLE NAMESWITH TREATMENT CONDITIONS ###

colData <- data.frame(row.names= colnames(gcounts), bucket= c("NN1","NN2","NN3","NU1","NU2","NU3","UN1","UN2","UN3","UU1","UU2","UU3"))
colData$treat = substr(colData$bucket, 1,2)
colData$treat_maternal = substr(colData$bucket, 1,1)
colData$treat_dev = substr(colData$bucket, 2,2)

### OUTLIERS ###
dds<-DESeqDataSetFromMatrix(countData=gcounts, colData=colData, design= ~ treat)
vsd=varianceStabilizingTransformation(dds, blind=T)
e=ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("treat"), force=T, outdir= "report_for_genes_Loaded")
# no outliers

### WALD TEST - FULL MODEL ###

dds<-DESeqDataSetFromMatrix(gcounts,

	colData = colData, 

	design = formula(~ treat_maternal + treat_dev))

rld=rlog(dds)
rld.df = assay(rld)

### Output rdl and trait data for WGCNA analysis, see WGCNA script for details ###

save(rld.df, colData, file="LoadedRNA_rld&traits_BM10_feb2019.Rdata")

### Wald test for maternal treatment ###

ddsM_NU<-DESeq(dds, minReplicatesForReplace=Inf) 
resM_NU<-results(ddsM_NU, contrast=c('treat_maternal', 'N', 'U')) #here is where the two contrasting conditions get defined
mcols(resM_NU,use.names=TRUE)
summary(resM_NU)

#out of 25332 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 928, 3.7% 
#LFC < 0 (down)   : 1516, 6% 
#outliers [1]     : 13, 0.051% 
#low counts [2]   : 4375, 17% 
#(mean count < 1)

#out of 18015 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 959, 5.3% 
#LFC < 0 (down)   : 1546, 8.6% 
#outliers [1]     : 10, 0.056% 
#low counts [2]   : 0, 0% 
#(mean count < 3)

### Wald test for developmental treatment ###

ddsD_NU<-DESeq(dds, minReplicatesForReplace=Inf) 
resD_NU<-results(ddsD_NU, contrast=c('treat_dev', 'N', 'U')) #here is where the two contrasting conditions get defined
mcols(resD_NU,use.names=TRUE)
summary(resD_NU)

#out of 25332 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 2214, 8.7% 
#LFC < 0 (down)   : 2335, 9.2% 
#outliers [1]     : 13, 0.051% 
#low counts [2]   : 7774, 31% 
#(mean count < 4)

#out of 18015 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 2297, 13% 
#LFC < 0 (down)   : 2398, 13% 
#outliers [1]     : 10, 0.056% 
#low counts [2]   : 699, 3.9% 
#(mean count < 4)

#save(rld.df,resM_NU,resD_NU, file="LoadedRNA_rld&res_gene_ALL_JUNE2019.Rdata")
save(rld.df,colData,resM_NU,resD_NU, file="LoadedRNA_rld&res_gene_BM3_JUNE2019.Rdata")

######## 

ll=load("LoadedRNA_rld&res_gene_BM3_JUNE2019.Rdata") 

head(rld.df)
head(resM_NU)

vals=(cbind(resM_NU$stat, resM_NU$pvalue, resM_NU$padj,resD_NU$stat, resD_NU$pvalue, resD_NU$padj)) #collect pvalues for each gene
colnames(vals)=c("stat_resM_NU","pval_resM_NU", "padj_resM_NU","stat_resD_NU","pval_resD_NU", "padj_resD_NU") 
rldpvals=as.data.frame(cbind(rld.df,vals)) #combine RLD with pvals
row.names(rldpvals)<-gsub("-tr","",row.names(rldpvals))
row.names(rldpvals)<-sub("transcript:","", row.names(rldpvals))

# add gene annotations #
annot=read.table("gene_info_table_header.txt",header=T,sep="\t",quote=NULL)
annot=as.data.frame(annot[-c(3:6, 11:19) ],)
genes = row.names(rldpvals)
genes2annot = match(genes,annot$spu_id)
rldpvals_annot=data.frame(GeneID=row.names(rldpvals), annot[genes2annot,], rldpvals)

# subset rld and pvals table to only those significant for each comparison #
rldpvals_M=rldpvals_annot[rldpvals_annot$padj_resM_NU<= 0.05 & !is.na(rldpvals_annot$padj_resM_NU),]
dim(rldpvals_M) #1811 18
rldpvals_D=rldpvals_annot[rldpvals_annot$padj_resD_NU<= 0.05 & !is.na(rldpvals_annot$padj_resD_NU), ]
dim(rldpvals_D) #3765 18

#save all data #
save(rldpvals,rldpvals_annot,rldpvals_M,rldpvals_D,file="rldpvals_all_bm3_JUNE2019.Rdata")

### venn diagram ###

rldpvals = as.data.frame(rldpvals)
maternal = row.names(rldpvals[rldpvals$padj_resM_NU<0.05 & !is.na(rldpvals$padj_resM_NU),])
larval = row.names(rldpvals[rldpvals$padj_resD_NU<0.05 & !is.na(rldpvals$padj_resD_NU),])
degs05 = union(maternal,larval)
length(degs05) #4931
candidates=list("Maternal"=maternal, "Developmental"=larval)

quartz()
prettyvenn=venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("coral3", "turquoise4"),
  alpha = 0.5,
  label.col = c("darkred", "white", "darkgreen"),
  cex = 1.5,
  fontfamily = "sans",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkgreen"),
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08, 0.08),
  cat.pos = 1
);
grid.draw(prettyvenn)

### principal coordinate calculation ###

dd.veg=vegdist(t(rld.df), "manhattan")
div.dd.veg=dd.veg/1000
head(div.dd.veg)

dd.pcoa=pcoa(div.dd.veg) 
head(dd.pcoa)
scores=dd.pcoa$vectors

conditions=colData
conditions$treat=as.factor(conditions$treat)
conditions$treat_maternal=as.factor(conditions$treat_maternal)
conditions$treat_dev=as.factor(conditions$treat_dev)

dd.pcoa
#$values
#   Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
#1    12.513728   0.33302308  0.274534304 0.3330231      0.2745343
#2     7.642850   0.20339627  0.183625213 0.5364194      0.4581595
#3     3.057074   0.08135674  0.138170668 0.6177761      0.5963302


### plotting PCoA ###

# plotting the first two principal axes
colorsMat <- c("coral3", "lightcoral", "mediumturquoise","turquoise4")
conditions$treat_c <- colorsMat[as.factor(conditions$treat)]

plot(scores[,1], scores[,2], col=conditions$treat_c,pch=19, xlab="PCo1", ylab="PCo2")
ordispider(scores,conditions$treat,label=F)
#ordiellipse(scores,conditions$treat,label=F)
legend("bottomright", legend = levels(conditions$treat), col=colorsMat, pch = 19)

# plotting the second and third principal axes  
plot(scores[,2], scores[,3], col=conditions$treat_c,pch=19, xlab="PCo2", ylab="PCo3")
ordispider(scores[,2:3],conditions$treat,label=F)
#ordiellipse(scores[,2:3],conditions$treat_maternal,label=F )
legend("topleft", legend = levels(conditions$treat), col=colorsMat, pch = 19)


ad=adonis(t(rld.df)~treat_maternal*treat_dev,data=conditions,method="manhattan")
#Call:
#adonis(formula = t(rld.df) ~ treat_maternal * treat_dev, data = conditions,method = "manhattan") 
#
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#                         Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#treat_maternal            1   7157886  7157886  3.4336 0.19049  0.002 ** 
#treat_dev                 1  11248146 11248146  5.3957 0.29934  0.001 ***
#treat_maternal:treat_dev  1   2492986  2492986  1.1959 0.06634  0.254    
#Residuals                 8  16677140  2084643         0.44382           
#Total                    11  37576158                  1.00000           
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


labs=c("Maternal_condition","Developmental_condition","Interaction","Residuals")
cols=c("coral2","turquoise4","gray30","grey80")
#cols = append(gg_color_hue(3), 'grey')
labs2 = paste(labs, round(ad$aov.tab$R2[1:4]*100, digits=1))
pie(ad$aov.tab$R2[1:4],labels=labs2,col=cols,main="GE")

############ OUTPUT FOR GO-MWU TESTS ############
ll=load("LoadedRNA_rld&res_gene_BM3_JUNE2019.Rdata")
#for maternal condition
resM_NU$lp=-log(resM_NU$pvalue,10)
resM_NU$lp[resM_NU$log2FoldChange<0]=-resM_NU$lp[resM_NU$log2FoldChange<0]

row.names(resM_NU)<-gsub("-tr","",row.names(resM_NU))
row.names(resM_NU)<-sub("transcript:","", row.names(resM_NU))

maternal_lpv=data.frame(cbind("gene"=row.names(resM_NU),"lp"=resM_NU$lp))
maternal_fc=data.frame(cbind("gene"=row.names(resM_NU),"fc"=resM_NU$log2FoldChange))
maternal_stat=data.frame(cbind("gene"=row.names(resM_NU),"stat"=resM_NU$stat))
write.csv(maternal_lpv,file="GO/maternal_wald_lpv_mrna.csv",quote=F, row.names=F)
write.csv(maternal_stat,file="GO/maternal_wald_stat_mrna.csv",quote=F, row.names=F)
write.csv(maternal_fc,file="GO/maternal_wald_fc_mrna.csv",quote=F, row.names=F)

##### output for fishers test for maternal condition
resM_NU=data.frame(resM_NU)
resM_NU$pvalue_go[resM_NU$pvalue<=0.01]<-1
resM_NU$pvalue_go[resM_NU$pvalue>0.01]<-0 
resM_NU$pvalue_go=as.character(resM_NU$pvalue_go) 
summary(resM_NU$pvalue_go)
#    0     1 
#15743  2262 
resM_NU=resM_NU[complete.cases(resM_NU),]
resM_NU$names=row.names(resM_NU)
maternal_fishers=data.frame(cbind("gene"=resM_NU$names,"pvalue"=resM_NU$pvalue_go))
maternal_fishers$pvalue=as.factor(maternal_fishers$pvalue)
write.csv(maternal_fishers,file="GO/maternal_wald_fishers_pval0.01.csv",quote=F, row.names=F)


#for larval condition
resD_NU$lp=-log(resD_NU$pvalue,10)
resD_NU$lp[resD_NU$log2FoldChange<0]=-resD_NU$lp[resD_NU$log2FoldChange<0]

row.names(resD_NU)<-gsub("-tr","",row.names(resD_NU))
row.names(resD_NU)<-sub("transcript:","", row.names(resD_NU))


dev_lpv=data.frame(cbind("gene"=row.names(resD_NU),"lp"=resD_NU$lp))
dev_fc=data.frame(cbind("gene"=row.names(resD_NU),"fc"=resD_NU$log2FoldChange))
dev_stat=data.frame(cbind("gene"=row.names(resD_NU),"stat"=resD_NU$stat))
write.csv(dev_lpv,file="GO/larval_wald_lpv_mrna.csv",quote=F, row.names=F)
write.csv(dev_stat,file="GO/larval_wald_stat_mrna.csv",quote=F, row.names=F)
write.csv(dev_fc,file="GO/larval_wald_fc_mrna.csv",quote=F, row.names=F)

##### output for fishers test for developmental condition
resD_NU=data.frame(resD_NU)
resD_NU$pvalue_go[resD_NU$pvalue<=0.01]<-1
resD_NU$pvalue_go[resD_NU$pvalue>0.01]<-0 
resD_NU$pvalue_go=as.character(resD_NU$pvalue_go) 
summary(resD_NU$pvalue_go)
#    0     1 
#14290  3715 
resD_NU=resD_NU[complete.cases(resD_NU),]
resD_NU$names=row.names(resD_NU)
larval_fishers=data.frame(cbind("gene"=resD_NU$names,"pvalue"=resD_NU$pvalue_go))
larval_fishers$pvalue=as.factor(larval_fishers$pvalue)
write.csv(larval_fishers,file="GO/larval_wald_fishers_pval0.01.csv",quote=F, row.names=F)





