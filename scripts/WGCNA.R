setwd("/Users/mariestrader/Dropbox/RRBS_Loaded/RNAseq/WGCNA")
library(WGCNA)
library(flashClust)
options(stringsAsFactors=FALSE)

lnames= load(file="LoadedRNA_rld&traits_BM10_feb2019.Rdata")
morph=read.csv("morph_WGCNA.csv")
oxy=read.csv("oxygenConsumpLoaded.csv")

colnames(morph)=c("sample","height","spiculeLength")
colnames(oxy)=c("sample","oxygenConsump")

allTraits=merge(morph, oxy, by='sample', all.x=T)

allTraits$treat = substr(allTraits$sample, 1,2)
allTraits$treat_maternal = substr(allTraits$sample, 1,1)
allTraits$treat_dev = substr(allTraits$sample, 2,2)

NN<-allTraits$treat=="NN"
NN[NN==T]<- 1
NN[NN==F]<- 0
allTraits$NN=NN

NU<-allTraits$treat=="NU"
NU[NU==T]<- 1
NU[NU==F]<- 0
allTraits$NU=NU

UN<-allTraits$treat=="UN"
UN[UN==T]<- 1
UN[UN==F]<- 0
allTraits$UN=UN

UU<-allTraits$treat=="UU"
UU[UU==T]<- 1
UU[UU==F]<- 0
allTraits$UU=UU

upwellingM<-allTraits$treat_maternal=="U"
upwellingM[upwellingM==T]<- 1
upwellingM[upwellingM==F]<- 0
allTraits$upwellingM=upwellingM

nonupwellingM<-allTraits$treat_maternal=="N"
nonupwellingM[nonupwellingM==T]<- 1
nonupwellingM[nonupwellingM==F]<- 0
allTraits$nonupwellingM=nonupwellingM

upwellingD<-allTraits$treat_dev=="U"
upwellingD[upwellingD==T]<- 1
upwellingD[upwellingD==F]<- 0
allTraits$upwellingD=upwellingD

nonupwellingD<-allTraits$treat_dev=="N"
nonupwellingD[nonupwellingD==T]<- 1
nonupwellingD[nonupwellingD==F]<- 0
allTraits$nonupwellingD=nonupwellingD

allTraits=as.data.frame(allTraits[-c(5:7)])
datTraits=allTraits
######################################
dim(rld.df) #15142    12
datExpr= as.data.frame(t(rld.df[,]))
dim(datExpr) #12 15142

save(datExpr, datTraits, file="LoadedRNA_bm10_SamplesAndTraits.RData")

################################# LOAD INTO CLUSTER, REST OF THIS PERFORMED THERE
library(WGCNA)
library(flashClust)
ll=load("LoadedRNA_bm10_SamplesAndTraits.RData")

gsg=goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK #If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:

#if (!gsg$allOK)
#	{if (sum(!gsg$goodGenes)>0)
#		printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse= ", ")));
#		if (sum(!gsg$goodSamples)>0)
#			printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse=", ")))
#		datExpr0= datExpr0[gsg$goodSamples, gsg$goodGenes]
#		}

#Now we cluster the samples to find outliers
sampleTree= flashClust(dist(datExpr), method="average")
sizeGrWindow(12,9)
par(cex=0.6)
par(mar= c(0,4,2,0))
plot(sampleTree, main= "Sample Clustering to Detect Outliers", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)
abline(h=30, col="red")

#form a data frame analogous to expression data that will hold the traits.
rownames(datTraits) <- datTraits$sample
datTraits$sample <- NULL
datTraits= as.data.frame(datTraits)
table(rownames(datTraits)==rownames(datExpr)) #should return TRUE if datasets align correctly, otherwise your names are out of order

#expression data is in datExpr, corresponding clinical traits are datTraits
sampleTree2=flashClust(dist(datExpr), method="average")
traitColors= numbers2colors(datTraits, signed= FALSE)
plotDendroAndColors(sampleTree2, traitColors, groupLabels= names(datTraits), main="Sample Dendrogram and Trait heatmap")
save(datExpr, datTraits, file="Loaded_bm10_SamplesAndTraits_ready.RData")

########################## Run this as a script "SoftThresh_TACC.R", ran for <5 mins
library(WGCNA)
library(flashClust)
ll=load("Loaded_bm10_SamplesAndTraits_ready.RData")
options(stringsAsFactors = FALSE);
###Step by step network construction and module detection

powers= c(c(1:10), seq(from =12, to=20, by=2)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function

sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
############################ 

#from this plot, we would choose a power of 9 becuase it's the lowest power for which the scale free topology index reaches 0.90

########################### Run this as script "TOM.R" or "TOM_merged.R"
library(WGCNA)
library(flashClust)
ll=load("Loaded_bm10_SamplesAndTraits_ready.RData")
options(stringsAsFactors = FALSE);

softPower = 9;
adjacency = adjacency(datExpr, power = softPower, type='signed');

TOM = TOMsimilarity(adjacency, TOMType='signed');
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
#pdf(file=Dendrogram_signed_BM10.pdf, width=20, height=20)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
#dev.off()

minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)
#dynamicMods
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
#2676 1863 1423 1395 1195  964  907  713  702  597  542  457  441  277  250  246 
#  17   18   19   20 
# 185  140   85   84 

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
#       black         blue        brown         cyan        green  greenyellow 
#         907         1863         1423          277         1195          542 
#      grey60    lightcyan   lightgreen  lightyellow      magenta midnightblue 
#         185          246          140           85          702          250 
#        pink       purple          red    royalblue       salmon          tan 
#         713          597          964           84          441          457 
#   turquoise       yellow 
#        2676         1395 
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
#pdf(file=Dendrogram_signed_BM10_colors.pdf, width=20, height=20)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
#sizeGrWindow(7, 6)
#pdf(file=ClusteringEigengenes.pdf, width=20, height=20)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
#dev.off()

MEDissThres = 0.2 
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

#sizeGrWindow(12, 9)
#pdf(file = "DendroAndColors_sft6_bm10.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "LoadedRNA_WGCNA_networkConstruct_signedsft6_bm10_merged.RData")

############################################# finished in less than 5 mins
library(WGCNA)
library(flashClust)
options(stringsAsFactors = FALSE);
ll=load("Loaded_bm10_SamplesAndTraits_ready.RData")
ll2= load(file="LoadedRNA_WGCNA_networkConstruct_signedsft6_bm10_merged.RData")

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors,softPower=9)$eigengenes #change softPower for host=6 or symbiont=10
MEs = orderMEs(MEs0)

#correlations of traits with eigengenes
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

#############################################           
#plotting massive table of all information - module membership, genes, gene names, etc.
#setwd("/Users/mariestrader/Dropbox/MEDIP/R_AnalysisMarch2018")
annot=read.table("gene_info_table_header.txt",header=T,sep="	\t",quote=NULL)
names(datExpr)<-gsub("-tr","",names(datExpr))
names(datExpr)<-sub("transcript:","", names(datExpr))

probes = names(datExpr)
probes2annot = match(probes,annot$spu_id)
datGS.Traits=data.frame(cor(datExpr,datTraits,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datME=moduleEigengenes(datExpr,moduleColors)$eigengenes
datKME=signedKME(datExpr, datME, outputColumnName="MM.")
datOutput=data.frame(ProbeID=names(datExpr), annot[probes2annot,],moduleColors,datKME,datGS.Traits)
dim(datOutput) #15142    44

meth.data=load("sigCpGdata.Rdata")
length(unique(sig_CpG$gene)) #143

methTF <- sapply(datOutput$spu_id, function(x) { # for each value in df1$x
  j <- which(sig_CpG$gene == x) # check if df2$y has a match
  ifelse(length(j) > 0, T, F) # if there is, give the location in the vector
}) # if not give NA
methTF.df=as.data.frame(methTF, row.names=names(methTF))
methTF.df$spu_id=row.names(methTF.df)
datOutput.m<- merge(datOutput, methTF.df, by="spu_id")

datOutput.meth=subset(datOutput.m, datOutput.m$methTF==T)

write.table(datOutput.meth,"Loaded_annotatedNetworkAnalysisResult_bm10_merged0.2_signed_methGenes.txt",row.names=F,sep="\t", col.names=T, quote=F)
write.table(datOutput.m,"Loaded_annotatedNetworkAnalysisResult_bm10_merged0.2_signed_meth.txt",row.names=F,sep="\t", col.names=T, quote=F)

###Making tables for GO analysis, categorical, interesting modules vs entire expression set

table(moduleColors)
#moduleColors
#      black        blue        cyan      grey60  lightgreen lightyellow 
#       2330        2405        1229         185         140          85 
#       pink         red   royalblue      salmon   turquoise      yellow 
#       1310        1421          84        1882        2676        1395 

table(datOutput.meth$moduleColors) #number of genes with sigCpG in each module
#     black       blue       cyan     grey60 lightgreen       pink        red 
#        19         26          9          2          2          7         11 
# royalblue     salmon  turquoise     yellow 
#         1         16         20         11 

#Percent of genes in module with sig CpG
#black 19/2330=0.008154506 - Dev. module
#blue 26/2405=0.01081081 - Maternal module
#cyan 9/1229=0.007323027 - Maternal module
#grey60 2/185=0.01081081
#lightgreen 2/140=0.01428571 - Maternal module
#pink 7/1310=0.005343511 - Dev module
#red 11/1421=0.007741027 - Maternal module, Dev. module
#royalblue 1/84=0.01190476 -Maternal module
#salmon 16/1882=0.008501594
#turquoise 20/2676=0.007473842 - Dev. module
#yellow 11/1395=0.007885305 - Dev. module

Perc.meth.mod=data.frame("Treatment"=c("dev","dev","dev","dev","dev","mat","mat","mat","mat","mat"), "PercGenesWithDMCpGs"=c("0.008154506","0.005343511","0.007741027","0.007473842","0.007885305","0.01081081","0.007323027","0.01428571","0.007741027","0.01190476"))
Perc.meth.mod$PercGenesWithDMCpGs=as.numeric(Perc.meth.mod$PercGenesWithDMCpGs)

boxplot(PercGenesWithDMCpGs~Treatment, data=Perc.meth.mod, ylab="PercGenesWithDMCpGs")
t.test(PercGenesWithDMCpGs~Treatment, data=Perc.meth.mod) #not significant
full=MCMCglmm(PercGenesWithDMCpGs~Treatment,data=Perc.meth.mod,nitt=50000)
summary(full) #not significant 

##############Go categorical
datOutput <- read.table("Loaded_annotatedNetworkAnalysisResult_bm10_merged0.2_signed.txt", header=T, sep="\t")

col="yellow"

#Categorical
tab=datOutput[,c(1,21)]

#tab2=rbind(tab,rest)
tab$moduleColors=as.character(tab$moduleColors)

tab$moduleColors[tab$moduleColors!=col]<-0
tab$moduleColors[tab$moduleColors==col]<-1 
tab$moduleColors=as.factor(tab$moduleColors) 
summary(tab) #do counts match table of module colors?
#tab$ProbeID<- sub("^", "gene:", tab$ProbeID )

write.csv(tab,file="GO_yellow_categorical.csv",quote=F,row.names=F)

####### output for genes with Sig DF CpGs
tab=datOutput[,c(1,45)]

tab$methTF[tab$methTF==T]<-1
tab$methTF[tab$methTF==F]<-0
tab$methTF=as.factor(tab$methTF) 
summary(tab) #do counts match table of module colors?

write.csv(tab,file="GO_meth_fishers.csv",quote=F,row.names=F)

#################
#Gene relationship to trait and important modules: Gene Significance and Module membership
#datOutput <- read.table("Loaded_annotatedNetworkAnalysisResult_bm10_merged0.2_signed_methGenes.txt", header=T, sep="\t")
length(unique(datOutput.m$spu_id)) #15112

library(scales)
CUT=0.01
P.TYPE='pvalue'
alpha = 1
XLIM=c(0,1)
YLIM=c(0,1)
text.x = 2.75
text.y=-3
CEX=.75
greys = grey.colors(2)
par(mfrow=c(1,1))

# Ran this for modules associated with spicule length
scatter=subset(datOutput.m, moduleColors=="blue")
scatter.sig=subset(datOutput.m, moduleColors=="blue" & methTF==T) #genes that has a Sig CpG in them

plot(abs(scatter$cor.upwellingM) ~ abs(scatter$MM.blue), col = alpha('blue', alpha), xlim = XLIM, ylim = YLIM, xlab = expression(paste("MM.blue")), ylab = expression(paste("cor.upwellingM")), axes = F, mgp=c(2.1,1,0), cex=CEX)
axis(1); axis(2, las=2);box()
points(abs(scatter.sig$cor.upwellingM) ~ abs(scatter.sig$MM.blue), col = 'black', pch = 19, xlim = XLIM, ylim = YLIM, cex=CEX)


#lm1 = lm(combo$log2FoldChange ~ combo$class)
#abline(lm1, col = 'blue')
#summary(lm1)
#lm2 = lm(combo_sig$log2FoldChange ~ combo_sig$class)
#abline(lm2, col = 'firebrick')
#summary(lm2)
#plot(density(na.omit(combo$class)), main='', xlab='', xlim=XLIM, axes=F, col='black');axis(1)
#lines(density(na.omit(combo_sig$class)), col='firebrick')


#Gene-trait significance correlation plots

par(mfrow=c(2,3))
module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
 abs(geneTraitSignificance[moduleGenes, 1]),
                 xlab = paste("ModMem in", module, "module"),
                 ylab = "Gene Sig for spicule Length",
                 main = paste("MM vs. GS\n"),
                 cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module,
                 plotPriority(tab$methTF=T, col= "black"))
  
  

###############################Pull out genes from interesting modules
#moduleColors
#      black        blue        cyan      grey60  lightgreen lightyellow 
#       2330        2405        1229         185         140          85 
#       pink         red   royalblue      salmon   turquoise      yellow 
#       1310        1421          84        1882        2676        1395

lnames= load(file="LoadedRNA_rld&traits_BM10_feb2019.Rdata")
row.names(rld.df)<-gsub("-tr","",row.names(rld.df))
row.names(rld.df)<-sub("transcript:","", row.names(rld.df))

names(datExpr)<-gsub("-tr","",names(datExpr))
names(datExpr)<-sub("transcript:","", names(datExpr))

head(rld.df)
vsd.wgcna=as.data.frame(rld.df)
row.names(vsd.wgcna) -> vsd.wgcna$X

cands=names(datExpr[moduleColors=="yellow"])
c.vsd=vsd.wgcna[vsd.wgcna$X %in% cands,]
head(c.vsd)
row.names(c.vsd)=c.vsd$X
c.vsd=c.vsd[,-13]
head(c.vsd)
length(c.vsd[,1])

#setwd("/Users/mariestrader/Documents/CompDevPaper/Analysis/WGCNA/RLD_BM10_sft14_merged0.1_signed")
write.csv(c.vsd,"rld_yellow.csv",quote=F)


####
##########To output ME by sample and plot boxplots
datOutput <- read.table("Loaded_annotatedNetworkAnalysisResult_bm10_merged0.2_signed.txt", header=T, sep="\t")
ll <-load("LoadedRNA_WGCNA_networkConstruct_signedsft6_bm10_merged.RData")

meout<-data.frame(cbind(row.names(datExpr),MEs))

meout$treat = substr(meout$row.names.datExpr., 1,2)
meout$treat_maternal = substr(meout$row.names.datExpr., 1,1)
meout$treat_dev = substr(meout$row.names.datExpr., 2,2)

meout$treat=as.factor(meout$treat)

boxplot(MEblack~treat, data=meout, ylab="Black Eigengene Expression")
boxplot(MEyellow~treat, data=meout, ylab="Yellow Eigengene Expression")
boxplot(MElightyellow~treat, data=meout, ylab="Lightyellow Eigengene Expression")
boxplot(MEcyan~treat, data=meout, ylab="Cyan Eigengene Expression")
boxplot(MEroyalblue~treat, data=meout, ylab="RoyalBlue Eigengene Expression")
boxplot(MElightgreen~treat, data=meout, ylab="Lightgreen Eigengene Expression")
boxplot(MEsalmon~treat, data=meout, ylab="Salmon Eigengene Expression")
boxplot(MEpink~treat, data=meout, ylab="pink Eigengene Expression")
boxplot(MEturquoise~treat, data=meout, ylab="Turquoise Eigengene Expression")
boxplot(MEblue~treat, data=meout, ylab="Blue Eigengene Expression")
boxplot(MEred~treat, data=meout, ylab="Red Eigengene Expression")

write.csv(meout,"MEbySample_forboxplots.csv",quote=F,row.names=F)

##############################heatmap of module expression with bar plot of trait of interest by sample...
ll <-load("LoadedRNA_WGCNA_networkConstruct_signedsft6_bm10_merged.RData")
ll2=load("Loaded_bm10_SamplesAndTraits_ready.RData")

sizeGrWindow(8,7);
which.module="black" #pick module of interest #black    blue   brown   green    magenta    pink     red
which.trait="spiculeLength" #change trait of interest here
datTraits2=datTraits[order((datTraits$spiculeLength),decreasing=T),]#change trait of interest here

trait=datTraits2[, paste(which.trait)]
genes=datExpr[,moduleColors==which.module ] #replace where says subgene below to plot all rather than just subset
genes=genes[rownames(datTraits2),]

#quartz()
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(genes) ),nrgcols=30,rlabels=F, clabels=rownames(genes), rcols=which.module,)
par(mar=c(5, 4.2, 0, 0.7))
barplot(trait, col=which.module, main="", cex.main=2,
        ylab="SpiculeLength",xlab="sample")#change trait of interest here

##############################heatmap of genes with sig CpGs in them

ll2=load("LoadedRNA_rld&res_gene_feb2019.Rdata")
datOutput <- read.table("Loaded_annotatedNetworkAnalysisResult_bm10_merged0.2_signed_methGenes.txt", header=T, sep="\t")
row.names(meth.rld.df)<-gsub("-tr","",row.names(meth.rld.df))
row.names(meth.rld.df)<-sub("transcript:","", row.names(meth.rld.df))
meth.rld.df=as.data.frame(meth.rld.df)
meth.rld.df$spu_id<-row.names(meth.rld.df)

combo=merge(meth.rld.df, datOutput, by="spu_id")
dim(combo) #124  57
row.names(combo) <- make.names(paste(combo$common_name,combo$moduleColors, sep="."), unique=T)

combo=combo[combo$moduleColors== "red", ]

combo=combo[2:13]
means=apply(combo,1,mean) # means of rows
explc=combo-means

library(pheatmap)
heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=0.9)(100)
#greenred=colorRampPalette(c("green","black","red"),bias=1)(100) #lower bias number gives more greens; higher bias gives more reds
#par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
#setwd("/Users/mariestrader/Documents/CompDevPaper/Analysis/WGCNA/GO/Plots")
pdf("SigDifMeth_heatmap_red.pdf",height=2,width=4)
pheatmap(explc,color=heat.colors,cluster_cols=F,border_color=NA,clustering_distance_rows="correlation")
dev.off()

##### plot only those that are also DEGs
combo=merge(meth.rld.df, datOutput, by="spu_id")
dim(combo) #124  57

resM_NU=as.data.frame(resM_NU)
row.names(resM_NU)<-gsub("-tr","",row.names(resM_NU))
row.names(resM_NU)<-sub("transcript:","", row.names(resM_NU))
resM_NU$spu_id<-row.names(resM_NU)

comboM=merge(combo, resM_NU, by="spu_id")
dim(comboM) #124  57

comboM=comboM[comboM$padj<= 0.05, ]
dim(comboM) #14  7

row.names(comboM) <- make.names(paste(comboM$spu_id,comboM$moduleColors, sep="."), unique=T)
#row.names(comboM)<-comboM$spu_id
comboM=comboM[2:13]
means=apply(comboM,1,mean) # means of rows
explc=comboM-means

library(pheatmap)
heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=0.8)(100)
#greenred=colorRampPalette(c("green","black","red"),bias=1)(100) #lower bias number gives more greens; higher bias gives more reds
#par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
#setwd("/Users/mariestrader/Documents/CompDevPaper/Analysis/WGCNA/GO/Plots")
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

#res_BL_PR=res_BL_PR[complete.cases(res_BL_PR), ]
comboD=comboD[comboD$padj<= 0.05, ]
dim(comboD) #30  7

row.names(comboD) <- make.names(paste(comboD$spu_id,comboD$moduleColors, sep="."), unique=T)
#row.names(comboD)<-comboD$spu_id
comboD=comboD[2:13]
means=apply(comboD,1,mean) # means of rows
explc=comboD-means

library(pheatmap)
heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=0.8)(100)
#greenred=colorRampPalette(c("green","black","red"),bias=1)(100) #lower bias number gives more greens; higher bias gives more reds
#par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
#setwd("/Users/mariestrader/Documents/CompDevPaper/Analysis/WGCNA/GO/Plots")
pdf("SigDifMeth_heatmap_devDEGs_spu.pdf",height=5,width=4.5)
pheatmap(explc,color=heat.colors,cluster_cols=F,border_color=NA,clustering_distance_rows="correlation")
dev.off()

meth.data=load("sigCpGdata.Rdata")
length(unique(sig_CpG$gene)) #143

library(dplyr)
sigMethSum <- 
  # split data by gene: 
  group_by(sig_CpG, gene) %>%
  # make summary table. to get total cpg sites for each gene just take the length of the column: 
  summarize(total_sig_cpg = length(id))
names(sigMethSum)[1]<-"spu_id"

megaData=merge(combo, sigMethSum, by="spu_id")
write.table(megaData,"sigCpG_ALLinfor.txt",quote=F,row.names=F, sep="\t")





