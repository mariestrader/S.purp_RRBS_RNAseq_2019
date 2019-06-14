

setwd("/Users/mariestrader/Dropbox/RRBS_Loaded/bismarkOutput/Ensembl")


############
library('methylKit')
library('genomation')
library('GenomicRanges')
library('reshape')
library('ggplot2')
library('vegan')
library('rgl')
library('ape')
library('adegenet')
library('VennDiagram')
library('dplyr')

file.list=list("NN1.cov",
                "NN2.cov", 
                "NN3.cov", 
                "NU1.cov", 
                "NU2.cov", 
                "NU3.cov", 
                "UN1.cov", 
                "UN2.cov", 
                "UN3.cov", 
                "UU1.cov", 
                "UU3_1.cov", 
                "UU3_2.cov")

myobj=methRead(file.list,
sample.id=list("NN1","NN2","NN3","NU1","NU2","NU3","UN1","UN2","UN3","UU1","UU3_1","UU3_2"),
           assembly="S.pur3.1",
           treatment=c(0,0,0,0,0,0,1,1,1,1,1,1),
           context="CpG",
           pipeline="bismarkCoverage"
           )
                  
#### filter sites for low coverage and high coverage ###
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=99.9)

### get methylation and coverage stats. Keep in mind the number in [[]] refers to the sample file ###
getMethylationStats(filtered.myobj[[1]],plot=FALSE,both.strands=FALSE)
#methylation statistics per base
#summary:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.00    0.00    0.00   22.69   27.27  100.00 
#percentiles:
#       0%       10%       20%       30%       40%       50%       60%       70%       80%       90%       95% 
#  0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000  11.53846  65.00000 100.00000 100.00000 
#      99%     99.5%     99.9%      100% 
#100.00000 100.00000 100.00000 100.00000 

getCoverageStats(filtered.myobj[[1]],plot=FALSE,both.strands=FALSE)
#read coverage statistics per base
#summary:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  10.00   13.00   18.00   27.51   29.00  684.00 
#percentiles:
#   0%   10%   20%   30%   40%   50%   60%   70%   80%   90%   95%   99% 99.5% 99.9%  100% 
#   10    11    12    14    16    18    21    26    33    48    70   176   251   474   684  

### unite sites across all samples ###
meth=unite(filtered.myobj, destrand=FALSE)
dim(meth) #245343     40
meth.df=as.data.frame(meth)
write.table(meth.df,file="meth_CpGs_loaded.tab",sep="\t", quote=F, row.names=F, col.names=T)

#### differential methylation analysis for maternal condition ###
myDiff=calculateDiffMeth(meth)
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
#328

myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
#356

myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
#684

### get percent methylation per site per sample ###
perc.meth=percMethylation(meth)
perc.meth=as.data.frame(perc.meth)

### save all files for meth analysis for maternal condition ###
save(filtered.myobj, meth, meth.df, perc.meth, myDiff, myDiff25p, myDiff25p.hyper, myDiff25p.hypo, file="methylKit_Spurp_ensembl_Maternal_CpGs.Rdata")

############################################# test by larval condition
file.list=list("NN1.cov",
                "NN2.cov", 
                "NN3.cov", 
                "NU1.cov", 
                "NU2.cov", 
                "NU3.cov", 
                "UN1.cov", 
                "UN2.cov", 
                "UN3.cov", 
                "UU1.cov", 
                "UU3_1.cov", 
                "UU3_2.cov")

myobj=methRead(file.list,
sample.id=list("NN1","NN2","NN3","NU1","NU2","NU3","UN1","UN2","UN3","UU1","UU3_1","UU3_2"),
           assembly="GCF_000002235.4_Spur_4.2_genomic.fa",
           treatment=c(0,0,0,1,1,1,0,0,0,1,1,1),
           context="CpG",
           pipeline="bismarkCoverage"
           )
           
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=99.9)

meth=unite(filtered.myobj, destrand=FALSE)
myDiff=calculateDiffMeth(meth)
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
#175

myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
#41 

myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
#216

### get percent methylation per site per sample ###
perc.meth=percMethylation(meth)
perc.meth=as.data.frame(perc.meth)

### save all files for meth analysis for maternal condition ###
save(myobj, filtered.myobj, meth, meth.df, perc.meth, myDiff, myDiff25p, myDiff25p.hyper, myDiff25p.hypo, file="methylKit_Spurp_ensembl_Developmental_CpGs.Rdata")

################################### Get methylation annotation and output for comparison

# read annotation file for genes #
meth.annot=read.table("meth_CpGs_loaded_annotated.tab")
length(unique(meth.annot$V41)) #9219
meth.annot$id <- paste(meth.annot$V1,meth.annot$V2,meth.annot$V3, sep=".")

perc.meth=percMethylation(meth, rowids=T)
perc.meth=as.data.frame(perc.meth)
perc.meth$id=row.names(perc.meth)

meth.annot.perc.meth=merge(perc.meth, meth.annot, by="id")
meth.annot.perc.meth=as.data.frame(meth.annot.perc.meth[-c(14:53, 55) ])
dim(meth.annot.perc.meth) #95197    15

meth.annot.perc.meth$median = apply(meth.annot.perc.meth[,2:13], 1, median)
summary(meth.annot.perc.meth)
hist(meth.annot.perc.meth$median) #Figure S4A

meth<-meth.annot.perc.meth$median>"1" 
meth.annot.perc.meth$meth=meth

summary(meth.annot.perc.meth$meth)
# values for median>1
#   Mode   FALSE    TRUE 
#logical   53257   41940 

dat_med <- 
  # split data by gene: 
  group_by(meth.annot.perc.meth, V41) %>%
  # make summary table. to get total cpg sites for each gene just take the length of the column: 
  summarize(total_cpg = length(id),
            # If methylated column is true/false, trues count as 1 and falses count as 0
            # so sum(methylated) should give you number of "true"s
            total_cpg_methylated = sum(meth),
            percent_methylated =total_cpg_methylated/total_cpg)

dat_med=as.data.frame(dat_med)

#filter for genes that have at least 5 CpGs represented
dat5_med=dat_med[dat_med$total_cpg>5,]
dim(dat5_med) #5034
hist(dat5_med$percent_methylated, breaks=20) #Figure S4B

save(meth.annot.perc.meth, dat_med, dat5_med, file="GS-CpGs.Rdata")

########################################################### multivariate analysis

ll=load("methylKit_Spurp_ensembl_Maternal_CpGs.Rdata")

dd.veg=vegdist(t(perc.meth), "manhattan")
div.dd.veg=dd.veg/100000
head(div.dd.veg)

dd.pcoa=pcoa(div.dd.veg) 
head(dd.pcoa)
scores=dd.pcoa$vectors
#$values
#   Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
#1    110.25115   0.13198560  0.274534304 0.1319856      0.2745343
#2     88.03524   0.10539014  0.183625213 0.2373757      0.4581595
#3     83.92882   0.10047420  0.138170668 0.3378499      0.5963302

colData <- data.frame(row.names= colnames(perc.meth), bucket= c("NN1","NN2","NN3","NU1","NU2","NU3","UN1","UN2","UN3","UU1","UU2","UU3"))
colData$treat = substr(colData$bucket, 1,2)
colData$treat_maternal = substr(colData$bucket, 1,1)
colData$treat_dev = substr(colData$bucket, 2,2)

conditions=colData
conditions$treat=as.factor(conditions$treat)
conditions$treat_maternal=as.factor(conditions$treat_maternal)
conditions$treat_dev=as.factor(conditions$treat_dev)

# plotting the significant maternal CDS's
# plotting the first two principal axes
colorsMat <- c("coral3", "lightcoral", "mediumturquoise","turquoise4")
conditions$treat_c <- colorsMat[as.factor(conditions$treat)]


plot(scores[,1], scores[,2], col=conditions$treat_c,pch=19, xlab="PCo1", ylab="PCo2")
ordispider(scores,conditions$treat,label=F)
#ordiellipse(scores,conditions$treat,label=F)
legend("bottomright", legend = levels(conditions$treat), col=colorsMat, pch = 19)

ad=adonis(t(perc.meth)~treat_maternal*treat_dev,data=conditions,method="manhattan") 
Call:
#adonis(formula = t(perc.meth) ~ treat_maternal * treat_dev, data = conditions,      method = "manhattan") 

#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#                         Df  SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)    
#treat_maternal            1 1.0732e+12 1.0732e+12 1.49028 0.12847  0.001 ***
#treat_dev                 1 8.1131e+11 8.1131e+11 1.12666 0.09713  0.088 .  
#treat_maternal:treat_dev  1 7.0795e+11 7.0795e+11 0.98313 0.08475  0.586    
#Residuals                 8 5.7608e+12 7.2011e+11         0.68965           
#Total                    11 8.3533e+12                    1.00000  
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

labs=c("Maternal_condition","Developmental_condition","Interaction","Residuals")
cols=c("coral2","turquoise4","gray30","grey80")
labs2 = paste(labs, round(ad$aov.tab$R2[1:4]*100, digits=1))
pie(ad$aov.tab$R2[1:4],labels=labs2,col=cols,main="CpG Meth")

####### Venn
ll=load("methylKit_Spurp_ensembl_Maternal_CpGs.Rdata")
myDiff25p$uniq <- paste(myDiff25p$chr,myDiff25p$start)
maternal=myDiff25p$uniq

ll=load("methylKit_Spurp_ensembl_Developmental_CpGs.Rdata")
myDiff25p$uniq <- paste(myDiff25p$chr,myDiff25p$start)
larval=myDiff25p$uniq

degs05 = union(maternal,larval)

length(degs05) #887

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


