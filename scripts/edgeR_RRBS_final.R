

############
library('methylKit')
library('edgeR')
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
library('gridExtra')

NN_1 <- read.delim("NN1.cov", header=F, nrows=6)
NN_2 <- read.delim("NN2.cov", header=F, nrows=6)
NN_3 <- read.delim("NN3.cov", header=F, nrows=6)
NU_1 <- read.delim("NU1.cov", header=F, nrows=6)
NU_2 <- read.delim("NU2.cov", header=F, nrows=6)
NU_3 <- read.delim("NU3.cov", header=F, nrows=6)
UN_1 <- read.delim("UN1.cov", header=F, nrows=6)
UN_2 <- read.delim("UN2.cov", header=F, nrows=6)
UN_3 <- read.delim("UN3.cov", header=F, nrows=6)
UU_1 <- read.delim("UU1.cov", header=F, nrows=6)
UU3_1 <- read.delim("UU3_1.cov", header=F, nrows=6)
UU3_2 <- read.delim("UU3_2.cov", header=F, nrows=6)

targets <- data.frame(row.names= c("NN1","NN2","NN3","NU1","NU2","NU3","UN1","UN2","UN3","UU1","UU3_1","UU3_2"))
targets$treat = substr(row.names(targets), 1,2)
targets$treat_maternal = substr(row.names(targets), 1,1)
targets$treat_dev = substr(row.names(targets), 2,2)

Sample <- row.names(targets)
files <- paste0(Sample, ".cov")
yall<- readBismark2DGE(files, sample.names=Sample)
dim(yall) #9354108      24

Methylation <- gl(2,1,ncol(yall), labels=c("Me","Un"))
Me <- yall$counts[, Methylation=="Me"]
Un <- yall$counts[, Methylation=="Un"]
Coverage <- Me + Un
head(Coverage)

HasCoverage <- rowSums(Coverage >= 10) == 12
table(HasCoverage)
#HasCoverage
#  FALSE    TRUE 
#9107424  246684 

y <- yall[HasCoverage ,, keep.lib.sizes=FALSE]

TotalLibSize <- y$samples$lib.size[Methylation=="Me"] + y$samples$lib.size[Methylation=="Un"]
y$samples$lib.size <- rep(TotalLibSize, each=2)
y$samples

Me <- y$counts[, Methylation=="Me"]
Un <- y$counts[, Methylation=="Un"]
M <- log2(Me + 2) - log2(Un + 2)
colnames(M) <- row.names(targets)
#plotMDS(M, col=rep(1:4, each=3), main="M-values")

#####
designSL <- model.matrix(~0+treat_maternal+treat_dev, data=targets)
design <- modelMatrixMeth(designSL)
y <- estimateDisp(y, design=design, trend="none")
fit <- glmFit(y, design)

contrM <- makeContrasts( Maternal = treat_maternalN - treat_maternalU, levels=design)
lrtM <- glmLRT(fit, contrast=contrM)
summary(decideTests(lrtM))
#       1*treat_maternalN -1*treat_maternalU
#Down                                    343
#NotSig                               246071
#Up                                      270
plotMD(lrtM)
topTags(lrtM)

out <- as.data.frame(topTags(lrtM, n=Inf, adjust.method="BH"))
keep <- out$FDR <= 0.05 
out = out[keep,]
out$id <- paste(out$Chr,out$Locus,  sep=".")

contrD <- makeContrasts( Dev = treat_devU, levels=design)
lrtD <- glmLRT(fit, contrast=contrD)
summary(decideTests(lrtD))
# 0 
save(out, lrtM, lrtD, file="edgeR_Spurp_ensembl_CpGs_lo.10.Rdata")

############# load methylKit results to compare
d3=load("methylKit_Spurp_ensembl_Maternal_CpGs_lo.10.Rdata")
myDiff25pM=getData(myDiff25p)

myDiff25pM=as.data.frame(myDiff25pM)
myDiff25pM$id <- paste(myDiff25pM$chr,myDiff25pM$start,  sep=".")

maternal_sigCpGs=merge(myDiff25pM, out, by="id")

#### make Venn for comparison between methods
edgeR=out$id
methylKit=myDiff25pM$id

degs05 = union(edgeR,methylKit)

length(degs05) #887

candidates=list("edgeR"=edgeR, "methylKit"=methylKit)

quartz()
prettyvenn=venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("coral3", "turquoise4"),
  alpha = 0.5,
  label.col = c("darkred", "black", "darkgreen"),
  cex = 1.5,
  fontfamily = "sans",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkgreen"),
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08, 0.08),
  cat.pos = 1, 
  main = "minimum depth 30"
);
grid.draw(prettyvenn)


