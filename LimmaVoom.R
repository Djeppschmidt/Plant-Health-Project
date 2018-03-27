

data<-readRDS("QPCR_Corrected_totalPS.rds")

otu<-otu_table(data)
#update...
#from: https://www.bioconductor.org/packages/3.7/bioc/vignettes/limma/inst/doc/usersguide.pdf
# pg 69.
dge <- DGEList(counts=counts)

keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

logCPM <- cpm(dge, log=TRUE, prior.count=3)
fit <- lmFit(logCPM, design)
#fit <- eBayes(fit, trend=TRUE) #pick one
#fit <- treat(fit, lfc=log2(1.2))
topTable(fit, coef=ncol(design))
#topTreat(fit, coef=ncol(design))# for second option (treat)

##################################################

library(limma)
#import phyloseq

#setwd("~/Desktop/PhD/PHP/")
bacarc<-readRDS("phyloseq_Bac_metadata_seqtab.rds")
bacarc2<-prune_samples(sample_data(bacarc)$Glyphosphate_Treatment!="undetermined", bacarc)
otutab<-otu_table(bacarc2)
sData<-sample_data(bacarc2)

library(edgeR)
?DGEList
otutab<-as.matrix(t(otu_table(bacarc2)))

dge<-DGEList(counts=otutab)
#do a scale normalization
dge<-calcNormFactors(dge)

head(Glyphosphate_Treatment)
Glyphosphate_Treatment<-sData$Glyphosphate_Treatment
Glyphosphate_Treatment<-lapply(Glyphosphate_Treatment, function(x) {
  gsub("no_spray", 0, x)})
Glyphosphate_Treatment<-lapply(Glyphosphate_Treatment, function(x) {
  gsub("spray", 1, x)})
Glyphosphate_Treatment<-as.numeric(Glyphosphate_Treatment)

design<-model.matrix(~Glyphosphate_Treatment)
v <- voom(dge, design, plot=TRUE)

fit<-lmFit(v, design)
