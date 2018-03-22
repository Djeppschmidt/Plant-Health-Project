

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
