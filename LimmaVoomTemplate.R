# limma-Voom ####
library(limma)
library(Glimma)
library(edgeR)
sv.RRCorn.otu<-as.matrix(t(otu_table(sv.RRcorn)))#prepare OTU table for analysis
#head(otutab)#inspect the OTU table
#?DGEList#figure out how to set up next step...
sv.RRCorn.dge<-DGEList(counts=sv.RRCorn.otu)#make OTU table into DGEList object
#do a scale normalization
sv.RRCorn.dge<-calcNormFactors(sv.RRCorn.dge)

#make df for model matrix!
sv.RRCorn.sam.df<-sample_data(sv.RRcorn)
attach(sv.RRCorn.sam.df)
head(sv.RRCorn.sam.df)
#
design<-model.matrix(~year*System.loc*Glyphosphate_Treatment*Soil_Zone*Sampling_date+Loc_plot_ID) #specify the model
sv.RRCorn.v <- voom(sv.RRCorn.dge, design, plot=TRUE) #make the voom object to test
sv.RRCorn.fit<-lmFit(sv.RRCorn.v, design)
sv.RRCorn.contr<-makeContrasts(categories , levels=colnames(coef(sv.RRCorn.fit))) #make treatment contrasts

?contrasts.fit
sv.RRCorn.tmp<-contrasts.fit(sv.RRCorn.fit, sv.RRCorn.contr) #run contrasts
?ebayes
sv.RRCorn.tmp<-eBayes(sv.RRCorn.tmp) #use bayesian statistics to determine significantly different taxa! ranks genes in order of evidence of differential expression!

sv.RRCorn.tmp2<-topTable(sv.RRCorn.tmp, coef=1, sort.by="P", n="INF") #make output table
sv.RRCorn.tmp2
