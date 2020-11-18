### Beltsville ###
library(phyloseq)
library(ggplot2)

bacarc<-readRDS("~/Desktop/PhD/PHP/bacarc_RAWcomDat.rds")
global.S<-subset_samples(bacarc, crop=="soy")
global.C<-subset_samples(bacarc, crop=="corn")


#Beltsville subset without rarefaction
Beltsville<-subset_samples(bacarc, Location=="Beltsville")
table(sample_data(Beltsville)$crop) #check to make sure formating is good...
Belt.C<-subset_samples(Beltsville, crop=="corn")
Belt.S<-subset_samples(Beltsville, crop=="soy")

#Stoneville subset without rarefaction
Stoneville<-subset_samples(bacarc, Location=="Stoneville")
Stone.S<-subset_samples(Stoneville, crop=="soy")
Stone.C<-subset_samples(Stoneville, crop=="corn")

head(sample_data(Belt.C))

# remove non-round-up ready
Belt.C<-subset_samples(Belt.C, genotype=="RR")
Belt.S<-subset_samples(Belt.S, genotype=="RR")

# remove non-round-up ready
Stone.C<-subset_samples(Stone.C, genotype=="RR")
Stone.S<-subset_samples(Stone.S, genotype=="RR")

###
###
### Start PCA!
###
###

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

library(DESeq2)
### Beltsville Corn PCA
head(sample_data(Belt.C))
ddsBelt.C = phyloseq_to_deseq2(Belt.C, ~Glyphosphate_Treatment) 
ddsBelt.C$group <- factor(paste(sample_data(Belt.C)$System.loc, sample_data(Belt.C)$Glyphosphate_Treatment))
ddsBelt.C$Sampling_date<-factor(ddsBelt.C$Sampling_date)
ddsBelt.C$Soil_Zone<-factor(ddsBelt.C$Soil_Zone)

levels(ddsBelt.C$group)<-sub(" ", "", levels(ddsBelt.C$group))
design(ddsBelt.C) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

#colData(ddsBelt.C)
#?design
#dim(ddsBelt.C)

#update model...
geoMeans = apply(counts(ddsBelt.C), 1, gm_mean)

ddsBelt.C = estimateSizeFactors(ddsBelt.C, geoMeans = geoMeans)
ddsBelt.C=estimateDispersions(ddsBelt.C)

#Belt.Cvst=getVarianceStabilizedData(ddsBelt.C, method="per-condition", sharingMode="maximum", fitType="local")
Belt.Cvst=getVarianceStabilizedData(ddsBelt.C) #not working for some reason... Error: is(cds, "CountDataSet") is not TRUE

ddsBelt.C<-DESeq(ddsBelt.C, test="Wald", fitType="parametric", parallel=T)

resBelt.CDS<-results(ddsBelt.C)

Belt.Cvst<-vst(ddsBelt.C, blind=T)
plotPCA(Belt.Cvst, intgroup="group")
plotPCA(Belt.Cvst, intgroup="System.loc")
plotPCA(Belt.Cvst, intgroup="Glyphosphate_Treatment")
plotPCA(Belt.Cvst, intgroup="Sampling_date")
plotPCA(Belt.Cvst, intgroup="Soil_Zone")
plotPCA(Belt.Cvst, intgroup="Glyphosphate_Treatment")


###
### Beltsville Soy PCA
###


ddsBelt.S = phyloseq_to_deseq2(Belt.S, ~Glyphosphate_Treatment) 
ddsBelt.S$group <- factor(paste(sample_data(Belt.S)$System.loc, sample_data(Belt.S)$Glyphosphate_Treatment))
ddsBelt.S$Sampling_date<-factor(ddsBelt.S$Sampling_date)
ddsBelt.S$Soil_Zone<-factor(ddsBelt.S$Soil_Zone)

levels(ddsBelt.S$group)<-sub(" ", "", levels(ddsBelt.S$group))
design(ddsBelt.S) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

#colData(ddsBelt.C)
#?design
#dim(ddsBelt.C)

#update model...
geoMeans = apply(counts(ddsBelt.S), 1, gm_mean)

ddsBelt.S = estimateSizeFactors(ddsBelt.S, geoMeans = geoMeans)
ddsBelt.S=estimateDispersions(ddsBelt.S)

#Belt.Cvst=getVarianceStabilizedData(ddsBelt.C, method="per-condition", sharingMode="maximum", fitType="local")
Belt.Svst=getVarianceStabilizedData(ddsBelt.S) #not working for some reason... Error: is(cds, "CountDataSet") is not TRUE

ddsBelt.S<-DESeq(ddsBelt.S, test="Wald", fitType="parametric", parallel=T)

resBelt.SDS<-results(ddsBelt.S)

Belt.Svst<-vst(ddsBelt.S, blind=T)
plotPCA(Belt.Svst, intgroup="group")
plotPCA(Belt.Svst, intgroup="System.loc")
plotPCA(Belt.Svst, intgroup="Glyphosphate_Treatment")
plotPCA(Belt.Svst, intgroup="Sampling_date")
plotPCA(Belt.Svst, intgroup="Soil_Zone")
plotPCA(Belt.Svst, intgroup="Glyphosphate_Treatment")

###
### Stoneville Soy PCA
###

ddsStone.S = phyloseq_to_deseq2(Stone.S, ~Glyphosphate_Treatment) 
ddsStone.S$group <- factor(paste(sample_data(Stone.S)$System.loc, sample_data(Stone.S)$Glyphosphate_Treatment))
ddsStone.S$Sampling_date<-factor(ddsStone.S$Sampling_date)
ddsStone.S$Soil_Zone<-factor(ddsStone.S$Soil_Zone)

levels(ddsStone.S$group)<-sub(" ", "", levels(ddsStone.S$group))
design(ddsStone.S) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

#colData(ddsBelt.C)
#?design
#dim(ddsBelt.C)

#update model...
geoMeans = apply(counts(ddsStone.S), 1, gm_mean)

ddsStone.S = estimateSizeFactors(ddsStone.S, geoMeans = geoMeans)
ddsStone.S=estimateDispersions(ddsStone.S)

#Belt.Cvst=getVarianceStabilizedData(ddsBelt.C, method="per-condition", sharingMode="maximum", fitType="local")
Stone.Svst=getVarianceStabilizedData(ddsStone.S) #not working for some reason... Error: is(cds, "CountDataSet") is not TRUE

ddsStone.S<-DESeq(ddsStone.S, test="Wald", fitType="parametric", parallel=T)

resStone.SDS<-results(ddsStone.S)

Stone.Svst<-vst(ddsStone.S, blind=T)
plotPCA(Stone.Svst, intgroup="group")
plotPCA(Stone.Svst, intgroup="System.loc")
plotPCA(Stone.Svst, intgroup="Glyphosphate_Treatment")
plotPCA(Stone.Svst, intgroup="Sampling_date")
plotPCA(Stone.Svst, intgroup="Soil_Zone")
plotPCA(Stone.Svst, intgroup="Glyphosphate_Treatment")

###
###
###

###
### Stoneville corn PCA
###

ddsStone.C = phyloseq_to_deseq2(Stone.C, ~Glyphosphate_Treatment) 
ddsStone.C$group <- factor(paste(sample_data(Stone.C)$System.loc, sample_data(Stone.C)$Glyphosphate_Treatment))
ddsStone.C$Sampling_date<-factor(ddsStone.C$Sampling_date)
ddsStone.C$Soil_Zone<-factor(ddsStone.C$Soil_Zone)

levels(ddsStone.C$group)<-sub(" ", "", levels(ddsStone.C$group))
design(ddsStone.C) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

#colData(ddsBelt.C)
#?design
#dim(ddsBelt.C)

#update model...
geoMeans = apply(counts(ddsStone.C), 1, gm_mean)

ddsStone.C = estimateSizeFactors(ddsStone.C, geoMeans = geoMeans)
ddsStone.C=estimateDispersions(ddsStone.C)

#Belt.Cvst=getVarianceStabilizedData(ddsBelt.C, method="per-condition", sharingMode="maximum", fitType="local")
Stone.Cvst=getVarianceStabilizedData(ddsStone.C) #not working for some reason... Error: is(cds, "CountDataSet") is not TRUE

ddsStone.C<-DESeq(ddsStone.C, test="Wald", fitType="parametric", parallel=T)

resStone.CDS<-results(ddsStone.C)

Stone.Cvst<-vst(ddsStone.C, blind=T)
plotPCA(Stone.Cvst, intgroup="group")
plotPCA(Stone.Cvst, intgroup="System.loc")
plotPCA(Stone.Cvst, intgroup="Glyphosphate_Treatment")
plotPCA(Stone.Cvst, intgroup="Sampling_date")
plotPCA(Stone.Cvst, intgroup="Soil_Zone")
plotPCA(Stone.Cvst, intgroup="Glyphosphate_Treatment")

###
### GLOBAL PCA Soy
###


ddsglobal.S = phyloseq_to_deseq2(global.S, ~Glyphosphate_Treatment) 
ddsglobal.S$group <- factor(paste(sample_data(global.S)$System.loc, sample_data(global.S)$Glyphosphate_Treatment))
ddsglobal.S$Sampling_date<-factor(ddsglobal.S$Sampling_date)
ddsglobal.S$Soil_Zone<-factor(ddsglobal.S$Soil_Zone)

levels(ddsglobal.S$group)<-sub(" ", "", levels(ddsglobal.S$group))
design(ddsglobal.S) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

#colData(ddsBelt.C)
#?design
#dim(ddsBelt.C)

#update model...
geoMeans = apply(counts(ddsglobal.S), 1, gm_mean)

ddsglobal.S = estimateSizeFactors(ddsglobal.S, geoMeans = geoMeans)
ddsglobal.S=estimateDispersions(ddsglobal.S)

#Belt.Cvst=getVarianceStabilizedData(ddsBelt.C, method="per-condition", sharingMode="maximum", fitType="local")
global.Svst=getVarianceStabilizedData(ddsglobal.S) #not working for some reason... Error: is(cds, "CountDataSet") is not TRUE

ddsglobal.S<-DESeq(ddsglobal.S, test="Wald", fitType="parametric", parallel=T)

resglobal.SDS<-results(ddsglobal.S)

global.Svst<-vst(ddsglobal.S, blind=T)
plotPCA(global.Svst, intgroup="group")
plotPCA(global.Svst, intgroup="System.loc")
plotPCA(global.Svst, intgroup="Glyphosphate_Treatment")
plotPCA(global.Svst, intgroup="Sampling_date")
plotPCA(global.Svst, intgroup="Soil_Zone")
plotPCA(global.Svst, intgroup="Glyphosphate_Treatment")


###
### GLOBAL PCA corn
###


ddsglobal.C = phyloseq_to_deseq2(global.C, ~Glyphosphate_Treatment) 
ddsglobal.C$group <- factor(paste(sample_data(global.C)$System.loc, sample_data(global.C)$Glyphosphate_Treatment))
ddsglobal.C$Sampling_date<-factor(ddsglobal.C$Sampling_date)
ddsglobal.C$Soil_Zone<-factor(ddsglobal.C$Soil_Zone)

levels(ddsglobal.C$group)<-sub(" ", "", levels(ddsglobal.C$group))
design(ddsglobal.C) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

#colData(ddsBelt.C)
#?design
#dim(ddsBelt.C)

#update model...
geoMeans = apply(counts(ddsglobal.C), 1, gm_mean)

ddsglobal.C = estimateSizeFactors(ddsglobal.C, geoMeans = geoMeans)
ddsglobal.C=estimateDispersions(ddsglobal.C)

#Belt.Cvst=getVarianceStabilizedData(ddsBelt.C, method="per-condition", sharingMode="maximum", fitType="local")
global.Cvst=getVarianceStabilizedData(ddsglobal.C) #not working for some reason... Error: is(cds, "CountDataSet") is not TRUE

ddsglobal.C<-DESeq(ddsglobal.C, test="Wald", fitType="parametric", parallel=T)

resglobal.CDS<-results(ddsglobal.C)

global.Cvst<-vst(ddsglobal.C, blind=T)
plotPCA(global.Cvst, intgroup="group")
plotPCA(global.Cvst, intgroup="System.loc")
plotPCA(global.Cvst, intgroup="Glyphosphate_Treatment")
plotPCA(global.Cvst, intgroup="Sampling_date")
plotPCA(global.Cvst, intgroup="Soil_Zone")
plotPCA(global.Cvst, intgroup="Glyphosphate_Treatment")


###
###
### Rarefied PCAs
###
###

rareBacarc<-readRDS("~/Desktop/PhD/PHP/bacarc_RAREcomDat.rds")


Rglobal.S<-subset_samples(rareBacarc, crop=="soy")
Rglobal.C<-subset_samples(rareBacarc, crop=="corn")


#Beltsville subset without rarefaction
RBeltsville<-subset_samples(rareBacarc, Location=="Beltsville")
#table(sample_data(Beltsville)$crop) #check to make sure formating is good...
RBelt.C<-subset_samples(RBeltsville, crop=="corn")
RBelt.S<-subset_samples(RBeltsville, crop=="soy")

#Stoneville subset without rarefaction
RStoneville<-subset_samples(rareBacarc, Location=="Stoneville")
RStone.S<-subset_samples(RStoneville, crop=="soy")
RStone.C<-subset_samples(RStoneville, crop=="corn")

#head(sample_data(Belt.C))

# remove non-round-up ready
RBelt.C<-subset_samples(RBelt.C, genotype=="RR")
RBelt.S<-subset_samples(RBelt.S, genotype=="RR")

# remove non-round-up ready
RStone.C<-subset_samples(RStone.C, genotype=="RR")
RStone.S<-subset_samples(RStone.S, genotype=="RR")

###
###
### Start PCA!
###
###

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

library(DESeq2)
### Beltsville Corn PCA
head(sample_data(Belt.C))
ddsRBelt.C = phyloseq_to_deseq2(RBelt.C, ~Glyphosphate_Treatment) 
ddsRBelt.C$group <- factor(paste(sample_data(RBelt.C)$System.loc, sample_data(RBelt.C)$Glyphosphate_Treatment))
ddsRBelt.C$Sampling_date<-factor(ddsRBelt.C$Sampling_date)
ddsRBelt.C$Soil_Zone<-factor(ddsRBelt.C$Soil_Zone)

levels(ddsRBelt.C$group)<-sub(" ", "", levels(ddsRBelt.C$group))
design(ddsRBelt.C) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

#colData(ddsBelt.C)
#?design
#dim(ddsBelt.C)

#update model...
geoMeans = apply(counts(ddsRBelt.C), 1, gm_mean)

ddsRBelt.C = estimateSizeFactors(ddsRBelt.C, geoMeans = geoMeans)
ddsRBelt.C=estimateDispersions(ddsRBelt.C)

#Belt.Cvst=getVarianceStabilizedData(ddsBelt.C, method="per-condition", sharingMode="maximum", fitType="local")
RBelt.Cvst=getVarianceStabilizedData(ddsRBelt.C) #not working for some reason... Error: is(cds, "CountDataSet") is not TRUE

ddsRBelt.C<-DESeq(ddsRBelt.C, test="Wald", fitType="parametric", parallel=T)

resRBelt.CDS<-results(ddsRBelt.C)

RBelt.Cvst<-vst(ddsRBelt.C, blind=T)
plotPCA(RBelt.Cvst, intgroup="group")
plotPCA(RBelt.Cvst, intgroup="System.loc")
plotPCA(RBelt.Cvst, intgroup="Glyphosphate_Treatment")
plotPCA(RBelt.Cvst, intgroup="Sampling_date")
plotPCA(RBelt.Cvst, intgroup="Soil_Zone")
plotPCA(RBelt.Cvst, intgroup="Glyphosphate_Treatment")


###
### Beltsville Soy PCA
###


ddsRBelt.S = phyloseq_to_deseq2(RBelt.S, ~Glyphosphate_Treatment) 
ddsRBelt.S$group <- factor(paste(sample_data(RBelt.S)$System.loc, sample_data(RBelt.S)$Glyphosphate_Treatment))
ddsRBelt.S$Sampling_date<-factor(ddsRBelt.S$Sampling_date)
ddsRBelt.S$Soil_Zone<-factor(ddsRBelt.S$Soil_Zone)

levels(ddsRBelt.S$group)<-sub(" ", "", levels(ddsRBelt.S$group))
design(ddsRBelt.S) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

#colData(ddsBelt.C)
#?design
#dim(ddsBelt.C)

#update model...
geoMeans = apply(counts(ddsRBelt.S), 1, gm_mean)

ddsRBelt.S = estimateSizeFactors(ddsRBelt.S, geoMeans = geoMeans)
ddsRBelt.S=estimateDispersions(ddsRBelt.S)

#Belt.Cvst=getVarianceStabilizedData(ddsBelt.C, method="per-condition", sharingMode="maximum", fitType="local")
RBelt.Svst=getVarianceStabilizedData(ddsRBelt.S) #not working for some reason... Error: is(cds, "CountDataSet") is not TRUE

ddsRBelt.S<-DESeq(ddsRBelt.S, test="Wald", fitType="parametric", parallel=T)

resRBelt.SDS<-results(ddsRBelt.S)

RBelt.Svst<-vst(ddsRBelt.S, blind=T)
plotPCA(RBelt.Svst, intgroup="group")
plotPCA(RBelt.Svst, intgroup="System.loc")
plotPCA(RBelt.Svst, intgroup="Glyphosphate_Treatment")
plotPCA(RBelt.Svst, intgroup="Sampling_date")
plotPCA(RBelt.Svst, intgroup="Soil_Zone")
plotPCA(RBelt.Svst, intgroup="Glyphosphate_Treatment")

###
### Stoneville Soy PCA
###

ddsRStone.S = phyloseq_to_deseq2(RStone.S, ~Glyphosphate_Treatment) 
ddsRStone.S$group <- factor(paste(sample_data(RStone.S)$System.loc, sample_data(RStone.S)$Glyphosphate_Treatment))
ddsRStone.S$Sampling_date<-factor(ddsRStone.S$Sampling_date)
ddsRStone.S$Soil_Zone<-factor(ddsRStone.S$Soil_Zone)

levels(ddsRStone.S$group)<-sub(" ", "", levels(ddsRStone.S$group))
design(ddsRStone.S) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

#colData(ddsBelt.C)
#?design
#dim(ddsBelt.C)

#update model...
geoMeans = apply(counts(ddsRStone.S), 1, gm_mean)

ddsRStone.S = estimateSizeFactors(ddsRStone.S, geoMeans = geoMeans)
ddsRStone.S=estimateDispersions(ddsRStone.S)

#Belt.Cvst=getVarianceStabilizedData(ddsBelt.C, method="per-condition", sharingMode="maximum", fitType="local")
RStone.Svst=getVarianceStabilizedData(ddsRStone.S) #not working for some reason... Error: is(cds, "CountDataSet") is not TRUE

ddsRStone.S<-DESeq(ddsRStone.S, test="Wald", fitType="parametric", parallel=T)

resRStone.SDS<-results(ddsRStone.S)

RStone.Svst<-vst(ddsRStone.S, blind=T)
plotPCA(RStone.Svst, intgroup="group")
plotPCA(RStone.Svst, intgroup="System.loc")
plotPCA(RStone.Svst, intgroup="Glyphosphate_Treatment")
plotPCA(RStone.Svst, intgroup="Sampling_date")
plotPCA(RStone.Svst, intgroup="Soil_Zone")
plotPCA(RStone.Svst, intgroup="Glyphosphate_Treatment")

###
###
###

###
### Stoneville corn PCA
###

ddsRStone.C = phyloseq_to_deseq2(RStone.C, ~Glyphosphate_Treatment) 
ddsRStone.C$group <- factor(paste(sample_data(RStone.C)$System.loc, sample_data(RStone.C)$Glyphosphate_Treatment))
ddsRStone.C$Sampling_date<-factor(ddsRStone.C$Sampling_date)
ddsRStone.C$Soil_Zone<-factor(ddsRStone.C$Soil_Zone)

levels(ddsRStone.C$group)<-sub(" ", "", levels(ddsRStone.C$group))
design(ddsRStone.C) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

#colData(ddsBelt.C)
#?design
#dim(ddsBelt.C)

#update model...
geoMeans = apply(counts(ddsRStone.C), 1, gm_mean)

ddsRStone.C = estimateSizeFactors(ddsRStone.C, geoMeans = geoMeans)
ddsRStone.C=estimateDispersions(ddsRStone.C)

#Belt.Cvst=getVarianceStabilizedData(ddsBelt.C, method="per-condition", sharingMode="maximum", fitType="local")
RStone.Cvst=getVarianceStabilizedData(ddsRStone.C) #not working for some reason... Error: is(cds, "CountDataSet") is not TRUE

ddsRStone.C<-DESeq(ddsRStone.C, test="Wald", fitType="parametric", parallel=T)

resRStone.CDS<-results(ddsRStone.C)

RStone.Cvst<-vst(ddsRStone.C, blind=T)
plotPCA(RStone.Cvst, intgroup="group")
plotPCA(RStone.Cvst, intgroup="System.loc")
plotPCA(RStone.Cvst, intgroup="Glyphosphate_Treatment")
plotPCA(RStone.Cvst, intgroup="Sampling_date")
plotPCA(RStone.Cvst, intgroup="Soil_Zone")
plotPCA(RStone.Cvst, intgroup="Glyphosphate_Treatment")

###
### GLOBAL PCA Soy
###


ddsRglobal.S = phyloseq_to_deseq2(Rglobal.S, ~Glyphosphate_Treatment) 
ddsRglobal.S$group <- factor(paste(sample_data(Rglobal.S)$System.loc, sample_data(Rglobal.S)$Glyphosphate_Treatment))
ddsRglobal.S$Sampling_date<-factor(ddsRglobal.S$Sampling_date)
ddsRglobal.S$Soil_Zone<-factor(ddsRglobal.S$Soil_Zone)

levels(ddsRglobal.S$group)<-sub(" ", "", levels(ddsRglobal.S$group))
design(ddsRglobal.S) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

#colData(ddsBelt.C)
#?design
#dim(ddsBelt.C)

#update model...
geoMeans = apply(counts(ddsRglobal.S), 1, gm_mean)

ddsRglobal.S = estimateSizeFactors(ddsRglobal.S, geoMeans = geoMeans)
ddsRglobal.S=estimateDispersions(ddsRglobal.S)

#Belt.Cvst=getVarianceStabilizedData(ddsBelt.C, method="per-condition", sharingMode="maximum", fitType="local")
Rglobal.Svst=getVarianceStabilizedData(ddsRglobal.S) #not working for some reason... Error: is(cds, "CountDataSet") is not TRUE

ddsRglobal.S<-DESeq(ddsRglobal.S, test="Wald", fitType="parametric", parallel=T)

resRglobal.SDS<-results(ddsRglobal.S)

Rglobal.Svst<-vst(ddsRglobal.S, blind=T)
plotPCA(Rglobal.Svst, intgroup="group")
plotPCA(Rglobal.Svst, intgroup="System.loc")
plotPCA(Rglobal.Svst, intgroup="Glyphosphate_Treatment")
plotPCA(Rglobal.Svst, intgroup="Sampling_date")
plotPCA(Rglobal.Svst, intgroup="Soil_Zone")
plotPCA(Rglobal.Svst, intgroup="Glyphosphate_Treatment")


###
### GLOBAL PCA corn
###


ddsRglobal.C = phyloseq_to_deseq2(Rglobal.C, ~Glyphosphate_Treatment) 
ddsRglobal.C$group <- factor(paste(sample_data(Rglobal.C)$System.loc, sample_data(Rglobal.C)$Glyphosphate_Treatment))
ddsRglobal.C$Sampling_date<-factor(ddsRglobal.C$Sampling_date)
ddsRglobal.C$Soil_Zone<-factor(ddsRglobal.C$Soil_Zone)

levels(ddsRglobal.C$group)<-sub(" ", "", levels(ddsRglobal.C$group))
design(ddsRglobal.C) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

#colData(ddsBelt.C)
#?design
#dim(ddsBelt.C)

#update model...
geoMeans = apply(counts(ddsRglobal.C), 1, gm_mean)

ddsRglobal.C = estimateSizeFactors(ddsRglobal.C, geoMeans = geoMeans)
ddsRglobal.C=estimateDispersions(ddsRglobal.C)

#Belt.Cvst=getVarianceStabilizedData(ddsBelt.C, method="per-condition", sharingMode="maximum", fitType="local")
Rglobal.Cvst=getVarianceStabilizedData(ddsRglobal.C) #not working for some reason... Error: is(cds, "CountDataSet") is not TRUE

ddsRglobal.C<-DESeq(ddsRglobal.C, test="Wald", fitType="parametric", parallel=T)

resRglobal.CDS<-results(ddsRglobal.C)

Rglobal.Cvst<-vst(ddsRglobal.C, blind=T)
plotPCA(Rglobal.Cvst, intgroup="group")
plotPCA(Rglobal.Cvst, intgroup="System.loc")
plotPCA(Rglobal.Cvst, intgroup="Glyphosphate_Treatment")
plotPCA(Rglobal.Cvst, intgroup="Sampling_date")
plotPCA(Rglobal.Cvst, intgroup="Soil_Zone")
plotPCA(Rglobal.Cvst, intgroup="Glyphosphate_Treatment")
?varianceStabilizingTransformation






###
###
### ADONIS ###
###
###

# the strata argument constrains the permutation analysis;
# therefore statistical comparisons cannot be made of strata categories
# e.g. if strata = system.loc, then systems cannot be compared
# thus: adonis(bray~Sampling_date+Glyphosphate_Treatment+Loc_plot_ID+crop+Soil_Zone, bacarc, strata="System.loc", permutations=999, method="bray")
# constrains the permutation to give better resolution within the strata

library(vegan)

?adonis
