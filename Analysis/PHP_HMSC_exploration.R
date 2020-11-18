# 
install.packages("Hmsc")
library(phyloseq)
library(Hmsc)
library(corrplot)
library(dada2)
library(phangorn)
library(msa)

getwd()
setwd("~/project/php-bacarc/new_preproc/")

fun<-readRDS("FUN_RAWcomDat2020.rds") # load fungal data
#fun<-readRDS("~/Desktop/PhD/PHP/PHP_Fungi_2020.rds")
bac<-readRDS("bacarc_RAWcomDat.rds") # load bacterial data
#bac<-readRDS("~/Desktop/PhD/PHP/bacarc_RAWcomDat.rds")

# normalize bacteria / subset to taxa of interest
bac.filt = prune_samples(sample_sums(bac)>=10000, bac)
sample_data(bac.filt)$Quantity..picograms.<-as.numeric(sample_data(bac.filt)$Quantity..picograms.)
bac.filt = prune_samples(!is.na(sample_data(bac.filt)$Quantity..picograms.), bac.filt)
bac.filt = prune_samples(sample_data(bac.filt)$Quantity..picograms.!=0, bac.filt)
meta<-as.data.frame(as.matrix(sample_data(bac.filt)))
meta$Quantity..picograms.<-as.numeric(meta$Quantity..picograms.)
normfactor<-2000000*(meta$Quantity..picograms./mean(meta$Quantity..picograms.))



#taxa_names(bac.filt)<-paste0("taxon", c(1:ntaxa(bac.filt)))
df<-data.frame(t(as.matrix(otu_table(transform_sample_counts(bac.filt, function(x) x/sum(x))))))
normfactor<-2000000*(meta$Quantity..picograms./mean(meta$Quantity..picograms.))
n<-normfactor
scaledf<-data.frame(mapply(`*`, df, n))
t.bac.filt<-bac.filt
scaledf<-round(scaledf)
rownames(scaledf)<-paste0("taxon", rownames(scaledf))
sample_names(t.bac.filt)<-paste0("X", sample_names(t.bac.filt))
otu_table(t.bac.filt)<-otu_table(scaledf, taxa_are_rows=TRUE) # woohoo!!
proteobacteria.PHP<-subset_taxa(t.bac.filt, Phylum=="Proteobacteria")
proteobacteria.PHP<-subset_taxa(t.bac.filt, Family!="Mitochondria")
Proteobacteria.Org6<-prune_samples(sample_data(proteobacteria.PHP)$System.loc=="Org_6", proteobacteria.PHP)
Proteobacteria.Org6<-tax_glom(Proteobacteria.Org6, taxrank="Genus", NArm=FALSE)
Proteobacteria.Org6.corn<-prune_samples(sample_data(Proteobacteria.Org6)$crop=="corn", Proteobacteria.Org6)
Proteobacteria.Org6.corn<-filter_taxa(Proteobacteria.Org6.corn, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
Proteobacteria.Org6.soy<-prune_samples(sample_data(Proteobacteria.Org6)$crop=="soy", Proteobacteria.Org6)
Proteobacteria.Org6.soy<-filter_taxa(Proteobacteria.Org6.soy, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

# normalize fungi / subset to taxa of interest
fun.filt<-prune_samples(sample_sums(fun)>=10000, fun)
fun.filt<-prune_samples() # match samples to bacterial dataset
# can't do relative abundance normalization; Do CSS normalization instead?
# or total DNA in each sample?
Agaricomycetes<-subset_taxa(fun.filt, Class=="Agaricomycetes")
Agaricomycetes.Org6<-prune_samples(sample_data(Agaricomycetes)$System.loc=="Org_6", Agaricomycetes)
Sordariomycetes<-subset_taxa(fun.filt, Class=="Sordariomycetes")
Sordariomycetes.Org6<-prune_samples(sample_data(Sordariomycetes)$System.loc=="Org_6", Sordariomycetes)

Sordariomycetes.Org6.corn<-prune_samples(sample_data(Sordariomycetes.Org6)$crop=="corn", Sordariomycetes.Org6)
Sordariomycetes.Org6.soy<-prune_samples(sample_data(Sordariomycetes.Org6)$crop=="soy", Sordariomycetes.Org6)

Agaricomycetes.Org6.corn<-prune_samples(sample_data(Agaricomycetes.Org6)$crop=="corn", Sordariomycetes.Org6)
Agaricomycetes.Org6.soy<-prune_samples(sample_data(Agaricomycetes.Org6)$crop=="soy", Sordariomycetes.Org6)

# Hmsc for bacteria alone
b.Ydat<-as.data.frame(as.matrix(otu_table(Proteobacteria.Org6.corn)))# dim(b.Ydat)
b.XDat1<-as.data.frame(as.matrix(sample_data(Proteobacteria.Org6.corn)), stringsAsFactors = TRUE)
b.XDat1<-b.XDat1[,c(1,4,5,6,8,9,10,11,12,14)]
rownames(b.XDat1)<-c(1:nrow(b.XDat1))
b.XDat1$Sample<-as.factor(c(1:nrow(b.XDat1)))
b.XFormula= ~Glyphosphate_Treatment + genotype # changed since last time
studyDesign = data.frame("Sample"=b.XDat1$Sample,"Loc_plot_ID"=b.XDat1$Loc_plot_ID, "Sampling_date"=b.XDat1$Sampling_date)
rL1 <- HmscRandomLevel(units=studyDesign$Sample)
rL2 <- HmscRandomLevel(units=studyDesign$Loc_plot_ID) # set random level to loc_plot_ID
rL3 <- HmscRandomLevel(units=studyDesign$Sampling_date) # second random level is sampling 
# final model should have pH and organic matter as random levels
# first do probit model ####
b.m1<-Hmsc(Y=as.matrix(t(b.Ydat)), XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="probit")
b.m1<-sampleMcmc(b.m.lp, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)

# examine correlation matrix for probit model
OmegaCor.probit=computeAssociations(b.m1)

hist(OmegaCor.probit[[1]]$support)
hist(OmegaCor.probit[[2]]$support)
hist(OmegaCor.probit[[3]]$support)

supportLevel=0.95
toPlot.probit=((OmegaCor.probit[[1]]$support>supportLevel) +
           (OmegaCor.probit[[1]]$support<(1-supportLevel))>0)*OmegaCor.probit[[1]]$mean

corrplot(toPlot.probit, method="color", 
         col=colorRampPalette(c("red", "white", "blue"))(200),
         tl.cex=0.6, tl.col="black",
         title=paste("Probit random effect level:", b.m1$rLNames[1]), mar=c(0,0,1,0))

# convert to conditional matrix
Y2<-as.matrix(t(b.Ydat))
Y2[Y2==0]<-NA

# run full model conditional on presence ####
b.m.lp<-Hmsc(Y=Y2, XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="lognormal poisson")
b.m.p<-Hmsc(Y=Y2, XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="poisson")


#b.m.lp<-Hmsc(Y=Y2, XData=b.XDat, XFormula=b.XFormula, distr="lognormal poisson")
#b.m.p<-Hmsc(Y=Y2, XData=b.XDat, XFormula=b.XFormula, distr="poisson")
b.m2.lp<-sampleMcmc(b.m.lp, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100) # lognormal poisson distribution
b.m2.p<-sampleMcmc(b.m.p, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100) # poisson distribution test

# evaluate model fits

preds.lp=computePredictedValues(b.m2.lp, expected = FALSE)
preds.p=computePredictedValues(b.m2.p, expected = FALSE)


evaluateModelFit(hM=b.m2.lp, predY=preds.lp)
evaluateModelFit(hM=b.m2.p, predY=preds.p)


mpost<-convertToCodaObject(b.m2)
#par(mfrow=c(1,1))
ess.beta=effectiveSize(mpost$Beta)
psrf.beta=gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
hist(ess.beta)
hist(psrf.beta)

preds=computePredictedValues(b.m2, expected = TRUE)
MF=evaluateModelFit(hM=b.m2, predY=preds)
hist(MF$SR2, xlim=c(0,1), main=paste0("mean =", round(mean(MF$SR2),2))) # Glyphosate has very little explanatory power on microbial abundance ...

hist(MF$C.SR2, main=paste0("mean =", round(mean(MF$C.SR2),2)))

head(b.m2$X)
VP = computeVariancePartitioning(b.m2, group=c(1,1,2), groupnames=c("Glyphosphate_treatment", "genotype"))
plotVariancePartitioning(b.m2, VP=VP)

postBeta=getPostEstimate(b.m2, parName="Beta")
plotBeta(b.m2, post=postBeta, param="Support", plotTree=F, supportLevel=0.95, split=0.4, spNamesNumbers=c(F,F))
plotBeta(b.m2, post=postBeta, param="Sign", plotTree=F, supportLevel=0.95, split=0.4, spNamesNumbers=c(F,F))
plotBeta(b.m2, post=postBeta, param="Mean", plotTree=F, supportLevel=0.95, split=0.4, spNamesNumbers=c(F,F))

hist(postBeta$mean[2,])
# let's look at seq counts at each level of glyphosate:
as.data.frame(t(as.matrix(otu_table(Proteobacteria.Org6.corn))))%>%group_by(c(b.XDat1$Plot_Loc_ID,b.XDat1$genotype,b.XDat1$Soil_Zone,b.XDat1$Glyphosphate_Treatment,b.XDat1$Sampling_date)) %>% summarise("pre"=mean(pre), "post"=mean(post))# %>% summarise("effect ratio"=pre/post)
b.XDat1$effect<-cumsum(!duplicated(b.XDat1[,c(4:6,9),])) # sampling date is no 3 ..

yd<-as.data.frame(as.matrix(otu_table(Proteobacteria.Org6.corn)))
yd<-yd[,b.XDat1$genotype!="Non_RR"]
xd<-b.XDat1[b.XDat1$genotype!="Non_RR",]
xd %>% 
  mutate(effect = group_indices(xd, .dots=c("Loc_plot_ID", "Soil_Zone"))) -> xd
xd %>% 
  mutate(ID = group_indices(xd, .dots=c("Loc_plot_ID", "Soil_Zone", "Glyphosphate_Treatment"))) -> xd

View(yd)
effectdf<-matrix(data=NA, nrow=nrow(yd),ncol=length(unique(xd$ID)))
yd[yd==0]<-1

(yd[,b.XDat1$ID==1 & b.XDat1$Sampling_date=="pre"]-yd[,b.XDat1$ID==1 & b.XDat1$Sampling_date=="post"])/yd[,b.XDat1$ID==1 & b.XDat1$Sampling_date=="pre"]

for(i in 1:length(unique(xd$ID))){
  if(length(xd$ID[xd$ID==i])!=2){effectdf[,i]<-c(rep(NA, 316))}
  if(length(xd$ID[xd$ID==i])==2){
  effectdf[,i]<-mean(yd[,xd$effect==i & xd$Sampling_date=="pre"]-yd[,xd$effect==i & xd$Sampling_date=="post"])/mean(yd[,xd$effect==i & xd$Glyphosphate_Treatment=="pre"]-yd[,xd$effect==i & xd$Sampling_date=="post"])#/yd[,xd$ID==i & xd$Sampling_date=="pre"]}
  
} # this function effectively sorts bc it goes in order from low to high
l<-c(1:length(unique(xd$ID)))
e<-l
for(i in 1:length(l)){e[i]<-unique(xd$effect[xd$ID==i])}
e2<-l
for(i in 1:length(l)){e2[i]<-unique(xd$Glyphosphate_Treatment[xd$ID==i])} # 1 = no spray; 2 = spray

xd$Glyphosphate_Treatment[xd$ID==2]
effectdf2<-matrix(data=NA, nrow=ntaxa(Proteobacteria.Org6.corn), ncol=15)
for(i in 1:ncol(effectdf2)){
  if(length(e2[e==i])==2){
  effectdf2[,i]<-(effectdf[,e==i & e2==1]/effectdf[,e==i & e2==2])
  }
  print(i)
}
View(effectdf)
View(effectdf2)


# now run with genotype
studyDesign = data.frame("Sample"=b.XDat1$Sample,"Loc_plot_ID"=b.XDat1$Loc_plot_ID, "Sampling_date"=b.XDat1$Sampling_date)

# species interactions
OmegaCor=computeAssociations(b.m2)

hist(OmegaCor[[1]]$support)
hist(OmegaCor[[2]]$support)
hist(OmegaCor[[3]]$support)

supportLevel=0.95
toPlot1=((OmegaCor[[1]]$support>supportLevel) +
          (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
toPlot2=((OmegaCor[[2]]$support>supportLevel) +
           (OmegaCor[[2]]$support<(1-supportLevel))>0)*OmegaCor[[2]]$mean
toPlot3=((OmegaCor[[3]]$support>supportLevel) +
           (OmegaCor[[3]]$support<(1-supportLevel))>0)*OmegaCor[[3]]$mean

corrplot(toPlot1, method="color", 
         col=colorRampPalette(c("red", "white", "blue"))(200),
         tl.cex=0.6, tl.col="black",
         title=paste("random effect level:", b.m2$rLNames[1]), mar=c(0,0,1,0))
corrplot(toPlot2, method="color", 
         col=colorRampPalette(c("red", "white", "blue"))(200),
         tl.cex=0.6, tl.col="black",
         title=paste("random effect level:", b.m2$rLNames[2]), mar=c(0,0,1,0))
corrplot(toPlot3, method="color", 
         col=colorRampPalette(c("red", "white", "blue"))(200),
         tl.cex=0.6, tl.col="black",
         title=paste("random effect level:", b.m2$rLNames[3]), mar=c(0,0,1,0))

# make glyphosate a random variable
b.XFormula= ~genotype # changed since last time
studyDesign = data.frame("Sample"=b.XDat1$Sample,"Loc_plot_ID"=b.XDat1$Loc_plot_ID, "Sampling_date"=b.XDat1$Sampling_date, "Glyphosphate_Treatment"=b.XDat1$Glyphosphate_Treatment)
rL1 <- HmscRandomLevel(units=studyDesign$Sample)
rL2 <- HmscRandomLevel(units=studyDesign$Loc_plot_ID) # set random level to loc_plot_ID
rL3 <- HmscRandomLevel(units=studyDesign$Sampling_date) # second random level is sampling 
rL4 <- HmscRandomLevel(units=studyDesign$Glyphosphate_Treatment)
# final model should have pH and organic matter as random levels
b.m.lp2<-Hmsc(Y=Y2, XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3, "Glyphosphate_Treatment"=rL4), distr="lognormal poisson")


#b.m.lp<-Hmsc(Y=Y2, XData=b.XDat, XFormula=b.XFormula, distr="lognormal poisson")
#b.m.p<-Hmsc(Y=Y2, XData=b.XDat, XFormula=b.XFormula, distr="poisson")
b.m2.lp2<-sampleMcmc(b.m.lp2, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100) # lognormal poisson distribution

OmegaCor.gly=computeAssociations(b.m2.lp2)
toPlot1.gly=((OmegaCor.gly[[4]]$support>supportLevel) +
           (OmegaCor.gly[[4]]$support<(1-supportLevel))>0)*OmegaCor.gly[[4]]$mean
corrplot(toPlot1.gly, method="color", 
         col=colorRampPalette(c("red", "white", "blue"))(200),
         tl.cex=0.6, tl.col="black",
         title=paste("random effect level:", b.m2.lp2$rLNames[4]), mar=c(0,0,1,0))


# Scratch ####

# make phylogenetic tree ####
seqs <- taxa_names(bac)
names(seqs) <- seqs # This propagates to the tip labels of the tree
mult <- msa(seqs, method="ClustalW", type="dna", order="input")
library("phangorn")
phang.align <- as.phyDat(mult, type="DNA", names=names(seqs))
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))