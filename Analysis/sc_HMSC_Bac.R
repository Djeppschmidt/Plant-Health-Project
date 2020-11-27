# scinet R script for analysis
# run computationally challenging functions
# output RDS files that are easy to use to build outputs

# load packages for analysis

library(Hmsc)
library(phyloseq)
library(igraph)
library(stringr)
library(Model.Microbiome)

# To Do List:
# return to see table merging for the metadata      [_] ***
        # filter archaea / mitochondria             [_]
        # data for all samples                      [_]
# sort code into proper files                       [X]
# Q1 Probit model                                   [X]
# Q1 Lognormal Model                                [X]
# Q1 (For Q2) Lognormal Model - Farming System      [X]
# Q1 (For Q2) Probit model - Farming System         [X]
# Q1 (For Q3) Lognormal Model + AM fungi            [X]
# Q1 (For Q3) Probit Model + AM fungi               [X]
  
# Q2 Probit for each plot                           [X]
# Q2 lognormal for each plot                        [X]
# Q2 Probit Model + AM Fungi for each plot (Q3)     [X]
# Q2 lognormal Model + AM Fungi for each plot (Q3)  [X]

# Q3 Predict Env. Gradient (Global)                 [_] # can do this in post
# Q3 Compute variation explained by AM (Global)     [X]

# Final Syntax Check                                [_] 
# pilot run functions on local computer             [_] ***
# move visualization script to visfile              [X]


# Analysis 1: Role of AMF abundance on microbial networks; Drivers of +/- interaciton networks ; Exploratory analysis of bacterial and fungal competition or symbiosis
    # use combined dataset
fb.RH<-readRDS("Data/fb_Rhizosphere_2020.RDS")
fb.RH<-tax_glom(fb.RH, taxrank="Genus")
    # Model 1: All edaphic factors as predictors (date/time as random)
    # Apply the model to whole dataset
Ydat<-as.data.frame(as.matrix(otu_table(fb.RH)))
XData<-as.data.frame(as.matrix(sample_data(fb.RH)), stringsAsFactors = TRUE)
XDat1<-XData[,c(2,3,6,8:10,15,47,75,81)] # subset data to what I need for this run!!
rownames(XDat1)<-c(1:nrow(XDat1))
XDat1$Sample<-as.factor(c(1:nrow(XDat1)))
XFormula= ~Glyphosphate_Treatment + genotype + System.loc + crop + pH + OM...
studyDesign = data.frame("Sample"=XDat1$Sample,"Sampling_date"=XDat1$Sampling_date, "year"=XDat1$year, "Plot"=XDat1$Loc_plot_ID)
rL1 <- HmscRandomLevel(units=studyDesign$Sample)
rL2 <- HmscRandomLevel(units=studyDesign$Loc_plot_ID)
rL3 <- HmscRandomLevel(units=studyDesign$Sampling_date) # set random level to loc_plot_ID
rL4 <- HmscRandomLevel(units=studyDesign$year) # second random level is sampling 
# final model should have pH and organic matter as random levels
# first do probit model ####
Q1.bfm1<-Hmsc(Y=as.matrix(t(Ydat)), XData=XDat1, XFormula=XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Sampling_date"=rL2, "year"=rL3), distr="probit") # chech XDat1 to make sure all data is there before running!
probitQ1<-NULL
probitQ1$bf<-sampleMcmc(Q1.bfm1, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)

# examine correlation matrix for probit model
# get outputs for HMSC analysis
probitQ1$OmegaCor.probit=computeAssociations(probitQ1$bf)
probitQ1$preds<-computePredictedValues(probitQ1$bf)
probitQ1$MF<-evaluateModelFit(probitQ1$bf, predY=probitQ1$preds)
probitQ1$VP<-computeVariancePartitioning(probit1$bf, group=c(1,1,1,2,1,3,3,4,4,4,4), groupnames=c("Expt", "System", "Edaphic", "Random"))
probitQ1$postBeta<-getPostEstimate(probit1$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal

saveRDS(probitQ1, "Data/ProbitModel1Dat.RDS") # plot on local machine

#probit1$Gradient.pH<-constructGradient()
#robit1$Gradient.OM<-constructGradient()
#probit1$predpHGradient<-predict()
#probit1$predOMGradient<-predict()
G.Y2<-as.matrix(t(Ydat)) # global OTU table
G.Y2[G.Y2==0]<-NA # remove zero counts

LogPoiQ1<-NULL
LogPoiQ1$bf<-Hmsc(Y=G.Y2, XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="poisson") # use same formulas as previous
LogPoiQ1$bf<-sampleMcmc(b.m.lp, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)

LogPoiQ1$OmegaCor.lp=computeAssociations(LogPoiQ1$bf)
LogPoiQ1$preds<-computePredictedValues(LogPoiQ1$bf)
LogPoiQ1$MF<-evaluateModelFit(LogPoiQ1$bf, predY=LogPoiQ1$preds)
LogPoiQ1$VP<-computeVariancePartitioning(LogPoiQ1$bf, group=c(1,1,1,2,1,3,3,4,4,4,4), groupnames=c("Expt", "System", "Edaphic", "Random"))
LogPoiQ1$postBeta<-getPostEstimate(LogPoiQ1$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal

saveRDS(LogPoiQ1, "Data/LogPoiModel1Dat.RDS")

# Q2: effect of farming system: same as above but model removes Farming system 

GQ2.XFormula= ~Glyphosphate_Treatment + genotype + crop + pH + OM...
# first do probit model ####
GQ2.bf<-Hmsc(Y=Ydat, XData=XDat, XFormula=GQ2.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Sampling_date"=rL2, "year"=rL3), distr="probit")
probitGQ2<-NULL
probitGQ2$bf<-sampleMcmc(GQ2.bf, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)

# examine correlation matrix for probit model
# get outputs for HMSC analysis
probitGQ2$OmegaCor.probit=computeAssociations(probitGQ2$bf)
probitGQ2$preds<-computePredictedValues(probitGQ2$bf)
probitGQ2$MF<-evaluateModelFit(probitGQ2$bf, predY=probitGQ2$preds)
probitGQ2$VP<-computeVariancePartitioning(probitGQ2$bf, group=c(1,1,1,1,2,2,3,3,3,3), groupnames=c("Expt", "Edaphic", "Random"))
probitGQ2$postBeta<-getPostEstimate(probitGQ2$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal

saveRDS(probitGQ2, "Data/ProbitModel2_GDat.RDS") # plot on local machine

LogPoiQ2<-NULL
LogPoiQ2$bf<-Hmsc(Y=G.Y2, XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="poisson") # use same formulas as previous
LogPoiGQ2$bf<-sampleMcmc(LogPoiQ2$bf, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)

LogPoiGQ2$OmegaCor.lp=computeAssociations(LogPoiGQ2$bf.lp)
LogPoiGQ2$preds<-computePredictedValues(LogPoiGQ2$bf.lp)
LogPoiGQ2$MF<-evaluateModelFit(LogPoiGQ2$bf.lp, predY=LogPoiGQ2$preds)
LogPoiGQ2$VP<-computeVariancePartitioning(LogPoiGQ2$bf.lp, group=c(1,1,1,1,2,2,3,3,3,3), groupnames=c("Expt", "Edaphic", "Random"))
LogPoiGQ2$postBeta<-getPostEstimate(LogPoiGQ2$bf.lp, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal

saveRDS(LogPoiGQ2, "Data/LogPoiModel2_GDat.RDS")

# Q3
GQ3.XFormula= ~Glyphosphate_Treatment + genotype + System.loc + crop + pH + OM... + NLFA.AM.Fungi

probitGQ3<-NULL
probitGQ3$bf<-Hmsc(Y=Ydat, XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Sampling_date"=rL2, "year"=rL3), distr="probit")
probitGQ3$bf<-sampleMcmc(probitGQ3$bf, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)

# examine correlation matrix for probit model
# get outputs for HMSC analysis
probitGQ3$OmegaCor.probit=computeAssociations(probitGQ3$bf)
probitGQ3$preds<-computePredictedValues(probitGQ3$bf)
probitGQ3$MF<-evaluateModelFit(probitGQ3$bf, predY=probitQ1$preds)
probitGQ3$VP<-computeVariancePartitioning(probitGQ3$bf, group=c(1,1,1,1,1,2,2,3,4,4,4,4), groupnames=c("Expt", "Edaphic","AMF", "Random"))
probitGQ3$postBeta<-getPostEstimate(probitGQ3$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal

saveRDS(probitGQ3, "Data/ProbitModel3_GDat.RDS") # plot on local machine

LogPoiGQ3<-NULL
LogPoiQ3$bf<-Hmsc(Y=G.Y2, XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="poisson") # use same formulas as previous
LogPoiGQ3$bf<-sampleMcmc(LogPoiGQ3$bf, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100) # check on previous sampleMCMC commands !!!!

#### Edit to make work for logpoisson model !!!!
LogPoiGQ3$OmegaCor.lp=computeAssociations(LogPoiGQ3$bf)
LogPoiGQ3$preds<-computePredictedValues(LogPoiGQ3$bf)
LogPoiGQ3$MF<-evaluateModelFit(LogPoiGQ3$bf, predY=LogPoiGQ3$preds)
LogPoiGQ3$VP<-computeVariancePartitioning(LogPoiGQ3$bf, group=c(1,1,1,1,1,2,2,3,4,4,4,4), groupnames=c("Expt", "Edaphic", "AMF", "Random"))
LogPoiGQ3$postBeta<-getPostEstimate(LogPoiGQ3$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal
saveRDS(LogPoiGQ3, "Data/LogPoiModel3_Gdat.RDS")



# Analysis 2: Effect of farming system on Competition and Cooperation of bacteria vs fungi
    # Use combined dataset
    # Separate datasets into objects for each Farming System and each plot
Org6.107<-subset_samples(fb.RH, Loc_plot_ID=="107") # plots: 107, 109, 110, 204, 205, 206, 308, 310, 311, 412, 413, 415, 
Org6.109<-subset_samples(fb.RH, Loc_plot_ID=="109")
Org6.110<-subset_samples(fb.RH, Loc_plot_ID=="110")
Org6.204<-subset_samples(fb.RH, Loc_plot_ID=="204")
Org6.205<-subset_samples(fb.RH, Loc_plot_ID=="205")
Org6.206<-subset_samples(fb.RH, Loc_plot_ID=="206")
Org6.308<-subset_samples(fb.RH, Loc_plot_ID=="308")
Org6.310<-subset_samples(fb.RH, Loc_plot_ID=="310")
Org6.311<-subset_samples(fb.RH, Loc_plot_ID=="311")
Org6.412<-subset_samples(fb.RH, Loc_plot_ID=="412")
Org6.413<-subset_samples(fb.RH, Loc_plot_ID=="413")
Org6.415<-subset_samples(fb.RH, Loc_plot_ID=="415")
# plots: 113, 114, 115, 209, 210, 211, 301, 302, 303, 409, 410, 411, 
Org_3.113<-subset_samples(fb.RH, Loc_plot_ID=="113")
Org_3.114<-subset_samples(fb.RH, Loc_plot_ID=="114")
Org_3.115<-subset_samples(fb.RH, Loc_plot_ID=="115")
Org_3.209<-subset_samples(fb.RH, Loc_plot_ID=="209")
Org_3.210<-subset_samples(fb.RH, Loc_plot_ID=="210")
Org_3.211<-subset_samples(fb.RH, Loc_plot_ID=="211")
Org_3.301<-subset_samples(fb.RH, Loc_plot_ID=="301")
Org_3.302<-subset_samples(fb.RH, Loc_plot_ID=="302")
Org_3.303<-subset_samples(fb.RH, Loc_plot_ID=="303")
Org_3.409<-subset_samples(fb.RH, Loc_plot_ID=="409")
Org_3.410<-subset_samples(fb.RH, Loc_plot_ID=="410")
Org_3.411<-subset_samples(fb.RH, Loc_plot_ID=="411")

NT.103<-subset_samples(fb.RH, Loc_plot_ID=="103") # plots: 103, 104, 117, 215, 216, 217, 304, 305, 306, 403, 407, 408, 
NT.104<-subset_samples(fb.RH, Loc_plot_ID=="104")
NT.117<-subset_samples(fb.RH, Loc_plot_ID=="117")
NT.215<-subset_samples(fb.RH, Loc_plot_ID=="215")
NT.216<-subset_samples(fb.RH, Loc_plot_ID=="216")
NT.217<-subset_samples(fb.RH, Loc_plot_ID=="217")
NT.304<-subset_samples(fb.RH, Loc_plot_ID=="304")
NT.305<-subset_samples(fb.RH, Loc_plot_ID=="305")
NT.306<-subset_samples(fb.RH, Loc_plot_ID=="306")
NT.403<-subset_samples(fb.RH, Loc_plot_ID=="403")
NT.407<-subset_samples(fb.RH, Loc_plot_ID=="407")
NT.408<-subset_samples(fb.RH, Loc_plot_ID=="408")

CT.101<-subset_samples(fb.RH, Loc_plot_ID=="101") # plots: 101, 102, 116, 212, 213, 214, 313, 314, 315, 404,405, 406,
CT.102<-subset_samples(fb.RH, Loc_plot_ID=="102") 
CT.116<-subset_samples(fb.RH, Loc_plot_ID=="116") 
CT.212<-subset_samples(fb.RH, Loc_plot_ID=="212") 
CT.213<-subset_samples(fb.RH, Loc_plot_ID=="213") 
CT.214<-subset_samples(fb.RH, Loc_plot_ID=="214") 
CT.313<-subset_samples(fb.RH, Loc_plot_ID=="313") 
CT.314<-subset_samples(fb.RH, Loc_plot_ID=="314") 
CT.315<-subset_samples(fb.RH, Loc_plot_ID=="315") 
CT.404<-subset_samples(fb.RH, Loc_plot_ID=="404") 
CT.405<-subset_samples(fb.RH, Loc_plot_ID=="405") 
CT.406<-subset_samples(fb.RH, Loc_plot_ID=="406") 


    # Format data for HMSC
    # Model 1: All edaphic factors as predictors (date/time as random)
    # Apply the model to each dataset
list.CT<-list(CT.101, CT.102, CT.116, CT.212, CT.213, CT.214, CT.313, CT.315, CT.404, CT.405, CT.406)
list.NT<-list(NT.103, NT.104, NT.117, NT.215, NT.216, NT.217, NT.304, NT.305, NT.306, NT.403, NT.407, NT.408)
list.Org3<-list(Org_3.113, Org_3.114, Org_3.115, Org_3.209, Org_3.210, Org_3.211, Org_3.301, Org_3.302, Org_3.303)
list.Org6<-list(Org_6.107, Org_6.109, Org_6.110, Org_6.204, Org_6.205, Org_6.206, Org_6.308, Org_6.310, Org_6.311, Org_6.412, Org_6.413, Org_6.415)
runHMSC<-function(list, method=NULL){
  out<-sapply(list, applyHMSC, simplify=F, USE.NAMES = T, method)
  out
}
applyHMSC<-function(ps, method){
  out<-NULL
  out$probit<-NULL
  out$LogPoi<-NULL
  Ydat<-as.data.frame(as.matrix(otu_table(ps)))
  XData<-as.data.frame(as.matrix(sample_data(ps)), stringsAsFactors = TRUE)
  XDat1<-XData[,c(2,3,6,8:10,15,47,75,81)] # subset data to what I need for this run!!
  rownames(XDat1)<-c(1:nrow(XDat1))
  if(is.null(method)){
    XFormula= ~Glyphosphate_Treatment + genotype + crop + pH + OM...}
  if(method=="AM"){
    XFormula= ~Glyphosphate_Treatment + genotype + crop + pH + OM... + NLFA.AM.Fungi
  }
  XDat1$Sample<-as.factor(c(1:nrow(XDat1)))
  studyDesign = data.frame("Sample"=b.XDat1$Sample,"Sampling_date"=XDat1$Sampling_date, "year"=XDat1$year)
  rL1 <- HmscRandomLevel(units=studyDesign$Sample)
  rL2 <- HmscRandomLevel(units=studyDesign$Sampling_date) # set random level to loc_plot_ID
  rL3 <- HmscRandomLevel(units=studyDesign$year)
  # first do probit model ####
  bf<-Hmsc(Y=as.matrix(Ydat), XData=XDat, XFormula=XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Sampling_date"=rL2, "year"=rL3), distr="probit") 
  out$probit$bf<-sampleMcmc(bf, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)
  
  # examine correlation matrix for probit model
  # get outputs for HMSC analysis
  out$probitGQ2$OmegaCor.probit=computeAssociations(out$probitGQ2$bf)
  out$probitGQ2$preds<-computePredictedValues(out$probitGQ2$bf)
  out$probitGQ2$MF<-evaluateModelFit(out$probitGQ2$bf, predY=out$probitGQ2$preds)
  if(is.null(method)){out$probitGQ2$VP<-computeVariancePartitioning(out$probitGQ2$bf, group=c(1,1,1,1,2,2,3,3,3,3), groupnames=c("Expt", "Edaphic", "Random"))}
  if(method=="AM"){out$probitGQ2$VP<-computeVariancePartitioning(out$probitGQ2$bf, group=c(1,1,1,1,2,2,3,4,4,4,4), groupnames=c("Expt", "Edaphic","AMF", "Random"))}
  out$probitGQ2$postBeta<-getPostEstimate(out$probitGQ2$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal
  
  Y2<-as.matrix(Ydat) # global OTU table
  Y2[Y2==0]<-NA # remove zero counts
  out$LogPoi$bf<-Hmsc(Y=Y2, XData=b.XDat, XFormula=XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="poisson") # use same formulas as previous
  out$LogPoi$bf<-sampleMcmc(out$LogPoiQ2$bf, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)
  out$LogPoi$OmegaCor.lp=computeAssociations(out$LogPoiGQ2$bf)
  out$LogPoi$preds<-computePredictedValues(out$LogPoiGQ2$bf)
  out$LogPoi$MF<-evaluateModelFit(out$LogPoiGQ2$bf, predY=out$LogPoiGQ2$preds)
  if(is.null(method)){out$LogPoi$VP<-computeVariancePartitioning(out$LogPoiGQ2$bf, group=c(1,1,1,1,2,2,3,3,3,3), groupnames=c("Expt", "Edaphic", "Random"))}
  if(method=="AM"){out$LogPoi$VP<-computeVariancePartitioning(out$LogPoiGQ2$bf, group=c(1,1,1,1,2,2,3,4,4,4,4), groupnames=c("Expt", "Edaphic","AMF", "Random"))}
  out$LogPoi$postBeta<-getPostEstimate(out$LogPoiGQ2$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal
  out
}

CTQ2<-runHMSC(list.CT)
NTQ2<-runHMSC(list.NT)
O3Q2<-runHMSC(list.Org3)
O6Q2<-runHMSC(list.Org6)

saveRDS(CTQ2, "Data/Q2CT.RDS")
saveRDS(NTQ2, "Data/Q2NT.RDS")
saveRDS(O3Q2, "Data/Q2O2.RDS")
saveRDS(O6Q2, "Data/Q2O6.RDS")

CTQ3<-runHMSC(list.CT, method="AM")
NTQ3<-runHMSC(list.NT,method="AM")
O3Q3<-runHMSC(list.Org3, method="AM")
O6Q3<-runHMSC(list.Org6,method="AM")

saveRDS(CTQ3, "Data/Q3CT.RDS")
saveRDS(NTQ3, "Data/Q3NT.RDS")
saveRDS(O3Q3, "Data/Q3O2.RDS")
saveRDS(O6Q3, "Data/Q3O6.RDS")

q()
# Analysis 3: Effect of AMF Abundance on bacterial and fungal networks
# format data for HMSC
# Model 1: All edaphic factors - AMF as predictors (date/time as random)
# Model 2: All edaphic factors + AMF as predictors (date/time as random)
# difference between Model 1 and Model 2 is effect of AMF on association matrix So:
# extract associations from Model 1 and Model 2
# subtract Model 2 from Model 1 (Model 2 has no effect of AMF, Model 1 has effect)
# Identify large effects of AMF, if any exist
# check the variation explained by AMF for each of those taxa
# AMF explanatory power:
# interactions with high AMF effect on both taxa are strongly supported
# interactions with only one taxon affected are not supported
# interactions with no affect of AMF are likely not due to AMF themselves
# Descriptive statistics on collections after filtering for support:
# number of +/- edges
# number of taxa in each category
# identity of taxa
# ID of high centrality score taxa in network
# effect of cropping system on AMF and selected highly central taxa

####### End of File #########




# format data for HMSC



# Analysis 4: Effect of Glyphosate on bacterial+fungal networks
# use fungal + bacterial dataset
# format data for HMSC
# Model 1: All edaphic factors - Glyphosate as predictors (date/time as random)
# Model 2: All edaphic factors + Glyphosate as predictors (date/time as random)
# difference between Model 1 and Model 2 is effect of AMF on association matrix So:
# extract associations from Model 1 and Model 2
# subtract Model 2 from Model 1 (Model 2 has no effect of AMF, Model 1 has effect)
# Identify large effects of AMF, if any exist
# check the variation explained by AMF for each of those taxa
# AMF explanatory power:
# interactions with high AMF effect on both taxa are strongly supported
# interactions with only one taxon affected are not supported
# interactions with no affect of AMF are likely not due to AMF themselves
# Descriptive statistics on collections after filtering for support:
# number of +/- edges
# number of taxa in each category
# identity of taxa
# ID of high centrality score taxa in network
# effect of cropping system on AMF and selected highly central taxa
# amf as predictor

# amf as random (extract cooccurence and covariation)

# env. gradients (?) do we need this?

# phylogenetic signal 

# import all PHP data

# compare PLFA bulk to rhizosphere

# total bacteria vs categorical variables

# import Fungal data

# normalize by PLFA to Bac data

# association matrix production for network analysis
