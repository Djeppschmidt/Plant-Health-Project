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

# Final Syntax Check                                [X] 
# pilot run functions on local computer             [X] ***
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
probitQ1$VP<-computeVariancePartitioning(probitQ1$bf, group=c(1,1,1,2,1,3,3), groupnames=c("Expt", "System", "Edaphic"))
probitQ1$postBeta<-getPostEstimate(probitQ1$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal

saveRDS(probitQ1, "Data/ProbitModel1Dat.RDS") # plot on local machine

#probit1$Gradient.pH<-constructGradient()
#robit1$Gradient.OM<-constructGradient()
#probit1$predpHGradient<-predict()
#probit1$predOMGradient<-predict()
G.Y2<-as.matrix(t(Ydat)) # global OTU table
G.Y2[G.Y2==0]<-NA # remove zero counts

LogPoiQ1<-NULL
LogPoiQ1$bf<-Hmsc(Y=G.Y2, XData=XDat1, XFormula=XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="poisson") # use same formulas as previous
LogPoiQ1$bf<-sampleMcmc(b.m.lp, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)

LogPoiQ1$OmegaCor.lp=computeAssociations(LogPoiQ1$bf)
LogPoiQ1$preds<-computePredictedValues(LogPoiQ1$bf)
LogPoiQ1$MF<-evaluateModelFit(LogPoiQ1$bf, predY=LogPoiQ1$preds)
LogPoiQ1$VP<-computeVariancePartitioning(LogPoiQ1$bf, group=c(1,1,1,2,1,3,3), groupnames=c("Expt", "System", "Edaphic"))
LogPoiQ1$postBeta<-getPostEstimate(LogPoiQ1$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal

saveRDS(LogPoiQ1, "Data/LogPoiModel1Dat.RDS")

# Q2: effect of farming system: same as above but model removes Farming system 

GQ2.XFormula= ~Glyphosphate_Treatment + genotype + crop + pH + OM...
# first do probit model ####
GQ2.bf<-Hmsc(Y=Ydat, XData=XDat1, XFormula=GQ2.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Sampling_date"=rL2, "year"=rL3), distr="probit")
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

LogPoiGQ2<-NULL
LogPoiGQ2$bf<-Hmsc(Y=G.Y2, XData=XDat1, XFormula=GQ2.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="poisson") # use same formulas as previous
LogPoiGQ2$bf<-sampleMcmc(LogPoiGQ2$bf, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)

LogPoiGQ2$OmegaCor.lp=computeAssociations(LogPoiGQ2$bf.lp)
LogPoiGQ2$preds<-computePredictedValues(LogPoiGQ2$bf.lp)
LogPoiGQ2$MF<-evaluateModelFit(LogPoiGQ2$bf.lp, predY=LogPoiGQ2$preds)
LogPoiGQ2$VP<-computeVariancePartitioning(LogPoiGQ2$bf.lp, group=c(1,1,1,1,2,2,3,3,3,3), groupnames=c("Expt", "Edaphic", "Random"))
LogPoiGQ2$postBeta<-getPostEstimate(LogPoiGQ2$bf.lp, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal

saveRDS(LogPoiGQ2, "Data/LogPoiModel2_GDat.RDS")

# Q3
GQ3.XFormula= ~Glyphosphate_Treatment + genotype + System.loc + crop + pH + OM... + NLFA.AM.Fungi

probitGQ3<-NULL
probitGQ3$bf<-Hmsc(Y=Ydat, XData=XDat1, XFormula=GQ3.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Sampling_date"=rL2, "year"=rL3), distr="probit")
probitGQ3$bf<-sampleMcmc(probitGQ3$bf, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)

# examine correlation matrix for probit model
# get outputs for HMSC analysis
probitGQ3$OmegaCor.probit=computeAssociations(probitGQ3$bf)
probitGQ3$preds<-computePredictedValues(probitGQ3$bf)
probitGQ3$MF<-evaluateModelFit(probitGQ3$bf, predY=probitGQ3$preds)
probitGQ3$VP<-computeVariancePartitioning(probitGQ3$bf, group=c(1,1,1,1,1,2,2,3), groupnames=c("Expt", "Edaphic","AMF"))
probitGQ3$postBeta<-getPostEstimate(probitGQ3$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal

saveRDS(probitGQ3, "Data/ProbitModel3_GDat.RDS") # plot on local machine

LogPoiGQ3<-NULL
LogPoiGQ3$bf<-Hmsc(Y=G.Y2, XData=XDat1, XFormula=GQ3.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="poisson") # use same formulas as previous
LogPoiGQ3$bf<-sampleMcmc(LogPoiGQ3$bf, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100) # check on previous sampleMCMC commands !!!!

#### Edit to make work for logpoisson model !!!!
LogPoiGQ3$OmegaCor.lp=computeAssociations(LogPoiGQ3$bf)
LogPoiGQ3$preds<-computePredictedValues(LogPoiGQ3$bf)
LogPoiGQ3$MF<-evaluateModelFit(LogPoiGQ3$bf, predY=LogPoiGQ3$preds)
LogPoiGQ3$VP<-computeVariancePartitioning(LogPoiGQ3$bf, group=c(1,1,1,1,1,2,2,3), groupnames=c("Expt", "Edaphic", "AMF"))
LogPoiGQ3$postBeta<-getPostEstimate(LogPoiGQ3$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal
saveRDS(LogPoiGQ3, "Data/LogPoiModel3_Gdat.RDS")



# Analysis 2: Effect of farming system on Competition and Cooperation of bacteria vs fungi
    # Use combined dataset
    # Separate datasets into objects for each Farming System and each plot
# automation / Q2-3 probit & poisson
# plots: 107, 109, 110, 204, 205, 206, 308, 310, 311, 412, 413, 415,
Org6.107s<-subset_samples(test, Loc_plot_ID=="107")  # 1 crop: soy
Org6.109c<-subset_samples(test, Loc_plot_ID=="109")  # 1 crop: corn
Org6.110c<-subset_samples(test, Loc_plot_ID=="110"&crop=="corn")  # 2 crops
Org6.110s<-subset_samples(test, Loc_plot_ID=="110"&crop=="soy")
Org6.204s<-subset_samples(test, Loc_plot_ID=="204")  # 1 crop: soy
Org6.205c<-subset_samples(test, Loc_plot_ID=="205")  # 1 crop: corn
Org6.206s<-subset_samples(test, Loc_plot_ID=="206"&crop=="soy")  # 2 crops 
Org6.206c<-subset_samples(test, Loc_plot_ID=="206"&crop=="corn")
Org6.308c<-subset_samples(test, Loc_plot_ID=="308"&crop=="corn")  # 2 crops
Org6.308s<-subset_samples(test, Loc_plot_ID=="308"&crop=="soy")
Org6.310c<-subset_samples(test, Loc_plot_ID=="310")  # 1 crop: corn
Org6.311s<-subset_samples(test, Loc_plot_ID=="311")  # 1 crop: soy
Org6.412c<-subset_samples(test, Loc_plot_ID=="412")  # 1 crop: corn
Org6.413s<-subset_samples(test, Loc_plot_ID=="413")  # 1 crop: soy
Org6.415c<-subset_samples(test, Loc_plot_ID=="415"&crop=="corn")  # 2 crops
Org6.415s<-subset_samples(test, Loc_plot_ID=="415"&crop=="soy")  # 2 crops
list.Org6c<-list(Org6.109c, Org6.110c, Org6.205c, Org6.206c, Org6.308c, Org6.310c, Org6.412c, Org6.415c)
list.Org6s<-list(Org6.107s, Org6.110s, Org6.204s, Org6.206s, Org6.308s, Org6.311s, Org6.413s, Org6.415s)
# plots: 113, 114, 115, 209, 210, 211, 301, 302, 303, 409, 410, 411, 
Org_3.113c<-subset_samples(test, Loc_plot_ID=="113"&crop=="corn") # 2 crops
Org_3.113s<-subset_samples(test, Loc_plot_ID=="113"&crop=="soy") 
Org_3.114s<-subset_samples(test, Loc_plot_ID=="114") # soy
Org_3.115c<-subset_samples(test, Loc_plot_ID=="115") # corn
Org_3.209c<-subset_samples(test, Loc_plot_ID=="209") # corn
Org_3.210s<-subset_samples(test, Loc_plot_ID=="210") # soy
Org_3.211s<-subset_samples(test, Loc_plot_ID=="211"&crop=="soy") # 2 crops
Org_3.211c<-subset_samples(test, Loc_plot_ID=="211"&crop=="corn")
Org_3.301c<-subset_samples(test, Loc_plot_ID=="301") # corn
Org_3.302s<-subset_samples(test, Loc_plot_ID=="302"&crop=="soy") # 2 crops
Org_3.302c<-subset_samples(test, Loc_plot_ID=="302"&crop=="corn")
Org_3.303s<-subset_samples(test, Loc_plot_ID=="303") # soy
Org_3.409s<-subset_samples(test, Loc_plot_ID=="409") # soy
Org_3.410c<-subset_samples(test, Loc_plot_ID=="410") # corn
Org_3.411c<-subset_samples(test, Loc_plot_ID=="411"&crop=="corn") # 2 crops
Org_3.411s<-subset_samples(test, Loc_plot_ID=="411"&crop=="soy") 
list.Org3c<-list(Org_3.113c, Org_3.115c, Org_3.209c, Org_3.211c, Org_3.301c, Org_3.302c, Org_3.410c, Org_3.411c)
list.Org3s<-list(Org_3.113s, Org_3.114s, Org_3.210s, Org_3.211s, Org_3.302s, Org_3.303s, Org_3.409s, Org_3.411s)
# plots: 103, 104, 117, 215, 216, 217, 304, 305, 306, 403, 407, 408, 
NT.103s<-subset_samples(test, Loc_plot_ID=="103"&crop=="soy") # 2 
NT.103c<-subset_samples(test, Loc_plot_ID=="103"&crop=="corn")
NT.104s<-subset_samples(test, Loc_plot_ID=="104") # s
NT.117c<-subset_samples(test, Loc_plot_ID=="117") # c
NT.215c<-subset_samples(test, Loc_plot_ID=="215") # c
NT.216c<-subset_samples(test, Loc_plot_ID=="216"&crop=="corn") # 2
NT.216s<-subset_samples(test, Loc_plot_ID=="216"&crop=="soy")
NT.217s<-subset_samples(test, Loc_plot_ID=="217") # s
NT.304s<-subset_samples(test, Loc_plot_ID=="304") # s
NT.305s<-subset_samples(test, Loc_plot_ID=="305"&crop=="soy") # 2
NT.305c<-subset_samples(test, Loc_plot_ID=="305"&crop=="corn")
NT.306c<-subset_samples(test, Loc_plot_ID=="306") # c
NT.403c<-subset_samples(test, Loc_plot_ID=="403") # c
NT.407s<-subset_samples(test, Loc_plot_ID=="407"&crop=="soy") # 2
NT.407c<-subset_samples(test, Loc_plot_ID=="407"&crop=="corn")
NT.408s<-subset_samples(test, Loc_plot_ID=="408") # s

list.NTc<-list(NT.103c, NT.117c, NT.215c, NT.216c, NT.305c, NT.306c, NT.403c, NT.407c)
list.NTs<-list(NT.103s, NT.104s, NT.216s, NT.217s, NT.304s, NT.305s, NT.407s, NT.408s)
# plots: 101, 102, 116, 212, 213, 214, 313, 314, 315, 404,405, 406,
CT.101s<-subset_samples(test, Loc_plot_ID=="101"&crop=="soy") # 2
CT.101c<-subset_samples(test, Loc_plot_ID=="101"&crop=="corn")
CT.102s<-subset_samples(test, Loc_plot_ID=="102") # S
CT.116c<-subset_samples(test, Loc_plot_ID=="116") # c
CT.212s<-subset_samples(test, Loc_plot_ID=="212") # s
CT.213s<-subset_samples(test, Loc_plot_ID=="213"&crop=="soy") # 2
CT.213c<-subset_samples(test, Loc_plot_ID=="213"&crop=="corn")
CT.214c<-subset_samples(test, Loc_plot_ID=="214") # c
CT.313c<-subset_samples(test, Loc_plot_ID=="313") # c
CT.314s<-subset_samples(test, Loc_plot_ID=="314"&crop=="soy") # 2
CT.314c<-subset_samples(test, Loc_plot_ID=="314"&crop=="corn")
CT.315s<-subset_samples(test, Loc_plot_ID=="315") # s
CT.404c<-subset_samples(test, Loc_plot_ID=="404") # c
CT.405s<-subset_samples(test, Loc_plot_ID=="405") # s
CT.406s<-subset_samples(test, Loc_plot_ID=="406"&crop=="soy") # 2
CT.406c<-subset_samples(test, Loc_plot_ID=="406"&crop=="corn")

# Format data for HMSC
# Model 1: All edaphic factors as predictors (date/time as random)
# Apply the model to each dataset
list.CTc<-list(CT.101c, CT.116c, CT.213c, CT.214c, CT.313c, CT.404c, CT.406c)
list.CTs<-list(CT.101s, CT.102s, CT.212s, CT.213s, CT.315s, CT.405s, CT.406s)
applyHMSC<-function(ps){
  out<-NULL
  
  Ydat<-as.data.frame(as.matrix(otu_table(ps)))
  XData<-as.data.frame(as.matrix(sample_data(ps)), stringsAsFactors = TRUE)
  XDat1<-XData[,c(2,3,6,8:10,15,47,75,81)] # subset data to what I need for this run!!
  rownames(XDat1)<-c(1:nrow(XDat1))
  XDat1$pH<-as.character(XDat1$pH)
  XDat1$OM...<-as.character(XDat1$OM...)
  XDat1$pH[is.na(XDat1$pH)]<-0
  XDat1$OM...[is.na(XDat1$OM...)]<-0
  # continue as usual:
  XFormula1= ~Glyphosphate_Treatment # consider genotype later crop not necessary, metadata not necessary ~Glyphosphate_Treatment + genotype + pH + OM...
  
  XFormula2= ~Glyphosphate_Treatment + NLFA.AM.Fungi # crop not necessary
  
  XDat1$Sample<-as.factor(c(1:nrow(XDat1)))
  studyDesign = data.frame("Sample"=XDat1$Sample,"Sampling_date"=XDat1$Sampling_date)#, "year"=XDat1$year)
  
  #rL3 <- HmscRandomLevel(units=studyDesign$year)
  # first do probit model ####
  #old<-Sys.time()
  run<-function(Ydat, XDat, XFormula, studyDesign, AM){
    out<-NULL
    out$probit<-NULL
    out$LogPoi<-NULL
    rL1 <- HmscRandomLevel(units=studyDesign$Sample)
    rL2 <- HmscRandomLevel(units=studyDesign$Sampling_date) # set random level to loc_plot_ID
    bf<-Hmsc(Y=as.matrix(Ydat), XData=XDat1, XFormula=XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Sampling_date"=rL2), distr="probit")#, "year"=rL3), distr="probit") 
    out$probit$bf<-sampleMcmc(bf, thin=2, samples=1000, transient=500, nChains=1, nParallel=1, verbose=100)
    
    # examine correlation matrix for probit model
    # get outputs for HMSC analysis
    out$probit$OmegaCor.probit=computeAssociations(out$probit$bf)
    print("probit1")
    out$probit$preds<-computePredictedValues(out$probit$bf)
    print("preds1")
    out$probit$MF<-evaluateModelFit(out$probit$bf, predY=out$probit$preds)
    print("MF1")
    print(out$probit$bf$X)
    if(AM==0){out$probit$VP<-computeVariancePartitioning(out$probit$bf, group=c(1,1), groupnames=c("Glyphosate"))} # removed last 3 from vector so that random accounts for removing the year from the model, and remove 1 from expt to account for removing crop from the model
    if(AM==1){out$probit$VP<-computeVariancePartitioning(out$probit$bf, group=c(1,1,2), groupnames=c("Glyphosate", "AMF"))}
    
    out$probit$postBeta<-getPostEstimate(out$probit$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal
    
    Y2<-as.matrix(Ydat) # global OTU table
    Y2[Y2==0]<-NA # remove zero counts
    out$LogPoi$bf<-Hmsc(Y=Y2, XData=XDat1, XFormula=XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Sampling_date"=rL2), distr="poisson") # use same formulas as previous
    out$LogPoi$bf<-sampleMcmc(out$LogPoi$bf, thin=2, samples=1000, transient=500, nChains=1, nParallel=1, verbose=100)
    print("bf Completed") # troubleshooting
    out$LogPoi$OmegaCor.lp=computeAssociations(out$LogPoi$bf)
    print("associations Completed") # troubleshooting
    out$LogPoi$preds<-computePredictedValues(out$LogPoi$bf)
    print("preds Completed") # troubleshooting
    out$LogPoi$MF<-evaluateModelFit(out$LogPoi$bf, predY=out$LogPoi$preds)
    print("model fit Completed") # troubleshooting
    if(AM==0){out$LogPoi$VP<-computeVariancePartitioning(out$LogPoi$bf, group=c(1,1), groupnames=c("Glyphosate"))}
    if(AM==1){out$LogPoi$VP<-computeVariancePartitioning(out$LogPoi$bf, group=c(1,1,2), groupnames=c("Glyphosate", "AMF"))}
    print("VP Completed") # troubleshooting
    out$LogPoi$postBeta<-getPostEstimate(out$LogPoi$bf, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal
    print("Completed")
    out
  }
  
  out$REF<-run(Ydat=Ydat, XDat=XDat1, XFormula=XFormula1, studyDesign=studyDesign, AM=0)
  out$AM<-run(Ydat=Ydat, XDat=XDat1, XFormula=XFormula2, studyDesign=studyDesign, AM=1)
  out$ps<-ps
  out
}
runHMSC<-function(list){
  out<-sapply(list, applyHMSC, simplify=F, USE.NAMES = T)
  out
}
CTc<-runHMSC(list.CTc)
CTs<-runHMSC(list.CTs)
NTc<-runHMSC(list.NTc)
NTs<-runHMSC(list.NTs)
O3c<-runHMSC(list.Org3c)
O3s<-runHMSC(list.Org3s)
O6c<-runHMSC(list.Org6c)
O6s<-runHMSC(list.Org6s)

saveRDS(CTc, "Data/CTc.RDS")
saveRDS(NTc, "Data/NTc.RDS")
saveRDS(O3c, "Data/O3c.RDS")
saveRDS(O6c, "Data/O6c.RDS")

saveRDS(CTs, "Data/CTs.RDS")
saveRDS(NTs, "Data/NTs.RDS")
saveRDS(O3s, "Data/O3s.RDS")
saveRDS(O6s, "Data/O6s.RDS")
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
# (Already finished as part of other functions!!)
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
