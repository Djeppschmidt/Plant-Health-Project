# scinet R script for analysis
# run computationally challenging functions
# output RDS files that are easy to use to build outputs

# load packages for analysis

library(Hmsc)
library(phyloseq)
library(igraph)
library(stringr)
library(Model.Microbiome)
setwd("~/Documents/GitHub/Plant-Health-Project/Analysis")

# import rhizosphere data
b.RH<-readRDS("Data/BacterialRhizospherePhyloseq18112020.rds") # check scinet for new annotation + phylogenetic tree
fun<-readRDS("Data/PHP_Fungi_2020.rds")

# subset fungi to only rhizosphere samples
f.RH<-subset_samples(fun, Soil_Zone=="rhizosphere"&Location=="Beltsville")

# To Do List:
# return to see table mergning for the metadata     [_]
# sort code into proper files                       [_]
# Q1 Probit model                                   [X]
# Q1 Lognormal Model                                [X]
# Q1 (For Q2) Lognormal Model - Farming System      [X]
# Q1 (For Q2) Probit model - Farming System         [X]
# Q1 (For Q3) Lognormal Model + AM fungi            [_]
# Q1 (For Q3) Probit Model + AM fungi               [_]
  
# Q2 Probit for each plot                           [_]
# Q2 lognormal for each plot                        [_]
# Q2 Probit Model + AM Fungi for each plot (Q3)     [_]
# Q2 lognormal Model + AM Fungi for each plot (Q3)  [_]

# Q3 Predict Env. Gradient (Global)                 [_]
# Q3 Compute variation explained by AM (Global)     [_]

# Final Syntax Check                                [_]

# reconcile fungal and bacterial sample abundance
    # create unique sample ID to match bac to fungal samples ...
# check that plates match
p_table<-data.frame("plate.16s"=sample_data(b.RH)$X16s_seq_plate, "plate.ITS"= sample_data(b.RH)$ITS_seq_plate, "M.Plates"=sample_data(b.RH)$ITS_seq_plate)
p_table$M.Plates<-as.character(p_table$M.Plates)
p_table$M.Plates[p_table$M.Plates=="228"]<-"P1"
p_table$M.Plates[p_table$M.Plates=="230"]<-"P2"
p_table$M.Plates[p_table$M.Plates=="233"]<-"P3"
p_table$M.Plates[p_table$M.Plates=="234"]<-"P4"
p_table$M.Plates[p_table$M.Plates=="258"]<-"P5"
p_table$M.Plates[p_table$M.Plates=="262"]<-"P6"
p_table$M.Plates[p_table$M.Plates=="263"]<-"P7"
p_table$M.Plates[p_table$M.Plates=="265"]<-"P8"
#View(p_table) # names of plates are aligned
#str_match(sample_names(b.RH), pattern="S..?")

    # measure fungal to bacterial proportion by NLFA
sample_data(b.RH)$NLFA.Ratio<-sample_data(b.RH)$NLFA.Fungi/(sample_data(b.RH)$NLFA.Gram.Positive+sample_data(b.RH)$NLFA.Gram.Negative)
    # measure fungal to bacterial proportion by PLFA
sample_data(b.RH)$PLFA.Ratio<-sample_data(b.RH)$PLFA.Fungi/(sample_data(b.RH)$PLFA.Gram.Positive+sample_data(b.RH)$PLFA.Gram.Negative)
    # define total counts per sample (as a function of anchor value * relative total PLFA)
hist(sample_sums(subset_taxa(b.RH, Kingdom=="Bacteria"))/(100*(sample_data(b.RH)$NLFA.Gram.Positive+sample_data(b.RH)$NLFA.Gram.Negative)), breaks=30, main="bacterial sequence counts per unit NLFA", xlab="Sequence counts per unit NLFA")

# get total fungal NLFA into fungi ... (including merge step)
sample_names(b.RH)<-paste(as.character(sample_data(b.RH)$ITS_seq_plate), str_match(sample_names(b.RH), pattern="S..?"), sep=".")
fb.RH<-merge_phyloseq(b.RH, f.RH)
fb.RH

f2.RH<-subset_taxa(fb.RH, Kingdom=="k__Fungi")
hist(sample_sums(f2.RH)/(sample_data(f2.RH)$NLFA.Fungi*100), breaks=30, main="fungal sequence counts per unit NLFA", xlab="Sequence counts per unit NLFA")

    # allocate proportion of counts to fungi and bacteria

trfm<-function(ps, val){
  scaled<-data.frame(mapply(`*`, data.frame(t(as.matrix(otu_table(transform_sample_counts(ps, function(x) x/sum(x)))))), val))# sample_data(ps)$val))
  
  rnames<-rownames(data.frame(t(as.matrix(otu_table(ps)))))
  cnames<-sample_names(ps)
  rownames(scaled)<-rnames
  colnames(scaled)<-cnames
  scaled<-round(scaled)
  
  p2<-ps
  otu_table(p2)<- otu_table(scaled, taxa_are_rows=T)
  p2
} # this function takes total bacteria as scaling factor, and assigns counts. Count values of less than 0.5 are rounded to 0; this way taxa that are relatively over sampled are removed from the dataset

# scratch space for transforming bacteria by NLFA
#df.s<-data.frame(mapply(`*`, data.frame(t(as.matrix(otu_table(transform_sample_counts(b.RH, function(x) x/sum(x)))))), sample_data(b.RH)$PLFA.Gram.Positive+sample_data(b.RH)$PLFA.Gram.Negative))

#dim(df.s)
#ns<-rownames(data.frame(t(as.matrix(otu_table(b.RH)))))
#rownames(df.s)<-ns
#b2.RH<-b.RH
#colnames(df.s)<-sample_names(b2.RH) # make sure order is right (it is)
#otu_table(b2.RH)<- otu_table(df.s, taxa_are_rows=T)
#sample_names(b2.RH)

b2.RH<-trfm(b.RH, 100*(sample_data(b.RH)$NLFA.Gram.Negative+sample_data(b.RH)$NLFA.Gram.Positive))
b2.RH<-prune_taxa(taxa_sums(b2.RH) > 0, b2.RH)    

f2.RH<-trfm(f2.RH, 100*sample_data(f2.RH)$NLFA.Fungi)
f2.RH<-prune_taxa(taxa_sums(f2.RH)>0, f2.RH)


    # merge OTU tables with names & make final phyloseq object for analysis

fb.RH<-merge_phyloseq(b.RH, f2.RH)
fb.RH
   
saveRDS(fb.RH, "Data/fb_Rhizosphere_2020.RDS")

# scratch space for species accumulation
sample_data(b2.RH)$NLFA.bac.tot<-sample_data(b2.RH)$NLFA.Gram.Negative+sample_data(b2.RH)$NLFA.Gram.Positive
acum.test<-subset_samples(b2.RH, Location=="Beltsville" & NLFA.bac.tot > 19000 & NLFA.bac.tot<21000)
acum.test<-subset_taxa(acum.test, taxa_sums(acum.test)>0)
acum.test<-tax_glom(acum.test, taxrank="Genus")
acum.test

plot(estimate_richness(acum.test, measures="Observed")$Observed ~ sample_sums(acum.test), col=as.factor(sample_data(acum.test)$System.loc))

out<-specpool(as.data.frame(as.matrix(otu_table(acum.test))))
rarecurve()

# Analysis 1: Role of AMF abundance on microbial networks; Drivers of +/- interaciton networks ; Exploratory analysis of bacterial and fungal competition or symbiosis
    # use combined dataset
fb.RH<-readRDS("Data/fb_Rhizosphere_2020.RDS")
fb.RH<-tax_glom(fb.RH, taxrank="Genus")
    # Model 1: All edaphic factors as predictors (date/time as random)
    # Apply the model to whole dataset
Ydat<-as.data.frame(as.matrix(otu_table(fb.RH)))
XData<-as.data.frame(as.matrix(sample_data(fb.RH)), stringsAsFactors = TRUE)
XDat1<-XData[,c(2,3,6,8:10,15,47,75,81)] # subset data to what I need for this run!!
rownames(b.XDat1)<-c(1:nrow(b.XDat1))
b.XDat1$Sample<-as.factor(c(1:nrow(b.XDat1)))
b.XFormula= ~Glyphosphate_Treatment + genotype + System.loc + crop + pH + OM...
studyDesign = data.frame("Sample"=b.XDat1$Sample,"Sampling_date"=b.XDat1$Sampling_date, "year"=b.XDat1$year)
rL1 <- HmscRandomLevel(units=studyDesign$Sample)
rL2 <- HmscRandomLevel(units=studyDesign$Loc_plot_ID)
rL3 <- HmscRandomLevel(units=studyDesign$Sampling_date) # set random level to loc_plot_ID
rL4 <- HmscRandomLevel(units=studyDesign$year) # second random level is sampling 
# final model should have pH and organic matter as random levels
# first do probit model ####
Q1.bfm1<-Hmsc(Y=as.matrix(t(b.Ydat)), XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Sampling_date"=rL2, "year"=rL3), distr="probit")
probitQ1<-NULL
probitQ1$bf<-sampleMcmc(b.m.lp, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)

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
bf.lp<-Hmsc(Y=G.Y2, XData=XDat1, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="lognormal poisson") # use same formulas as previous
LogPoiQ1$bf.lp<-Hmsc(Y=G.Y2, XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="poisson") # use same formulas as previous

#### Edit to make work for logpoisson model !!!!
LogPoiQ1$OmegaCor.lp=computeAssociations(LogPoiQ1$bf.lp)
LogPoiQ1$preds<-computePredictedValues(LogPoiQ1$bf.lp)
LogPoiQ1$MF<-evaluateModelFit(LogPoiQ1$bf.lp, predY=LogPoiQ1$preds)
LogPoiQ1$VP<-computeVariancePartitioning(LogPoiQ1$bf.lp, group=c(1,1,1,2,1,3,3,4,4,4,4), groupnames=c("Expt", "System", "Edaphic", "Random"))
LogPoiQ1$postBeta<-getPostEstimate(LogPoiQ1$bf.lp, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal

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

#probit1$Gradient.pH<-constructGradient()
#robit1$Gradient.OM<-constructGradient()
#probit1$predpHGradient<-predict()
#probit1$predOMGradient<-predict()

LogPoiQ2<-NULL
bf.lp<-Hmsc(Y=G.Y2, XData=XDat1, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="lognormal poisson") # use same formulas as previous
LogPoiQ2$bf.lp<-Hmsc(Y=G.Y2, XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="poisson") # use same formulas as previous


LogPoiGQ2$OmegaCor.lp=computeAssociations(LogPoiGQ2$bf.lp)
LogPoiGQ2$preds<-computePredictedValues(LogPoiGQ2$bf.lp)
LogPoiGQ2$MF<-evaluateModelFit(LogPoiGQ2$bf.lp, predY=LogPoiGQ2$preds)
LogPoiGQ2$VP<-computeVariancePartitioning(LogPoiGQ2$bf.lp, group=c(1,1,1,1,2,2,3,3,3,3), groupnames=c("Expt", "Edaphic", "Random"))
LogPoiGQ2$postBeta<-getPostEstimate(LogPoiGQ2$bf.lp, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal

saveRDS(LogPoiGQ2, "Data/LogPoiModel2_GDat.RDS")

# construct gradients post if necessary (all data is in the objects created in the main analysis)
# is this necessary? 
LogPoiQ1$Gradient.pH<-constructGradient(bf.lp, focalVariable="pH", non.focalVariables=list("OM..."=list(1), "System.loc"=list(3, "NT")))
LogPoiQ1$Gradient.OM<-constructGradient(bf.lp, focalVariable="OM...", non.focalVariables=list("pH"=list(1), "System.loc"=list(3, "NT")))
LogPoiQ1$predpHGradient<-predict(bf.lp,)
LogPoiQ1$predOMGradient<-predict()

# To do on local machine:
    # extract + associations
    # extract - associations
    # filter for large effects in each
    # filter for bacterial-fungal interactions [examine seperately unique interactions]
    # Construct + / - association networks
    # Identify taxa that are central (either bacteria or fungi)
    # identify hubs
    # Explore each of the hubs and central taxa individually
    # (should the networks include taxa that only strongly associate with others of same type? / compare levels between networks? are there more fungi vs bacteria than within each group?)

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
runHMSC<-function(list, method){
  out<-lapply(list, method)
  out
}
applyHMSC<-function(ps, method){
  out<-NULL
  XDat<-
  Ydat
  
}
# for local computer
    # Make positive and negative association networks
    # collect descriptive stats on networks ( can get from model 1)
    # Run stats on difference among networks in farming systems


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

# format data for HMSC

# probit model (update to reflect the actual analysis?)
RH.Ydat<-as.data.frame(as.matrix(otu_table(RH)))# dim(b.Ydat)
RH.XDat1<-as.data.frame(as.matrix(sample_data(RH)), stringsAsFactors = TRUE)
RH.XDat1<-b.XDat1[,c()] # this needs to be updated!!
rownames(b.XDat1)<-c(1:nrow(b.XDat1))
RH.XDat1$Sample<-as.factor(c(1:nrow(b.XDat1)))
RH.XFormula1= ~Glyphosphate_Treatment + genotype + System.loc + Loc_plot_ID + crop + pH + OM... # changed since last time
RH.XFormula2= ~Glyphosphate_Treatment + genotype +System.loc + Loc_plot_ID + crop + pH + OM...+ NLFA.AM.Fungi
studyDesign = data.frame("Sample"=b.XDat1$Sample,"Loc_plot_ID"=b.XDat1$Loc_plot_ID, "Sampling_date"=b.XDat1$Sampling_date)
rL1 <- HmscRandomLevel(units=studyDesign$Sample)
rL2 <- HmscRandomLevel(units=studyDesign$Loc_plot_ID) # set random level to loc_plot_ID
rL3 <- HmscRandomLevel(units=studyDesign$Sampling_date)
rL4 <- HmscRandomLevel(units=studyDesign$)# second random level is sampling 
# final model should have pH and organic matter as random levels
# first do probit model ####
b.m1<-Hmsc(Y=as.matrix(t(b.Ydat)), XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="probit")
b.m1<-sampleMcmc(b.m.lp, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)

OmegaCor.lp=computeAssociations(bf.lp)
probit3$OmegaCor.probit=computeAssociations(probit1$bfm1)
probit3$preds<-computePredictedValues(probit1$bfmq)
probit3$MF<-evaluateModelFit(probit1$bgm1, predY=probit1$preds)
probit3$VP<-computeVariancePartitioning(probit1$bfm1, group=c(), groupnames=c())
probit3$postBeta<-getPostEstimate(probit1$bfm1, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal
probit3$GradientNT.AM<-constructGradient(bf.m1, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "NT-MD")))
probit3$GradientCT.AM<-constructGradient(bf.lp, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "CT-MD")))
probit3$GradientO3.AM<-constructGradient(bf.lp, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "Org_3")))
probit13$GradientO6.AM<-constructGradient(bf.lp, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "Org_6")))
probit3$predAMGradientNT<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)# rename things!!
probit3$predAMGradientCT<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)
probit3$predAMGradientO3<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)
probit3$predAMGradientO6<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)

Y2<-as.matrix(t(b.Ydat))
Y2[Y2==0]<-NA
# lognormal poisson:
b.m1<-Hmsc(Y=as.matrix(t(b.Ydat)), XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="lognormal poisson")
b.m1<-sampleMcmc(b.m.lp, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)

OmegaCor.lp=computeAssociations(bf.lp)
probit1$OmegaCor.probit=computeAssociations(probit1$bfm1)
probit1$preds<-computePredictedValues(probit1$bfmq)
probit1$MF<-evaluateModelFit(probit1$bgm1, predY=probit1$preds)
probit1$VP<-computeVariancePartitioning(probit1$bfm1, group=c(), groupnames=c())
probit1$postBeta<-getPostEstimate(probit1$bfm1, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal
probit3$GradientNT.AM<-constructGradient(bf.lp, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "NT-MD")))
probit3$GradientCT.AM<-constructGradient(bf.lp, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "CT-MD")))
probit3$GradientO3.AM<-constructGradient(bf.lp, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "Org_3")))
probit13$GradientO6.AM<-constructGradient(bf.lp, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "Org_6")))
probit3$predAMGradientNT<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)# rename things!!
probit3$predAMGradientCT<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)
probit3$predAMGradientO3<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)
probit3$predAMGradientO6<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)

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
