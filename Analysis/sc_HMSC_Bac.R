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
    # Model 1: All edaphic factors as predictors (date/time as random)
    # Apply the model to whole dataset
Ydat<-as.data.frame(as.matrix(otu_table(fb.RH)))
XData<-as.data.frame(as.matrix(otu_table(fb.RH)), stringsAsFactors = TRUE)

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
    # Format data for HMSC
    # Model 1: All edaphic factors as predictors (date/time as random)
    # Apply the model to each dataset
    # Make positive and negative association networks
    # collect descriptive stats on networks
    # Run stats on difference among networks in farming systems

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

# Analysis 3: Effect of AMF Abundance on bacterial networks
# use bacterial dataset only
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
RH.XDat1<-b.XDat1[,c(1,4,5,6,8,9,10,11,12,14)] # this needs to be updated!!
rownames(b.XDat1)<-c(1:nrow(b.XDat1))
RH.XDat1$Sample<-as.factor(c(1:nrow(b.XDat1)))
RH.XFormula= ~Glyphosphate_Treatment + genotype # changed since last time
studyDesign = data.frame("Sample"=b.XDat1$Sample,"Loc_plot_ID"=b.XDat1$Loc_plot_ID, "Sampling_date"=b.XDat1$Sampling_date)
rL1 <- HmscRandomLevel(units=studyDesign$Sample)
rL2 <- HmscRandomLevel(units=studyDesign$Loc_plot_ID) # set random level to loc_plot_ID
rL3 <- HmscRandomLevel(units=studyDesign$Sampling_date)
rL4 <- HmscRandomLevel(units=studyDesign$)# second random level is sampling 
# final model should have pH and organic matter as random levels
# first do probit model ####
b.m1<-Hmsc(Y=as.matrix(t(b.Ydat)), XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="probit")
b.m1<-sampleMcmc(b.m.lp, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)

OmegaCor.probit.RH=computeAssociations(b.m1)

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
