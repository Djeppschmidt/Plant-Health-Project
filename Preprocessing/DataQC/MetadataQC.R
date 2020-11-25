# reconcile metadata
library(dplyr)
PLFA<-as.data.frame(read.csv("/Users/dietrich/Desktop/PhD/PHP/Hmsc_reproducible/PLFA.csv")) # from excel sheets provided by Jude Maul
Chem<-as.data.frame(read.csv("/Users/dietrich/Desktop/PhD/PHP/key.csv")) # sheet labeled key
CFU<-as.data.frame(read.csv("/Users/dietrich/Desktop/PhD/PHP/CFU.csv"))
meta<-as.data.frame(as.matrix(sample_data(bac))) # from workflow output on scinet (prepared by Ryan?)

# need to build out unique analyses for each because of dataset structure?
# plot level chem data
# PLFA = all with bulk and rhizosphere
# nlfa = no bulk

# interesting comparisons: 
# rhizosphere across samples (nlfa + plfa in plot ID for chem data)
# rhizosphere plfa vs bulk

# paper on achoring and fungal prediction:
# compare plfa anchoring to QPCR (correlation)
# Use nlfa amf to predict bacterial abundance
      # probit model to identify potential symbionts (drill into details?)
# use plfa or nlfa to merge bacterial with fungal datasets
# hurdle probit and hurdle models for interactions
      # what are the appropriate taxa to include?

# subset to BARC
meta<-meta[meta$System.loc=="Org_6" | meta$System.loc=="Org_3" | meta$System.loc=="CT-MD" | meta$System.loc=="NT-MD",]

# data reconciliation tasks: 
# fill chem data to plot level [X]
# fill PLFA data to sample level [X]
# PLFA subset to rhizosphere [X]
# fill NLFA data to rhizosphere sample level [X]
# make phyloseq with rhizosphere [X]
# make phyloseq with all plfa [X]


# prepare phyloseq objects:

# all:
Belt<-subset_samples(bac, Site=="Beltsville")
Belt # 768 samples
# rhizosphere:
RH<-subset_samples(Belt, Soil_Zone=="rhizosphere")
RH # 384 samples
# subset PLFA to same data as meta
View(PLFA)

# make data match formats between data tables
PLFA$TimeOfyear<-as.character(PLFA$TimeOfyear)
PLFA$TimeOfyear[PLFA$TimeOfyear=="Pre"]<-"pre"
PLFA$TimeOfyear[PLFA$TimeOfyear=="Post"]<-"post"
PLFA$genotype<-as.character(PLFA$genotype)
PLFA$genotype[PLFA$genotype=="RR2"]<-"RR"
PLFA$genotype[PLFA$genotype=="AF2"]<-"Non_RR"
PLFA$genotype[PLFA$genotype=="allen"]<-"RR"
PLFA$genotype[PLFA$genotype=="usg56"]<-"Non_RR"
PLFA$crop<-as.character(PLFA$crop)
PLFA$crop[PLFA$crop=="Corn"]<-"corn"
PLFA$crop[PLFA$crop=="Soy"]<-"soy"
PLFA$crop[PLFA$crop=="Soybean"]<-"soy"
PLFA$Glyphosphate.Treatment<-as.character(PLFA$Glyphosphate.Treatment)
PLFA$Glyphosphate.Treatment[PLFA$Glyphosphate.Treatment=="+"]<-"spray"
PLFA$Glyphosphate.Treatment[PLFA$Glyphosphate.Treatment=="-"]<-"no_spray"
PLFA$Bl.Rh<-as.character(PLFA$Bl.Rh)
PLFA$Bl.Rh[PLFA$Bl.Rh=="Bulk"]<-"bulk"
PLFA$Bl.Rh[PLFA$Bl.Rh=="Rhizo"]<-"rhizosphere"
PLFA$system.ID<-as.character(PLFA$system.ID)
PLFA$system.ID[PLFA$system.ID=="CT"]<-"CT-MD"
PLFA$system.ID[PLFA$system.ID=="NT"]<-"NT-MD"
PLFA$system.ID[PLFA$system.ID=="Org 6"]<-"Org_6"
PLFA$system.ID[PLFA$system.ID=="Org 3"]<-"Org_3"
PLFA$Bl.Rh[PLFA$Bl.Rh==""]<-NA

# extract PLFA markers
NLFA<-NLFA[NLFA$MarkerType=="nlfa",]
NLFA2013<-NLFA[NLFA$year==2013,]
NLFA2014<-NLFA[NLFA$year==2014,]

PLFA<-PLFA[PLFA$MarkerType=="plfa",]
PLFA2013<-PLFA[PLFA$year==2013,]
PLFA2014<-PLFA[PLFA$year==2014,]
PLFA$FSPplot<-as.factor(PLFA$FSPplot)

meta2013<-meta[meta$year==2013,]
meta2014<-meta[meta$year==2014,]
# Merge PLFA
PLFA2013$Merge<-paste0(PLFA2013$FSPplot, PLFA2013$Glyphosphate.Treatment, PLFA2013$TimeOfyear, PLFA2013$system.ID, PLFA2013$genotype, PLFA2013$crop, sep="")
meta2013$Merge<-paste0(meta2013$Loc_plot_ID, meta2013$Glyphosphate_Treatment, meta2013$Sampling_date, meta2013$System.loc, meta2013$genotype, meta2013$crop, sep="")

PLFA2014$Merge<-paste0(PLFA2014$FSPplot, PLFA2014$Glyphosphate.Treatment, PLFA2014$TimeOfyear, PLFA2014$Bl.Rh, PLFA2014$system.ID, PLFA2014$genotype, PLFA2014$crop, sep="")
meta2014$Merge<-paste0(meta2014$Loc_plot_ID, meta2014$Glyphosphate_Treatment, meta2014$Sampling_date, meta2014$Soil_Zone, meta2014$System.loc, meta2014$genotype, meta2014$crop, sep="")

m.PLFA2013<-merge(meta2013, PLFA2013, by="Merge", all.x=TRUE)
m.PLFA2014<-merge(meta2014, PLFA2014, by="Merge", all.x=TRUE)

PLFA2<-rbind(m.PLFA2013,m.PLFA2014)
colnames(PLFA2)
PLFA2<-PLFA2[,c(2,5:15,20,39,42,52:60)]
names(PLFA2)[17]<-"PLFA.General.FAME"
names(PLFA2)[18]<-"PLFA.AM.Fungi"
names(PLFA2)[19]<-"PLFA.Gram.Negative"
names(PLFA2)[20]<-"PLFA.Fungi"
names(PLFA2)[21]<-"PLFA.Gram.Positive"
names(PLFA2)[22]<-"PLFA.Anaerobe"
names(PLFA2)[23]<-"PLFA.Actinomycetes"
names(PLFA2)[24]<-"PLFA.Protozoa"
colnames(PLFA2)
View(PLFA2) # PLFA2 is a reconciled dataset!
RH.PLFA2<-PLFA2[PLFA2$Soil_Zone=="rhizosphere",]

# Merge NLFA
NLFA2013$Merge<-paste0(NLFA2013$FSPplot, NLFA2013$Glyphosphate.Treatment, NLFA2013$TimeOfyear, NLFA2013$system.ID, NLFA2013$genotype, NLFA2013$crop, sep="")

NLFA2014$Merge<-paste0(NLFA2014$FSPplot, NLFA2014$Glyphosphate.Treatment, NLFA2014$TimeOfyear, NLFA2014$Bl.Rh, NLFA2014$system.ID, NLFA2014$genotype, NLFA2014$crop, sep="")

m.NLFA2013<-merge(meta2013, NLFA2013, by="Merge", all.x=TRUE)
m.NLFA2014<-merge(meta2014, NLFA2014, by="Merge", all.x=TRUE)

NLFA2<-rbind(m.NLFA2013,m.NLFA2014)

# select columns and change rownames
colnames(NLFA2)
NLFA2<-NLFA2[,c(2,5:15,20,39,42,52:60)]
names(NLFA2)[17]<-"NLFA.General.FAME"
names(NLFA2)[18]<-"NLFA.AM.Fungi"
names(NLFA2)[19]<-"NLFA.Gram.Negative"
names(NLFA2)[20]<-"NLFA.Fungi"
names(NLFA2)[21]<-"NLFA.Gram.Positive"
names(NLFA2)[22]<-"NLFA.Anaerobe"
names(NLFA2)[23]<-"NLFA.Actinomycetes"
names(NLFA2)[24]<-"NLFA.Protozoa"
colnames(NLFA2)
NLFA2<-NLFA2[NLFA2$Soil_Zone=="rhizosphere",] # reconciled!!


View(NLFA2)
# Chem dataset reconciliation

View(Chem)
Chem<-Chem[-c(96:101),]
Chem$Glyphosphate_Treatment<-Chem$Microplot.treatment
Chem$Glyphosphate_Treatment<-as.character(Chem$Glyphosphate_Treatment)
Chem$Glyphosphate_Treatment[Chem$Glyphosphate_Treatment=="allen"]<-"no_spray"
Chem$Glyphosphate_Treatment[Chem$Glyphosphate_Treatment=="allen+gly"]<-"spray"
Chem$Glyphosphate_Treatment[Chem$Glyphosphate_Treatment=="usg56"]<-"no_spray"
Chem$Glyphosphate_Treatment[Chem$Glyphosphate_Treatment=="RR2"]<-"no_spray"
Chem$Glyphosphate_Treatment[Chem$Glyphosphate_Treatment=="RR2+gly"]<-"spray"
Chem$Glyphosphate_Treatment[Chem$Glyphosphate_Treatment=="AF2"]<-"no_spray"
Chem$Microplot.treatment<-as.character(Chem$Microplot.treatment)
Chem$Microplot.treatment[Chem$Microplot.treatment=="allen"]<-"RR"
Chem$Microplot.treatment[Chem$Microplot.treatment=="allen+gly"]<-"RR"
Chem$Microplot.treatment[Chem$Microplot.treatment=="usg56"]<-"Non_RR"
Chem$Microplot.treatment[Chem$Microplot.treatment=="RR2"]<-"RR"
Chem$Microplot.treatment[Chem$Microplot.treatment=="RR2+gly"]<-"RR"
Chem$Microplot.treatment[Chem$Microplot.treatment=="AF2"]<-"Non_RR"

Chem$crop<-as.character(Chem$crop)
Chem$crop[Chem$crop=="soybean"]<-"soy"
Chem$system<-as.character(Chem$system)
Chem$system[Chem$system=="ct"]<-"CT-MD"
Chem$system[Chem$system=="nt"]<-"NT-MD"
Chem$system[Chem$system=="Org6"]<-"Org_6"
Chem$system[Chem$system=="Org3"]<-"Org_3"
Chem$Plot.number<-as.character(Chem$Plot.number)

Chem$MergeC<-paste0(Chem$Plot.number, Chem$Glyphosphate_Treatment, Chem$system, Chem$Microplot.treatment, Chem$crop,sep="")
meta$MergeC<-paste0(meta$Loc_plot_ID, meta$Glyphosphate_Treatment, meta$System.loc, meta$genotype, meta$crop, sep="")

setdiff(Chem$MergeC, meta$MergeC) # should be = 0

Chem2<-merge(meta, Chem, by="MergeC", all.x=TRUE)
RH.Chem2<-Chem2[Chem2$Soil_Zone=="rhizosphere",]


# merge into phyloseq
RH.meta<-left_join(NLFA2, RH.PLFA2)
RH.meta<-left_join(RH.meta, RH.Chem2)

B.meta<-left_join(PLFA2,Chem2) # 
rownames(B.meta)<-B.meta$X.SampleID

rownames(RH.meta)<-RH.meta$X.SampleID
RH<-merge_phyloseq(RH, sample_data(RH.meta))

Belt<-merge_phyloseq(Belt, sample_data(B.meta))

saveRDS(RH, "~/Documents/GitHub/Plant-Health-Project/Analysis/Data/BacterialRhizospherePhyloseq18112020.rds")
saveRDS(Belt, "~/Documents/GitHub/Plant-Health-Project/Analysis/Data/Bacteria18112020.rds")
# ignore cfu for now

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

# pilot analysis:

library(phyloseq)
library(Hmsc)
library(corrplot)
setwd("~/Hmsc_reproducible")
Proteobacteria.Org6.corn<-readRDS("Proteobacteria_Org6_corn.rds") # load data, in phyloseq object

# Run Hmsc on bacteria

b.Ydat<-as.data.frame(as.matrix(otu_table(Proteobacteria.Org6.corn)))#, stringsAsFactors=TRUE) # extract taxon count table
b.Ydat<-as.matrix(t(b.Ydat))
rownames(b.Ydat)<-c(1:nrow(b.Ydat)) # I did this in case the rownames were causing the error; didn't change anything.
b.XDat1<-as.data.frame(as.matrix(sample_data(Proteobacteria.Org6.corn)), stringsAsFactors = TRUE) # extract sample data
b.XDat1<-b.XDat1[,c(1,4,5,6,8,9,10,11,12,14)] # subset to sample data of interest
rownames(b.XDat1)<-c(1:nrow(b.XDat1)) 
b.XDat1$Sample<-as.factor(c(1:nrow(b.XDat1)))
b.XFormula= ~Glyphosphate_Treatment # +Soil_Zone+genotype+Sampling_date # treatments of interest
studyDesign = data.frame("Sample"=b.XDat1$Sample,"Loc_plot_ID"=b.XDat1$Loc_plot_ID)
#rownames(studyDesign)<-b.XDat1$X.SampleID
rL1 <- HmscRandomLevel(units=studyDesign$Sample)
rL2 <- HmscRandomLevel(units=studyDesign$Loc_plot_ID) # set random level to loc_plot_ID
#rL3 <- HmscRandomLevel(units=studyDesign$Sampling_date) # set random level to micro-plot ; same output if using levels()
#rL2 <- HmscRandomLevel(units=studyDesign$Sampling_date) # second random level is sampling date (before vs after glyphosate application)


b.m<-Hmsc(Y=b.Ydat, XData=b.XDat1, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1,"Loc_plot_ID"=rL2)) # Full Hmsc object with all random levels
b.m2<-sampleMcmc(b.m, thin=2, samples=1000, transient=500, nChains=1, nParallel=1, verbose=100) # runs for about 3 hours on one 4 ghz core with all factors in xformula; 10 minutes when only using one factor.
preds=computePredictedValues(b.m2, expected = TRUE) 
OmegaCor=computeAssociations(b.m2)

# produces error: 
# Error in PiNew[, r] <- sapply(dfPiNew[, r], function(s) which(rowNames ==  : 
# incorrect number of subscripts on matrix

# or 

# Error in Eta[[r]][as.character(dfPiNew[, r]), ] : 
#  no 'dimnames' attribute for array

# depends on whether I use two random factors or 1

# Run model with no random levels
b.m3<-Hmsc(Y=as.matrix(t(b.Ydat)), XData=b.XDat, XFormula=b.XFormula) # No random levels
b.m3<-sampleMcmc(b.m3, thin=2, samples=1000, transient=500, nChains=1, nParallel=1, verbose=100) # runs in a couple minutes on one 4 ghz core

preds3=computePredictedValues(b.m3, expected = TRUE) # works as far as I can tell
