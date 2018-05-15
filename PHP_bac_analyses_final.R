# session configuration ####

setwd("~/Desktop/PhD/PHP/")

library(ggplot2)
library(data.table)
library(phyloseq)
library(gplots)
library(VennDiagram)

require(devtools)
install_version("vegan", version = "2.4-6", repos = "http://cran.us.r-project.org")
library(vegan)
library(DESeq2)

# creation of phyloseq object ####

# post-processing and error correction of metadata.
# correct version is imported below
##
#phpmeta2<-read.table("Master_Label_Table_PHPv1.csv", header = T, na.strings = "", sep = ",")
#phpmeta2[2]<-NULL
#phpmeta3<-phpmeta2[,-1]
#rownames(phpmeta3)<-phpmeta2[,1]
#saveRDS(phpmeta2, file = "fungi_metadata.rds")

#phpmeta$System.loc<-as.character(phpmeta$System.loc)
#phpmeta$System.loc[phpmeta$System.loc=="CT-FF"]<-"CT" #other transformations done as well
#saveRDS(phpmeta, file = "phpmeta_corrected.rds")

#phpmeta<-readRDS("phpmeta_corrected.rds")
#phpmeta$year<-as.factor(phpmeta$year)
#seqtab.new<-readRDS("fungi_seqs2.rds")
#phpfuntax.new<-readRDS(file = "fungi_taxa2.rds")
#phpfuntax.new <- cbind( phpfuntax.new, "seq" = rownames(phpfuntax.new), "fasta_name" = paste0("phptaxa", seq(1:nrow(phpfuntax.new))))



#ps.new <- phyloseq(otu_table(seqtab.new, taxa_are_rows=FALSE), 
#                   sample_data(phpmeta), 
#                   tax_table(phpfuntax.new))
#ntaxa(ps.new)

#taxa_names(ps.new)[1:5]
sam<-read.csv("~/Documents/GitHub/Plant-Health-Project/Master_SampleData_Table_PHPv3_2.csv", colClasses = "factor")
rownames(sam)<-sam$X.SampleID


ps.new<-readRDS("~/Desktop/PhD/PHP/bacarc_RAWcomDat.rds")
sample_data(ps.new)<-sam

#phpbactax.new<-tax_table(ps.new)

#phpbactax.new <- cbind( phpbactax.new, "seq" = rownames(phpbactax.new), "fasta_name" = paste0("phpbactaxa", seq(1:nrow(phpbactax.new))))

#if (identical(taxa_names(ps.new), as.character(phpbactax.new[,"seq"]))) {
#  taxa_names(ps.new) <- as.character(phpbactax.new[,"fasta_name"])
#}

#taxa_names(ps.new)[1:5]

#ps.new <- subset_taxa(ps.new, Kingdom == "k__Fungi")
#ntaxa(ps.new)
unique(as.data.frame(as.matrix(tax_table(ps.new)))$Kingdom)
table(as.data.frame(as.matrix(tax_table(ps.new)))$Kingdom)

ps.new<-subset_taxa(ps.new, Kingdom!="Eukaryota")
###
###
# Subseting by location ####


# NMDS for all samples, all sites ####
# At this point in the workflow there is no filtering of low abundance taxa
# You might want to do that later, but you will need to build a filtered
# phyloseq object 
# remove Urbana sites first
ps.newREL <- subset_samples(ps.new, Location != "Urbana")

#ps.newREL  = transform_sample_counts(ps.newREL, function(x) x / sum(x) )

set.seed(197)
all.ord.nmds.bray <- ordinate(ps.newREL, k = 4, method="NMDS", distance="bray", maxit = 30000, sratmax = 0.99999999, sfgrmin = 1e-10)#, sratmax = 0.9999999)#, k = 2, sfgrmin = 1e-8, sratmax = 0.999999, maxit = 30000)
all.plot<-plot_ordination(ps.newREL, all.ord.nmds.bray, type="samples", color="Location")#, shape="Location") 
all.plot + geom_point(size = 2.5)#+ geom_polygon(geom_point(size=5) + ggtitle("samples"))
ggsave("all_sites_nmds_reltrans_101602_taxa.pdf", plot = last_plot(), device = pdf, height = 6, width = 8)

sam<-read.csv("~/Documents/GitHub/Plant-Health-Project/Master_SampleData_Table_PHPv3_2.csv", colClasses = "factor")

?read.csv
cols<-colnames(sample_data(ps.new))
sam<-sample_data(ps.new)
cols<-colnames(sam)
sam[cols]<-lapply(sam[cols], factor)
head(sample_data(ps.new))


# FSP specific ####
# Should be run from unfiltered data
# Make FSP-RR-Soy only phyloseq object from unfiltered original 
fsp.RRsoy <- subset_samples(ps.new, Location == "Beltsville" & genotype == "RR" & crop == "soy") #  & Sampling_date == "pre" & crop == "soy"  & Glyphosphate_Treatment == "no_spray")
fsp.RRsoy <- prune_taxa(taxa_sums(fsp.RRsoy) > 10, fsp.RRsoy)
ntaxa(fsp.RRsoy)
nsamples(fsp.RRsoy)


# DESeq routines
# make DESeq object
fspDS <- phyloseq_to_deseq2(fsp.RRsoy, ~ Glyphosphate_Treatment)

# Add column for system+gly
fspDS$group <- factor(paste(fspDS$System.loc, fspDS$Glyphosphate_Treatment))
levels(fspDS$group)<-sub(" ", "", levels(fspDS$group))

# Set the model to something meaningful

design(fspDS) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

# Perform variance stabilization
# Start the clock
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


ptm <- proc.time()
geoMeans = apply(counts(fspDS), 1, gm_mean)
fspDS = estimateSizeFactors(fspDS, geoMeans=geoMeans)
fspDS = estimateDispersions(fspDS)
fspVST = getVarianceStabilizedData(fspDS)
# Stop the clock
proc.time() - ptm

# use VST counts
# backup original FSP-Soy-RR phyloseq object
fsp.RRsoy0 <- fsp.RRsoy
fspVST0 <- fspVST
fspVST[fspVST < 0.0] <- 0.0

otu_table(fsp.RRsoy) <- otu_table(fspVST, taxa_are_rows = TRUE)

fsp.ord.pca.bray <- ordinate(fsp.RRsoy, method="PCoA", distance="bray")

#fsp.RRsoy.eigen <- eigen(fsp.ord.pca.bray)

fsp.plot<-plot_ordination(fsp.RRsoy, fsp.ord.pca.bray, type="samples", color="System.loc", shape="Glyphosphate_Treatment") 
fsp.plot + ggtitle("FSP-Soy") +
  geom_point(size = 3)

fsp.RRsoy.axisvals <- fsp.ord.pca.bray$vectors

bac.fsp.soy.qiime <- merge(as.data.frame(sample_data(fsp.RRsoy)), fsp.RRsoy.axisvals, by = "row.names")
write.table(bac.fsp.soy.qiime, sep = "\t", file = "fsp.soy.qiime.txt", row.names = F, quote = F)

# FSP Soy PERMANOVA
metadata <- as(sample_data(fsp.RRsoy), "data.frame")
fsp.dist <- phyloseq::distance(fsp.RRsoy, "bray")
adonis(fsp.dist ~ System.loc + 
         year + 
         Glyphosphate_Treatment + 
         Sampling_date + 
         Soil_Zone + 
         Glyphosphate_Treatment + 
         Glyphosphate_Treatment:Soil_Zone + 
         Glyphosphate_Treatment:Sampling_date +
         Glyphosphate_Treatment:Sampling_date:System.loc:year, 
       as(sample_data(fsp.RRsoy), "data.frame"))

# FSP-RR-Corn ####
# phyloseq object from unfiltered original
fsp.RRcorn <- subset_samples(ps.new, Location == "Beltsville" & genotype == "RR" & crop == "corn") #  & Sampling_date == "pre" & crop == "corn"  & Glyphosphate_Treatment == "no_spray")
fsp.RRcorn <- prune_taxa(taxa_sums(fsp.RRcorn) > 10, fsp.RRcorn)
ntaxa(fsp.RRcorn)
nsamples(fsp.RRcorn)

# DESeq routines
# inherits phyloseq object from FSP-only 

fspDC <- phyloseq_to_deseq2(fsp.RRcorn, ~ Glyphosphate_Treatment)

fspDC$group <- factor(paste(fspDC$System.loc, fspDC$Glyphosphate_Treatment))

levels(fspDC$group)<-sub(" ", "", levels(fspDC$group))

design(fspDC) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

# Start the clock
ptm <- proc.time()
geoMeans = apply(counts(fspDC), 1, gm_mean)
fspDC = estimateSizeFactors(fspDC, geoMeans=geoMeans)
fspDC = estimateDispersions(fspDC)
fspVSTC = getVarianceStabilizedData(fspDC)
# Stop the clock
proc.time() - ptm

fsp.RRcorn0 <- fsp.RRcorn
fspVSTC0 <- fspVSTC
fspVSTC[fspVSTC < 0.0] <- 0.0

otu_table(fsp.RRcorn) <- otu_table(fspVSTC, taxa_are_rows = TRUE)
fsp.ord.pca.bray <- ordinate(fsp.RRcorn, method="PCoA", distance="bray")

fsp.plot<-plot_ordination(fsp.RRcorn, fsp.ord.pca.bray, type="samples", color="System.loc", shape="Glyphosphate_Treatment") 
fsp.plot + ggtitle("FSP-corn") +
  geom_point(size = 3)

fsp.RRcorn.axisvals <- fsp.ord.pca.bray$vectors

bac.fsp.corn.qiime <- merge(as.data.frame(sample_data(fsp.RRcorn)), fsp.RRcorn.axisvals, by = "row.names")
write.table(bac.fsp.corn.qiime, sep = "\t", file = "fsp.corn.qiime.txt", row.names = F, quote = F)

# Stoneville specific ####
# Should be run from unfiltered data
# Make sv-RR-Soy only phyloseq object from unfiltered original 
sv.RRsoy <- subset_samples(ps.new, Location == "Stoneville" & genotype == "RR" & crop == "soy") #  & Sampling_date == "pre" & crop == "soy"  & Glyphosphate_Treatment == "no_spray")
sv.RRsoy <- prune_taxa(taxa_sums(sv.RRsoy) > 10, sv.RRsoy)
ntaxa(sv.RRsoy)
nsamples(sv.RRsoy)

# DESeq routines
# make DESeq object
svDS <- phyloseq_to_deseq2(sv.RRsoy, ~ Glyphosphate_Treatment)

# Add column for system+gly
svDS$group <- factor(paste(svDS$System.loc, svDS$Glyphosphate_Treatment))
levels(svDS$group)<-sub(" ", "", levels(svDS$group))

# Set the model to something meaningful
design(svDS) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

# Perform variance stabilization
# Start the clock
ptm <- proc.time()
geoMeans = apply(counts(svDS), 1, gm_mean)
svDS = estimateSizeFactors(svDS, geoMeans=geoMeans)
svDS = estimateDispersions(svDS)
svVSTC = getVarianceStabilizedData(svDS)
# Stop the clock
proc.time() - ptm

# use VSTC counts
# backup original sv-Soy-RR phyloseq object
sv.RRsoy0 <- sv.RRsoy
svVSTC0 <- svVSTC
svVSTC[svVSTC < 0.0] <- 0.0

otu_table(sv.RRsoy) <- otu_table(svVSTC, taxa_are_rows = TRUE)
sv.ord.pca.bray <- ordinate(sv.RRsoy, method="PCoA", distance="bray")
sv.RRsoy.eigen <- eigen(sv.ord.pca.bray)

sv.plot<-plot_ordination(sv.RRsoy, sv.ord.pca.bray, type="samples", color="year", shape="Glyphosphate_Treatment") 
sv.plot + ggtitle("sv-Soy Relative transform NMDS (k=2) 16092 taxa") +
  geom_point(size = 3)

sv.RRsoy.axisvals <- sv.ord.pca.bray$vectors

bac.sv.soy.qiime <- merge(as.data.frame(sample_data(sv.RRsoy)), sv.RRsoy.axisvals, by = "row.names")
write.table(bac.sv.soy.qiime, sep = "\t", file = "sv.soy.qiime.txt", row.names = F, quote = F)

# sv Soy PERMANOVA
metadata <- as(sample_data(sv.RRsoy), "data.frame")
sv.dist <- phyloseq::distance(sv.RRsoy, "bray")
adonis(sv.dist ~ System.loc + 
         year + 
         Glyphosphate_Treatment + 
         Sampling_date + 
         Soil_Zone + 
         Glyphosphate_Treatment + 
         Glyphosphate_Treatment:Soil_Zone + 
         Glyphosphate_Treatment:Sampling_date +
         Glyphosphate_Treatment:Sampling_date:System.loc:year, 
       as(sample_data(sv.RRsoy), "data.frame"))

# sv-RR-Corn ####
# phyloseq object from unfiltered original
sv.RRcorn <- subset_samples(ps.new, Location == "Stoneville" & genotype == "RR" & crop == "corn") #  & Sampling_date == "pre" & crop == "corn"  & Glyphosphate_Treatment == "no_spray")
sv.RRcorn <- prune_taxa(taxa_sums(sv.RRcorn) > 10, sv.RRcorn)
ntaxa(sv.RRcorn)
nsamples(sv.RRcorn)

# DESeq routines
# inherits phyloseq object from sv-only 

svDC <- phyloseq_to_deseq2(sv.RRcorn, ~ Glyphosphate_Treatment)

svDC$group <- factor(paste(svDC$System.loc, svDC$Glyphosphate_Treatment))

levels(svDC$group)<-sub(" ", "", levels(svDC$group))

design(svDC) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

# Start the clock
ptm <- proc.time()
geoMeans = apply(counts(svDC), 1, gm_mean)
svDC = estimateSizeFactors(svDC, geoMeans=geoMeans)
svDC = estimateDispersions(svDC)
svVSTC = getVarianceStabilizedData(svDC)
# Stop the clock
proc.time() - ptm

sv.RRcorn0 <- sv.RRcorn
svVSTC0 <- svVSTC
svVSTC[svVSTC < 0.0] <- 0.0

otu_table(sv.RRcorn) <- otu_table(svVSTC, taxa_are_rows = TRUE)
sv.ord.pca.bray <- ordinate(sv.RRcorn, method="PCoA", distance="bray")

sv.plot<-plot_ordination(sv.RRcorn, sv.ord.pca.bray, type="samples", color="year", shape="Glyphosphate_Treatment") 
sv.plot + ggtitle("sv-corn Relative transform NMDS (k=2) 16092 taxa") +
  geom_point(size = 3)

sv.RRcorn.axisvals <- sv.ord.pca.bray$vectors

bac.sv.corn.qiime <- merge(as.data.frame(sample_data(sv.RRcorn)), sv.RRcorn.axisvals, by = "row.names")
write.table(bac.sv.corn.qiime, sep = "\t", file = "sv.corn.qiime.txt", row.names = F, quote = F)

