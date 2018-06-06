# Session configuration ####

setwd("~/Documents/PHP/fungi")
install.packages("doMC")
install.packages("foreach")
install.packages("data.table")
install.packages("pheatmap")
library(ggplot2)
library(data.table)
library(phyloseq)
library(gplots)
library(VennDiagram)
library(vegan)
library(DESeq2)
library(doMC)
library(foreach)
library(pheatmap)
options(warnings=-1)

# creation of phyloseq object ####

# post-processing and error correction of metadata ####
# correct version is imported below
# this can be disregarded, but included for historical considerations

#phpmeta2<-read.table("Master_Label_Table_PHPv1.csv", header = T, na.strings = "", sep = ",")
#phpmeta2[2]<-NULL
#phpmeta3<-phpmeta2[,-1]
#rownames(phpmeta3)<-phpmeta2[,1]
#saveRDS(phpmeta2, file = "fungi_metadata.rds")

#phpmeta$System.loc<-as.character(phpmeta$System.loc)
#phpmeta$System.loc[phpmeta$System.loc=="CT-FF"]<-"CT" #other transformations done as well
#saveRDS(phpmeta, file = "phpmeta_corrected.rds")

# Import output from dada2 and metadata ####
# create phyloseq object and change taxa names to something shorter

setwd("/Users/maullabmacbookpro/Documents/GitHub/Plant-Health-Project/Glyphosate/outputRaw")


#not normalized first round.
ps.new <- readRDS("/Users/maullabmacbookpro/Desktop/PhD/PHP/bacarc_RAWcomDat.rds")
ntaxa(ps.new)
sample_data(ps.new)
taxa_names(ps.new)[1:5]

# input taxa key!

seqs<-taxa_names(ps.new)
name<-paste0("Seq", seq(ntaxa(ps.new)))
nameKey<-data.frame(name, seqs)
head(nameKey$name)


tax<-as.data.frame(as.matrix(tax_table(ps.new)))
tax$seq<-name
tax$seqs<-seqs
tax_table(ps.new)<-tax_table(as.matrix(tax))
taxa_names(ps.new)<-name
taxa_names(ps.new)[1:5]

ntaxa(ps.new)
ps.new <- subset_taxa(ps.new, Kingdom != "Eukaryota")
ps.new <- subset_samples(ps.new, Location != "Urbana")
ntaxa(ps.new)

phpmeta<-as.data.frame(as.matrix(sample_data(ps.new)))
# Fix Stoneville System.loc names
phpmeta$System.loc <- as.character(phpmeta$System.loc)
phpmeta$Herbicide_History <- as.character(phpmeta$Herbicide_History)
phpmeta$System.loc[phpmeta$Location == "Stoneville"] <- paste(phpmeta$System.loc[phpmeta$Location == "Stoneville"], phpmeta$Herbicide_History[phpmeta$Location == "Stoneville"], sep = "_")
phpmeta$System.loc <- as.factor(phpmeta$System.loc)
phpmeta$Herbicide_History <- as.factor(phpmeta$Herbicide_History)
phpmeta$year<-as.factor(phpmeta$year)
# Make alternate sampleID and variables
phpmeta$group <- factor(paste(phpmeta$System.loc, phpmeta$Glyphosphate_Treatment, sep = "_"))
phpmeta$q_id <- factor(paste(phpmeta$System.loc, phpmeta$Glyphosphate_Treatment, phpmeta$Soil_Zone, phpmeta$year, sep = "_"))


#phpmeta$System.loc[phpmeta$Location == "Stoneville"] <- paste(phpmeta$System.loc[phpmeta$Location == "Stoneville"], phpmeta$Herbicide_History[phpmeta$Location == "Stoneville"], sep = "_")

#taxa_names(ps.new)[1:5]


# Run standard workflow in parallel ####
# create subsets for combinations of sample data values
# this version is for Location, Soil_Zone, crop
# Genotype is hardcoded for this version. See line BLAH

unique.meta<-as.data.frame(unique(phpmeta[c("Location", "crop")]))
#unique.meta <- unique.meta[(unique.meta$Location == "Stoneville"),]

registerDoMC(cores = 8)

ptm <- proc.time()
foreach(a=1:nrow(unique.meta), .packages = c("phyloseq")) %dopar% {
  # get variable values
  loc <- as.character(unique.meta["Location"][a,])
  crop <- as.character(unique.meta["crop"][a,])

  # make a name to recycle
  thing <- paste(loc, crop, sep = "_")

  # make folder tree to save output
  out.dir <- paste(getwd(),loc,thing, sep = "/")
  dir.create(out.dir, recursive = T)

  # determine the samples to work on
  #loc.list <- as.character(get_variable(ps.new, "Location")) == loc
  loc.row <- row.names(sample_data(ps.new)[sample_data(ps.new)$Location %in% loc, ])
  crop.row <- row.names(sample_data(ps.new)[sample_data(ps.new)$crop %in% crop, ])
  genotype.row <- row.names(phpmeta[phpmeta$genotype %in% "RR", ])
  common.row <- Reduce(intersect, list(loc.row, crop.row, genotype.row))

  # prune phyloseq object to desired samples
  ps.0 <- prune_samples(common.row, ps.new) #what is common.row?
  ps.0 <- prune_taxa(taxa_sums(ps.0) > 10, ps.0)
  saveRDS(ps.0, file.path(out.dir, paste(thing,".rds", sep = "")))


  # perform ordinations
  #psVST[psVST < 0.0] <- 0.0
  #otu_table(ps.0) <- otu_table(psVST, taxa_are_rows = TRUE)

  bray <- ordinate(ps.0, method="PCoA", distance="bray")
  ps.plot<-plot_ordination(ps.0, bray, type="samples", color="System.loc") +
    ggtitle(paste(thing, "samples:", nsamples(ps.0), "taxa:", ntaxa(ps.0), sep = " ")) +
    geom_point(size = 3)+theme_classic()
  ggsave(filename = file.path(out.dir, paste(thing,"pdf", sep = ".")), plot = ps.plot, width = 20, height = 14, units = "cm", device = "pdf", dpi = 300)

  # save vectors from PCoA for use in qiime2 longitudinal test
  axisvals <- bray$vectors
  qiime <- merge(as.data.frame(sample_data(ps.0)), axisvals, by = "row.names")
  write.table(qiime, sep = "\t", file = file.path(out.dir, paste(thing,"txt", sep = ".")), row.names = F, quote = F)

  # perform PERMANOVA
  out <- adonis(otu_table(ps.0) ~ System.loc * year * Soil_Zone * Glyphosphate_Treatment * Sampling_date, strata = sample_data(ps.0)$Loc_plot_ID, as(sample_data(ps.0), "data.frame"))
#out <- adonis(t(otu_table(ps.0)) ~ System.loc * year * Soil_Zone * Glyphosphate_Treatment * Sampling_date, strata = sample_data(ps.0)$Loc_plot_ID, as(sample_data(ps.0), "data.frame"))
  saveRDS(out, file.path(out.dir, paste(thing,"perm", "rds", sep = ".")))
  write.table(as.matrix(cbind(rownames(out$aov.tab), out$aov.tab)),
              file = file.path(out.dir, paste(thing,"perm", "txt", sep = ".")),
              row.names = F,
              col.names = c("Variable", colnames(out$aov.tab)),
              quote = F, sep = "\t")
}
proc.time() - ptm

# Unix and Perl one-liner to combine all PERMANOVA output.
# Execute from top level directory for each site. Edit as necessary.
# tail -n +1 Stoneville_*/*.perm.txt | perl -p -e 's/^.*\/(\w+)\.perm\.txt.*$/$1/g' > all_permanova.txt

setwd("/Users/maullabmacbookpro/Documents/GitHub/Plant-Health-Project/Glyphosate/outputRare")


#not normalized first round.
ps.new <- readRDS("/Users/maullabmacbookpro/Desktop/PhD/PHP/bacarc_RAREcomDat.rds")
ntaxa(ps.new)
sample_data(ps.new)
taxa_names(ps.new)[1:5]

# input taxa key!

seqs<-taxa_names(ps.new)
name<-paste0("Seq", seq(ntaxa(ps.new)))
nameKey<-data.frame(name, seqs)
head(nameKey$name)


tax<-as.data.frame(as.matrix(tax_table(ps.new)))
tax$seq<-name
tax$seqs<-seqs
tax_table(ps.new)<-tax_table(as.matrix(tax))
taxa_names(ps.new)<-name
taxa_names(ps.new)[1:5]

ntaxa(ps.new)
ps.new <- subset_taxa(ps.new, Kingdom != "Eukaryota")
ps.new <- subset_samples(ps.new, Location != "Urbana")
ntaxa(ps.new)

phpmeta<-as.data.frame(as.matrix(sample_data(ps.new)))
# Fix Stoneville System.loc names
phpmeta$System.loc <- as.character(phpmeta$System.loc)
phpmeta$Herbicide_History <- as.character(phpmeta$Herbicide_History)
phpmeta$System.loc[phpmeta$Location == "Stoneville"] <- paste(phpmeta$System.loc[phpmeta$Location == "Stoneville"], phpmeta$Herbicide_History[phpmeta$Location == "Stoneville"], sep = "_")
phpmeta$System.loc <- as.factor(phpmeta$System.loc)
phpmeta$Herbicide_History <- as.factor(phpmeta$Herbicide_History)
phpmeta$year<-as.factor(phpmeta$year)
# Make alternate sampleID and variables
phpmeta$group <- factor(paste(phpmeta$System.loc, phpmeta$Glyphosphate_Treatment, sep = "_"))
phpmeta$q_id <- factor(paste(phpmeta$System.loc, phpmeta$Glyphosphate_Treatment, phpmeta$Soil_Zone, phpmeta$year, sep = "_"))


#phpmeta$System.loc[phpmeta$Location == "Stoneville"] <- paste(phpmeta$System.loc[phpmeta$Location == "Stoneville"], phpmeta$Herbicide_History[phpmeta$Location == "Stoneville"], sep = "_")

#taxa_names(ps.new)[1:5]


# Run standard workflow in parallel ####
# create subsets for combinations of sample data values
# this version is for Location, Soil_Zone, crop
# Genotype is hardcoded for this version. See line BLAH

unique.meta<-as.data.frame(unique(phpmeta[c("Location", "crop")]))
#unique.meta <- unique.meta[(unique.meta$Location == "Stoneville"),]

registerDoMC(cores = 8)

ptm <- proc.time()
foreach(a=1:nrow(unique.meta), .packages = c("phyloseq")) %dopar% {
  # get variable values
  loc <- as.character(unique.meta["Location"][a,])
  crop <- as.character(unique.meta["crop"][a,])
  
  # make a name to recycle
  thing <- paste(loc, crop, sep = "_")
  
  # make folder tree to save output
  out.dir <- paste(getwd(),loc,thing, sep = "/")
  dir.create(out.dir, recursive = T)
  
  # determine the samples to work on
  #loc.list <- as.character(get_variable(ps.new, "Location")) == loc
  loc.row <- row.names(sample_data(ps.new)[sample_data(ps.new)$Location %in% loc, ])
  crop.row <- row.names(sample_data(ps.new)[sample_data(ps.new)$crop %in% crop, ])
  genotype.row <- row.names(phpmeta[phpmeta$genotype %in% "RR", ])
  common.row <- Reduce(intersect, list(loc.row, crop.row, genotype.row))
  
  # prune phyloseq object to desired samples
  ps.0 <- prune_samples(common.row, ps.new) #what is common.row?
  ps.0 <- prune_taxa(taxa_sums(ps.0) > 10, ps.0)
  saveRDS(ps.0, file.path(out.dir, paste(thing,".rds", sep = "")))
  
  
  # perform ordinations
  #psVST[psVST < 0.0] <- 0.0
  #otu_table(ps.0) <- otu_table(psVST, taxa_are_rows = TRUE)
  
  bray <- ordinate(ps.0, method="PCoA", distance="bray")
  ps.plot<-plot_ordination(ps.0, bray, type="samples", color="System.loc") +
    ggtitle(paste(thing, "samples:", nsamples(ps.0), "taxa:", ntaxa(ps.0), sep = " ")) +
    geom_point(size = 3)+theme_classic()
  ggsave(filename = file.path(out.dir, paste(thing,"pdf", sep = ".")), plot = ps.plot, width = 20, height = 14, units = "cm", device = "pdf", dpi = 300)
  
  # save vectors from PCoA for use in qiime2 longitudinal test
  axisvals <- bray$vectors
  qiime <- merge(as.data.frame(sample_data(ps.0)), axisvals, by = "row.names")
  write.table(qiime, sep = "\t", file = file.path(out.dir, paste(thing,"txt", sep = ".")), row.names = F, quote = F)
  
  # perform PERMANOVA
  out <- adonis(otu_table(ps.0) ~ System.loc * year * Soil_Zone * Glyphosphate_Treatment * Sampling_date, strata = sample_data(ps.0)$Loc_plot_ID, as(sample_data(ps.0), "data.frame"))
  #out <- adonis(t(otu_table(ps.0)) ~ System.loc * year * Soil_Zone * Glyphosphate_Treatment * Sampling_date, strata = sample_data(ps.0)$Loc_plot_ID, as(sample_data(ps.0), "data.frame"))
  saveRDS(out, file.path(out.dir, paste(thing,"perm", "rds", sep = ".")))
  write.table(as.matrix(cbind(rownames(out$aov.tab), out$aov.tab)),
              file = file.path(out.dir, paste(thing,"perm", "txt", sep = ".")),
              row.names = F,
              col.names = c("Variable", colnames(out$aov.tab)),
              quote = F, sep = "\t")
}
proc.time() - ptm
