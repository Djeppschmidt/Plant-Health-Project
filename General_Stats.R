#import dataset
#scp dietrich.schmidt@scinet-login.bioteam.net:/project/php-bacarc/new_preproc/bacarc_RAWcomDat.rds ~/Desktop/PhD/PHP/bacarc_RAWcomDat.rds

bacarc<-readRDS~("/Desktop/PhD/PHP/bacarc_RAWcomDat.rds")
#general stats

#rarecaction
set.seed(397652)
rare.com<-rarefy_even_depth(ps, sample.size = 20000, rngseed = FALSE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)


library(vegan)
richness<-estimate_richness(rare.com, split=TRUE, measures= "Observed")
S<-specnumber(otu_table(rare.com))
H<-diversity(otu_table(rare.com))
evenness<-H/log(S)


#source: https://github.com/joey711/phyloseq/issues/848
find.top.taxa <- function(x,taxa){
  require(phyloseq)
  top.taxa <- tax_glom(x, taxa)
  otu <- otu_table(t(top.taxa))
  if (taxa_are_rows(otu)){
    otu <- t(otu)
    }
  tax <- tax_table(top.taxa)
  j<-apply(otu,1,which.max)
  k <- j[!duplicated(j)]
  l <- data.frame(tax[k,])
  m <- data.frame(otu[,k])
  colnames(m) = l[,taxa]
  n <- colnames(m)[apply(m,1,which.max)]
  m[,taxa] <- n
  return(m)
}
#use this object to do variance decomposition
maxabund<-find.top.taxa(rare.com, "Genus")


#bray-curtis distance matrix for the general dataset

#adonis for general dataset

#NMS ordination for general dataset

#bray-curtis distance matrix for each subset

#adonis for each subset

#NMS ordination for each subset
