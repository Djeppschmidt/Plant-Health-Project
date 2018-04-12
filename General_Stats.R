#import dataset
#scp dietrich.schmidt@scinet-login.bioteam.net:/project/php-bacarc/new_preproc/bacarc_RAWcomDat.rds ~/Desktop/PhD/PHP/bacarc_RAWcomDat.rds
library(phyloseq)
bacarc<-readRDS("~/Desktop/PhD/PHP/bacarc_RAWcomDat.rds")
#general stats

#rarecaction
set.seed(397652)
rare.com<-rarefy_even_depth(bacarc, sample.size = 20000, rngseed = FALSE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
saveRDS(rare.com, "~/Desktop/PhD/PHP/bacarc_RAREcomDat.rds")
rare.com

install.packages("vegan", dependencies = TRUE)
library(vegan)
richness<-estimate_richness(rare.com, split=TRUE, measures= "Observed")
S<-specnumber(otu_table(rare.com))
H<-diversity(otu_table(rare.com))
evenness<-H/log(S)

divdf<-data.frame(richness, S, H, evenness)
head(divdf)

with(divdf, plot(evenness~richness$Observed))
with(divdf, boxplot(richness$Observed~System.loc, par(las=2), main="Richness"))
with(divdf, boxplot(evenness~System.loc, par(las=2), main="Evenness"))
with(divdf, boxplot(evenness~Glyphosphate_Treatment+System.loc, par(las=2), main="Evenness"))
with(divdf, boxplot(richness$Observed~Glyphosphate_Treatment+System.loc, par(las=2), cex.axis=0.5 ,main="Richness"))
with(divdf, boxplot(evenness~Glyphosphate_Treatment+System.loc, par(las=2), cex.axis=0.5 ,main="Evenness"))



with(divdf, summary(aov(richness$Observed~Glyphosphate_Treatment+System.loc)))

with(divdf, summary(aov(richness$Observed~Loc_plot_ID*System.loc/Sampling_date+Glyphosphate_Treatment)))
#this gives significant Glhyphosphate treatment to Richenss ->
with(divdf, summary(aov(richness$Observed~System.loc/Sampling_date+Glyphosphate_Treatment)))
with(divdf, summary(aov(evenness~System.loc/Sampling_date+Glyphosphate_Treatment)))


#with(divdf, summary(aov(richness$Observed~System.loc/Sampling_date/Glyphosphate_Treatment)))

?boxplot

#regress data?
plot(divdf$S~System.loc)
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

#not sure this is worthwhile... the values in PERMANOVA are pretty robust to order...

maxabund<-find.top.taxa(rare.com, "Genus")

library(nlme)
library(lme4)
install.packages("lme4", dependencies = TRUE)


model<-as.formula(RB41~(1|Location)+(1|Sampling_date)+(1|Loc_plot_ID)+(1|Soil_Zone)+
                  (1|genotype)+(1|crop)+(1|Herbicide_History)+(1|Glyphosphate_Treatment)+(1|System.loc)+(1|year)+
                  (1|Location))
lme(otu~...)
for(i in as.data.frame(maxabund)){
  a<-nlme(otu_table(rare.com)$i~model, data=test)
  a}

attach(sample_data(rare.com))
attach(maxabund)

test<-data.frame(sample_data(rare.com), maxabund)
head(test)

#no interactions

#a<-lme(RB41~Location)
#a<-lme(RB41~(1|model.matrix(Location)), data=test)
       
a<-lmer(RB41~(1|Location)+(1|Sampling_date)+(1|Loc_plot_ID)+(1|Soil_Zone)+
          (1|genotype)+(1|crop)+(1|Herbicide_History)+(1|Glyphosphate_Treatment)+(1|System.loc)+(1|year)+
          (1|Location), data=test)
#add interactions!
B<-lmer(RB41~(1|Location)
        +(1|Sampling_date)
        +(1|Loc_plot_ID)+(1|Loc_plot_ID:Sampling_date)
        +(1|Soil_Zone)+(1|Soil_Zone:Sampling_date)+(1|Soil_Zone:crop)+(1|Soil_Zone:Herbicide_History)+(1|Soil_Zone:Glyphosphate_Treatment)+(1|Soil_Zone:System.loc)
        +(1|genotype)+(1|genotype:Sampling_date)+(1|genotype:Soil_Zone)+(1|genotype:genotype)+(1|genotype:Herbicide_History)+(1|genotype:System.loc)
        +(1|crop)+(1|crop:Sampling_date)+(1|crop:Soil_Zone)+(1|crop:genotype)+(1|crop:Herbicide_History)+(1|crop:System.loc)
        +(1|Herbicide_History)
        +(1|Glyphosphate_Treatment)+(1|Glyphosphate_Treatment:Sampling_date)
        +(1|System.loc)+(1|System.loc:Herbicide_History)+(1|System.loc:Sampling_date)+(1|System.loc:Soil_Zone)+(1|System.loc:crop)+(1|System.loc:genotype)+(1|System.loc:Glyphosphate_Treatment)
        +(1|year)
        +(1|Location), data=test)
B
#a<-lme(RB41~Location, maxabund, random = )
?lme
head(sample_data(rare.com))
#bray-curtis distance matrix for the general dataset
bray<-vegdist(otu_table(rare.com), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE)
#adonis for general dataset
perm<-adonis(bray~System.loc+Sampling_date+Glyphosphate_Treatment+Loc_plot_ID+crop+Soil_Zone, test, permutations=999, method=bray)
#this one? ->
adonis(bray~Sampling_date+Sampling_date*Glyphosphate_Treatment%in%Loc_plot_ID+System.loc, test, permutations=999, method=bray)

adonis(bray~(Sampling_date+Sampling_date*Glyphosphate_Treatment)%in%Loc_plot_ID+System.loc, test, permutations=999, method=bray)


adonis(bray~Loc_plot_ID+System.loc, test, permutations=999, method=bray)

adonis(bray~Loc_plot_ID+System.loc+Sampling_date, test, permutations=999, method=bray)

adonis(bray~System.loc+Sampling_date+Glyphosphate_Treatment+Glyphosphate_Treatment*Sampling_date
             +Loc_plot_ID+crop+Soil_Zone, test, permutations=999, method=bray)

?adonis
#NMS ordination for general dataset
ord1<-metaMDS(bray, k=2, trymax=20, autotransform=TRUE)
Ord1<-ordinate(rare.com, "NMDS", "bray")
#plot
?plot_ordination
library(ggplot2)
head(sample_data(rare.com))
p1<-plot_ordination(rare.com, Ord1, color="Location")+geom_point(size=2)+ggtitle("samples")+theme_classic() + scale_shape_identity()

p1
plot_ordination(rare.com, Ord1, color="System.loc")+geom_point(size=0.005)+theme_classic() + scale_shape_identity()
plot_ordination(rare.com, Ord1, color="crop")+geom_point(size=0.001)+theme_classic() + scale_shape_identity()
plot_ordination(rare.com, Ord1, color="Glyphosphate_Treatment")+geom_point(size=2)+theme_classic() + scale_shape_identity()
plot_ordination(rare.com, Ord1, color="Herbicide_History")+geom_point(size=2)+theme_classic() + scale_shape_identity()


#bray-curtis distance matrix for each subset

#adonis for each subset

#NMS ordination for each subset
