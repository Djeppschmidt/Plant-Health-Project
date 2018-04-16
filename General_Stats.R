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

with(divdf, plot(evenness~richness$Observed), col=as.factor(sample_data(rare.com)$Quantity..picograms))

sample_data(rare.com)$Quantity..picograms
with(divdf, boxplot(richness$Observed~System.loc, par(las=2), main="Richness"))
with(divdf, boxplot(evenness~System.loc, par(las=2), main="Evenness"))
with(divdf, boxplot(evenness~Glyphosphate_Treatment+System.loc, par(las=2), main="Evenness"))
with(divdf, boxplot(richness$Observed~Glyphosphate_Treatment+System.loc, par(las=2), cex.axis=0.5 ,main="Richness"))
with(divdf, boxplot(evenness~Glyphosphate_Treatment+System.loc, par(las=2), cex.axis=0.5 ,main="Evenness"))


# aov expects a balanced design, which this is certainly not!!
#will need to rerun as lme in nlme!!

with(divdf, summary(aov(richness$Observed~Glyphosphate_Treatment+System.loc)))

with(divdf, summary(aov(richness$Observed~Loc_plot_ID*System.loc/Sampling_date+Glyphosphate_Treatment)))
#this gives significant Glhyphosphate treatment to Richenss ->
with(divdf, summary(aov(richness$Observed~System.loc/Sampling_date+Glyphosphate_Treatment)))
with(divdf, summary(aov(evenness~System.loc/Sampling_date+Glyphosphate_Treatment)))

#control for year
with(divdf, summary(aov(richness$Observed~(System.loc/Sampling_date+Glyphosphate_Treatment)%in%year)))
with(divdf, summary(aov(evenness~(System.loc/Sampling_date+Glyphosphate_Treatment)%in%year)))
?aov

#with(divdf, summary(aov(richness$Observed~System.loc/Sampling_date/Glyphosphate_Treatment)))
#lme tests
?lme
#model -> ~ fixed + effects + (1|random) + (|effects)

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

adonis(bray~(Sampling_date+Sampling_date*Glyphosphate_Treatment)%in%(Loc_plot_ID+System.loc), test, permutations=999, method=bray)


adonis(bray~Loc_plot_ID+System.loc, test, permutations=999, method=bray)

adonis(bray~Loc_plot_ID+System.loc+Sampling_date, test, permutations=999, method=bray)

adonis(bray~System.loc+Sampling_date+Glyphosphate_Treatment+Glyphosphate_Treatment*Sampling_date
             +Loc_plot_ID+crop+Soil_Zone, test, permutations=999, method=bray)

adonis(bray~System.loc+Glyphosphate_Treatment
       +Loc_plot_ID+crop+Soil_Zone, test, strata=Sampling_date, permutations=999, method=bray)

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

?adonis
#make subsets
#by place:
?subset_samples
head(sample_data(rare.com))
urbana<-subset_samples(rare.com, Location=="Urbana")
head(sample_data(urbana)$Location)
table(sample_data(urbana)$crop)

Beltsville<-subset_samples(rare.com, Location=="Beltsville")
table(sample_data(Beltsville)$crop)#check to make sure formating is good...
Belt.C<-subset_samples(Beltsville, crop=="corn")
Belt.S<-subset_samples(Beltsville, crop=="soy")
Stoneville<-subset_samples(rare.com, Location=="Stoneville")

Stone.S<-subset_samples(Stoneville, crop=="soy")
Stone.C<-subset_samples(Stoneville, crop=="corn")

table(sample_data(Stone.S)$Location) #check subset
table(sample_data(Stone.S)$crop) #check subset

#calculate richness for each set
#richness<-estimate_richness(rare.com, split=TRUE, measures= "Observed")
#S<-specnumber(otu_table(rare.com))
#H<-diversity(otu_table(rare.com))
#evenness<-H/log(S)

Urb.C.rich<-estimate_richness(urbana, split=TRUE, measures= "Observed")
Belt.C.rich<-estimate_richness(Belt.C, split=TRUE, measures= "Observed")
Belt.S.rich<-estimate_richness(Belt.S, split=TRUE, measures= "Observed")
Stone.C.rich<-estimate_richness(Stone.C, split=TRUE, measures= "Observed")
Stone.S.rich<-estimate_richness(Stone.S, split=TRUE, measures= "Observed")
#calculate evenness for each set

Urb.Sval<-specnumber(otu_table(urbana))
Urb.Hval<-diversity(otu_table(urbana))
Urb.C.Even<-Urb.Hval/log(Urb.Sval)

Belt.C.Sval<-specnumber(otu_table(Belt.C))
Belt.C.Hval<-diversity(otu_table(Belt.C))
Belt.C.Even<-Belt.C.Hval/log(Belt.C.Sval)

Belt.S.Sval<-specnumber(otu_table(Belt.S))
Belt.S.Hval<-diversity(otu_table(Belt.S))
Belt.S.Even<-Belt.S.Hval/log(Belt.S.Sval)

Stone.C.Sval<-specnumber(otu_table(Stone.C))
Stone.C.Hval<-diversity(otu_table(Stone.C))
Stone.C.Even<-Stone.C.Hval/log(Stone.C.Sval)

Stone.S.Sval<-specnumber(otu_table(Stone.S))
Stone.S.Hval<-diversity(otu_table(Stone.S))
Stone.S.Even<-Stone.S.Hval/log(Stone.S.Sval)

#add richness and evenness to sample data
sample_data(urbana)$richness<-Urb.C.rich
sample_data(urbana)$evenness<-Urb.C.Even
sample_data(urbana)$richness
sample_data(Belt.C)$richness<-Belt.C.rich
sample_data(Belt.C)$evenness<-Belt.C.Even

sample_data(Belt.S)$richness<-Belt.S.rich
sample_data(Belt.S)$evenness<-Belt.S.Even

sample_data(Stone.C)$richness<-Stone.C.rich
sample_data(Stone.C)$evenness<-Stone.C.Even

sample_data(Stone.S)$richness<-Stone.S.rich
sample_data(Stone.S)$evenness<-Stone.S.Even
head(sample_data(Stone.S))

#remove  undetermined from the plot
urbana_test<-prune_samples(sample_data(urbana)$Glyphosphate_Treatment!="undetermined", urbana)
table(sample_data(urbana_test)$Glyphosphate_Treatment) 

Belt.C<-prune_samples(sample_data(Belt.C)$Glyphosphate_Treatment!="undetermined", Belt.C)
Belt.S<-prune_samples(sample_data(Belt.S)$Glyphosphate_Treatment!="undetermined", Belt.S)

Stone.C<-prune_samples(sample_data(Stone.C)$Glyphosphate_Treatment!="undetermined", Stone.C)
Stone.S<-prune_samples(sample_data(Stone.S)$Glyphosphate_Treatment!="undetermined", Stone.S)

table(sample_data(Stone.C)$Glyphosphate_Treatment)


?prune_samples
#richness plots and stats

with(sample_data(urbana_test), boxplot(richness$Observed~System.loc+Glyphosphate_Treatment, par(las=2), cex.axis=0.5, main="Urbana Corn Richness"))
 

with(sample_data(Belt.C), boxplot(richness$Observed~System.loc+Glyphosphate_Treatment, par(las=2), cex.axis=0.5, main="Beltsville Corn Richness"))


with(sample_data(Belt.S), boxplot(richness$Observed~System.loc+Glyphosphate_Treatment, par(las=2), cex.axis=0.5, main="Beltsville Soy Richness"))


with(sample_data(Stone.C), boxplot(richness$Observed~System.loc+Glyphosphate_Treatment, par(las=2), cex.axis=0.5, main="Stoneville Corn Richness"))

with(sample_data(Stone.S), boxplot(richness$Observed~System.loc+Glyphosphate_Treatment, par(las=2), cex.axis=0.5, main="Stoneville Soy Richness"))

####
with(sample_data(urbana_test), boxplot(evenness~System.loc+Glyphosphate_Treatment, par(las=2), cex.axis=0.5, main="Urbana Corn Evenness"))


with(sample_data(Belt.C), boxplot(evenness~System.loc+Glyphosphate_Treatment, par(las=2), cex.axis=0.5, main="Beltsville Corn Evenness"))


with(sample_data(Belt.S), boxplot(evenness~System.loc+Glyphosphate_Treatment, par(las=2), cex.axis=0.5, main="Beltsville Soy Evenness"))


with(sample_data(Stone.C), boxplot(evenness~System.loc+Glyphosphate_Treatment, par(las=2), cex.axis=0.5, main="Stoneville Corn Evenness"))

with(sample_data(Stone.S), boxplot(evenness~System.loc+Glyphosphate_Treatment, par(las=2), cex.axis=0.5, main="Stoneville Soy Evenness"))

#bray-curtis distance matrix for each subset

#adonis for each subset

#NMS ordination for each subset
