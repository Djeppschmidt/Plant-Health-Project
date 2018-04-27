rare.com<-readRDS("~/Desktop/PhD/PHP/bacarc_RAREcomDat.rds")

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
