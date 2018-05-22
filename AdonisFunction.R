# Make function ####
# explanatory variables hardcoded in... ####
#Input Arguments:
# ps = phyloseq object with otu table and supplemental data
# n=subset value for significance table; take top and bottom n percentage of slopes from coefficients table
# indics<- list/vector of species of interest to look at coefficients

# make list of pseudomonas function ####
#make list of species from each ps object to look at pseudomonas coefficients
# makes input object for the permanova function for the subset I'm interested in...

makepseudo<- function(ps){
  require(phyloseq)
  udo<-subset_taxa(ps, Genus=="Pseudomonas")
  udo2<-as.matrix(tax_table(udo))
  indics<-rownames(udo2)
  indics
}
# test # function #


# permanova function hardcoded ####


permanova<-function(ps, n, indics) {
  require(phyloseq)
  require(vegan)
  require(stats)
  out=adonis(otu_table(ps) ~ System.loc * Glyphosphate_Treatment * Sampling_date * Soil_Zone, strata = sample_data(ps)$Loc_plot_ID, as(sample_data(ps), "data.frame"))

  tcoeffs<-data.frame(t(out$coefficients))
  tcoeffs$ID<-rownames(tcoeffs)

  syst=rbind(tcoeffs[which(tcoeffs$System.loc1> quantile(tcoeffs$System.loc1, prob=1-n/100)),], tcoeffs[which(tcoeffs$System.loc1< quantile(tcoeffs$System.loc1, prob=n/100)),]) #filter taxa based on steepest n% slopes in system.loc
  
  glyph=rbind(tcoeffs[which(tcoeffs$Glyphosphate_Treatment1> quantile(tcoeffs$Glyphosphate_Treatment1, prob=1-n/100)),], tcoeffs[which(tcoeffs$Glyphosphate_Treatment1< quantile(tcoeffs$Glyphosphate_Treatment1, prob=n/100)),]) #filter taxa based on steepest 5% slopes in Glyphosphate Treatment
  
  date=rbind(tcoeffs[which(tcoeffs$Sampling_date1> quantile(tcoeffs$Sampling_date1, prob=1-n/100)),], tcoeffs[which(tcoeffs$Sampling_date1< quantile(tcoeffs$Sampling_date1, prob=n/100)),]) #filter taxa based on steepest 5% slopes in Sampling Date

a<-NULL
b<-NULL
c<-NULL
d<-NULL
e<-NULL
ind<-as.data.frame(a,b,c,d,e)

for(i in indics){
  first<-tcoeffs[which(tcoeffs$ID==i),]
  ind<-rbind(ind, first)
}

return(list("out"=out, "syst"=syst, "glyph"=glyph, "date"=date, "sppInt"=ind))

}

# test permanova function ####
test.P<-NULL
test2 <- permanova(sv.RRsoy, n=5, indics)

test.P <- permanova(sv.RRsoy, n=5, formula, strata, indics)
test2$out

# parse indicators ####
coef <- coefficients(anova.sv.RRcorn$out)["Glyphosphate_Treatment1:Sampling_date1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
?rev


# try to add names
coef <- coefficients(anova.sv.RRcorn$out)[which(taxa_names(anova.sv.RRcorn)==sv.corn.pseudo)]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
?rev

coef <- coefficients(anova.sv.RRcorn$out)["Glyphosphate_Treatment1:Sampling_date1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
?rev

coef <- coefficients(anova.sv.RRcorn$out)["Glyphosphate_Treatment1:Sampling_date1",]

rownames(anova.sv.RRcorn$sppInt)

?sort
###
barplot(sort(anova.sv.RRcorn$sppInt$Glyphosphate_Treatment1.Soil_Zone1), horiz = T, las=1, main="Species of Interest") # figure out how to do in ggplot2

barplot(anova.sv.RRcorn$sppInt$Glyphosphate_Treatment1.Soil_Zone1, horiz = T, las=1, main="Species of Interest") # figure out how to do in ggplot2
anova.sv.RRcorn


db[order("anova.sv.RRcorn.sppInt.Glyphosphate_Treatment1.Soil_Zone1"),]


db<-data.frame(rownames(anova.sv.RRcorn$sppInt),anova.sv.RRcorn$sppInt$Glyphosphate_Treatment1.Soil_Zone1)
barplot(sort(db$anova.sv.RRcorn.sppInt.Glyphosphate_Treatment1.Soil_Zone1), horiz = T, las=1, main="Species of Interest")

# explanatory variables NOT hardcoded in... ####

permanova2<-function(ps , n, factors, env, strata, indics) {
  formula<-paste0()
  envdat<-paste0()
  out<-adonis(otu_table(ps) ~ formula+envdat, strata, as(sample_data(ps), "data.frame"))
  
  tcoeffs<-data.frame(t(out$coefficients))
  tcoeffs$ID<-rownames(tcoeffs)
  
  for(i in factors){
    i<-rbind(tcoeffs[which(tcoeffs$i> quantile(tcoeffs$i, prob=1-n/100)),], tcoeffs[which(tcoeffs$i< quantile(tcoeffs$i, prob=n/100)),])
    
  }
  #syst<-rbind(tcoeffs[which(tcoeffs$System.loc1> quantile(tcoeffs$System.loc1, prob=1-n/100)),], tcoeffs[which(tcoeffs$System.loc1< quantile(tcoeffs$System.loc1, prob=n/100)),]) #filter taxa based on steepest n% slopes in system.loc
  
  #slopes in Sampling Date
  
  a<-NULL
  b<-NULL
  c<-NULL
  d<-NULL
  e<-NULL
  ind<-as.data.frame(a,b,c,d,e)
  

for(i in indics){
  first<-tcoeffs[which(tcoeffs$ID==i),]
  ind<-rbind(ind, first)
}

return(list(factors, ind))
 
}


# taxa interpretation workspace ####

# scratch ####
indics<-c("TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGTGTGCGCAGGCGGCCGCGCAAGTCGAGTGTGAAAGCCCCGGGCTTAACTTGGGAATTGCGCTCGAAACTACGTGGCTGGAGTGTGGCAGAGGAAGGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAATGGCGAAGGCAGCCTTCTGGGCCAACACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGG", "GACGAACCGTCCAAACGTTATTCGGAATCACTGGGCTTACAGAGTTCGTAGGCGGTCTGGAAAGTTGGGTGTGAAATCCCTCGGCTCAACCGAGGAACTGCGCTTGAAACTACCAGACTCGAGGGAGATAGAGGAAAGCGGAACTGATGGTGGAGCGGTGAAATGCGTTGATATCATCAGGAACACCGGTGGCGAAGGCGGCTTTCTGGGTCTTACCTGACGCTGAGGAACGAAAGCCAGGGGAGCGAACGGG")

indics<-c("TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGTGTGCGCAGGCGGCCGCGCAAGTCGAGTGTGAAAGCCCCGGGCTTAACTTGGGAATTGCGCTCGAAACTACGTGGCTGGAGTGTGGCAGAGGAAGGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAATGGCGAAGGCAGCCTTCTGGGCCAACACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGG")
a<-NULL
b<-NULL
c<-NULL
d<-NULL
e<-NULL
ind<-as.data.frame(a,b,c,d,e)

for(i in indics){
  first<-tcoeffs[which(tcoeffs$ID==i),]
  ind<-rbind(ind, first)
}

head(tcoeffs$ID)

form<-as.formula(System.loc + Glyphosphate_Treatment + Sampling_date + Glyphosphate_Treatment:Sampling_date)
strata = sample_data(sv.RRsoy)$Loc_plot_ID

test2<-permanova2(sv.RRsoy, n=5, formula=System.loc + Glyphosphate_Treatment + Sampling_date + Glyphosphate_Treatment:Sampling_date, strata, indics)


n=5
permanova<-function(ps , n, indics) {
  out=adonis(otu_table(sv.RRsoy) ~ System.loc + Glyphosphate_Treatment + Sampling_date + Glyphosphate_Treatment:Sampling_date, strata = sample_data(sv.RRsoy)$Loc_plot_ID, as(sample_data(sv.RRsoy), "data.frame"))
  
  tcoeffs<-data.frame(t(out$coefficients))
  tcoeffs$ID<-rownames(tcoeffs)
  syst=rbind(tcoeffs[which(tcoeffs$System.loc1> quantile(tcoeffs$System.loc1, prob=1-n/100)),], tcoeffs[which(tcoeffs$System.loc1< quantile(tcoeffs$System.loc1, prob=n/100)),]) #filter taxa based on steepest n% slopes in system.loc
  
  glyph=rbind(tcoeffs[which(tcoeffs$Glyphosphate_Treatment1> quantile(tcoeffs$Glyphosphate_Treatment1, prob=1-n/100)),], tcoeffs[which(tcoeffs$Glyphosphate_Treatment1< quantile(tcoeffs$Glyphosphate_Treatment1, prob=n/100)),]) #filter taxa based on steepest 5% slopes in Glyphosphate Treatment
  
  date=rbind(tcoeffs[which(tcoeffs$Sampling_date1> quantile(tcoeffs$Sampling_date1, prob=1-n/100)),], tcoeffs[which(tcoeffs$Sampling_date1< quantile(tcoeffs$Sampling_date1, prob=n/100)),]) #filter taxa based on steepest 5% slopes in Sampling Date
  
  a<-NULL
  b<-NULL
  c<-NULL
  d<-NULL
  e<-NULL
  ind<-as.data.frame(a,b,c,d,e)
  
for(i in indics){
    first<-tcoeffs[which(tcoeffs$ID==i),]
    ind<-rbind(ind, first)
  }

  return(list("out"=out, "syst"=syst, "glyph"=glyph, "date"=date, "sppInt"=ind))
  
}
first<-tcoeffs[which(tcoeffs$ID==i),]
agh<-list("out"=out, "syst"=syst, "glyph"=glyph, "date"=date)

anova.sv.RRcorn$glyph$Glyphosphate_Treatment1.Soil_Zone1
rownames(anova.sv.RRcorn$glyph)

