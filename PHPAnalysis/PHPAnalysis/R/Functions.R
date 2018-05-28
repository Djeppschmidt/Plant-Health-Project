

# pick indicators bacteria

makepseudo<- function(ps){
  require(phyloseq)
  udo<-subset_taxa(ps, Genus=="Pseudomonas")
  udo2<-as.matrix(tax_table(udo))
  indics<-rownames(udo2)
  indics
}

# pick indicators fungi

makefusarium<- function(ps){
  require(phyloseq)
  udo<-subset_taxa(ps, Genus=="g__Fusarium")
  udo2<-as.matrix(tax_table(udo))
  indics<-rownames(udo2)
  indics
}

#permanova function

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
  rownames(ind)<-indics # testing...
  ind<-data.frame(ind, factor(indics)) #not sure if this will work
  return(list("out"=out, "syst"=syst, "glyph"=glyph, "date"=date, "sppInt"=ind))

}

# deseq diff abundance test

diffDESeq<-function(ps){
  L1 <-phyloseq_to_deseq2(ps, ~ group + Sampling_date + group:Sampling_date)

  diagdds = DESeq(L1, test="Wald", fitType="parametric")

  res = results(diagdds)

  return(res)
}

# extract p values from deseq test

extractP<-function(ps, pvalue){
  if (identical(rownames(tax_table(ps)), rownames(pvalue))) {
    T1<-ps
    new_tax<-as.data.frame(as.matrix(tax_table(T1)))
    new_tax$DESeq_padj<-pvalue$padj
    tax_table(T1)<-tax_table(as.matrix(new_tax))
    T1
  }
  else {return("rownames out of order")}
}
