

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

  year=rbind(tcoeffs[which(tcoeffs$year1> quantile(tcoeffs$year1, prob=1-n/100)),], tcoeffs[which(tcoeffs$year1< quantile(tcoeffs$year1, prob=n/100)),])

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
  return(list("out"=out,"year"=year, "syst"=syst, "glyph"=glyph, "date"=date, "sppInt"=ind))

}


# modified function for sv samples (with no second system)
permanova2<-function(ps, n, indics) {
  require(phyloseq)
  require(vegan)
  require(stats)
  out=adonis(otu_table(ps) ~ year * Glyphosphate_Treatment * Sampling_date * Soil_Zone, strata = sample_data(ps)$Loc_plot_ID, as(sample_data(ps), "data.frame"))

  tcoeffs<-data.frame(t(out$coefficients))
  tcoeffs$ID<-rownames(tcoeffs)

  year=rbind(tcoeffs[which(tcoeffs$year1> quantile(tcoeffs$year1, prob=1-n/100)),], tcoeffs[which(tcoeffs$year1< quantile(tcoeffs$year1, prob=n/100)),]) #filter taxa based on steepest n% slopes in system.loc

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
  return(list("out"=out, "year"=year, "glyph"=glyph, "date"=date, "sppInt"=ind))

}

# deseq diff abundance test

diffDESeqF<-function(ps){
  L1 <-phyloseq_to_deseq2(ps, ~ Glyphosphate_Treatment:Sampling_date + Glyphosphate_Treatment + Sampling_date + Loc_plot_ID)

  diagdds = DESeq(L1, test="Wald", fitType="parametric")

  res = results(diagdds)

  return(res)
}

diffDESeqB<-function(ps){
  L1 <-phyloseq_to_deseq2(ps, ~ Glyphosphate_Treatment + Sampling_date + Glyphosphate_Treatment:Sampling_date)
  geoMeans = apply(counts(L1), 1, gm_mean)
  L1 = estimateSizeFactors(L1, geoMeans = geoMeans)
  #L1 = DESeq(L1, fitType="local")
  L1 = DESeq(L1, test="Wald", fitType="parametric")

  res = results(L1)

  return(res)
}


# extract p values from deseq test and permanova output
?identical

extractP<-function(ps, DESeqP, PERMV){
  P1<-as.data.frame(as.matrix(t(PERMV$out$coefficients)))
  if (identical(rownames(tax_table(ps)), rownames(DESeqP)) && identical(rownames(tax_table(ps)), rownames(P1))) {
    T1<-ps
    new_tax<-as.data.frame(as.matrix(tax_table(T1)))
    new_tax$DESeq_padj<-DESeqP$padj
    new_tax$PERM_coeff<-P1$`Glyphosphate_Treatment1:Sampling_date1`
    tax_table(T1)<-tax_table(as.matrix(new_tax))
    T1
  }
  else {return("rownames out of order")}
}

#plot differential abundance plot

plotCoeffF<-function(ps, title, xlab, ylab){
  require(ggplot2)
  L1<-subset_taxa(ps, Genus=="g__Fusarium")
  pL1<-as.data.frame(as.matrix(tax_table(L1)))
  pL1$DESeq_padj<-as.numeric(as.character(pL1$DESeq_padj))
  pL1$PERM_coeff<-as.numeric(as.character(pL1$PERM_coeff))
  pL1$DESeq_padj[is.na(pL1$DESeq_padj)]<-1
  pL1$color<-cut(pL1$DESeq_padj, breaks=c(-Inf, 0.05, Inf), labels=c("Significant", "Not Significant"))
  p=ggplot(pL1, aes(x=reorder(fasta_name, -PERM_coeff), PERM_coeff, fill=color))
  p= p + geom_bar(stat="identity")
  p= p + theme_classic()
  p= p + coord_flip()
  p= p + theme(axis.text.x=element_text(hjust=1, size=10), axis.text.y=element_text(hjust=1, size=4))
  p= p + labs(title = title, x=xlab, y=ylab)
  return(p)
}

plotCoeffB<-function(ps, title, xlab, ylab){
  require(ggplot2)
  L1<-subset_taxa(ps, Genus=="Pseudomonas")
  pL1<-as.data.frame(as.matrix(tax_table(L1)))
  pL1$DESeq_padj<-as.numeric(as.character(pL1$DESeq_padj))
  pL1$PERM_coeff<-as.numeric(as.character(pL1$PERM_coeff))
  pL1$DESeq_padj[is.na(pL1$DESeq_padj)]<-1
  pL1$otu<-rownames(pL1)
  pL1$color<-cut(pL1$DESeq_padj, breaks=c(-Inf, 0.05, Inf), labels=c("Significant", "Not Significant"))
  p=ggplot(pL1, aes(x=reorder(otu, -PERM_coeff), PERM_coeff, fill=color))
  p= p + geom_bar(stat="identity")
  p= p + theme_classic()
  p= p + coord_flip()
  p= p + theme(axis.text.x=element_text(hjust=1, size=10), axis.text.y=element_text(hjust=1, size=4))
  p= p + labs(title = title, x=xlab, y=ylab)
  return(p)
}

# DESeq alternative mean calculation
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
