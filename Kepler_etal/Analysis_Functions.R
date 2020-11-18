#analysis functions

#model<-formula(y~a*b+c+d+f)
GenStat<-function(Input, model, Output){
  set.seed(397652)
  dir.create(Output)
  require(phyloseq)
  this1<-readRDS(Input) #import
  #alpha diversity
  pdf(file=paste0(Output,"DivPlot.pdf",sep=""))

  plot_richness(this1, x="Sampling_date", title="Sampling_date")
  plot_richness(this1, x="Glyphosphate_Treatment", title="Glyphosphate_Treatment")
  plot_richness(this1, x="System.loc", title="System.loc")
  plot_richness(this1, x="crop", title="crop")
  plot_richness(this1, x="Soil_Zone", title="Soil_Zone")

  dev.off()
  #define factors
  a<-sample_data(this1)$Sampling_date
  b<-sample_data(this1)$Glyphosphate_Treatment
  c<-sample_data(this1)$System.loc
  d<-sample_data(this1)$crop
  f<-sample_data(this1)$Soil_Zone
  y<-unlist(estimate_richness(this1, split=TRUE, measures= "Observed"), recursive=T, use.names=T)
  #use factors in model
  sink(paste0(Output,DivStats.txt,sep=""), append=TRUE)


  summary(aov(model))
  y<-phyloseq::distance(this1, method="bray")
  sampledf<-data.frame(sample_data(this1))
  adonis(model, data=sampledf)
  beta<-betadisper(y, a)#test dispersions of all groups
  permutest(beta)

  beta<-betadisper(y, b)
  permutest(beta)

  beta<-betadisper(y, c)
  permutest(beta)

  beta<-betadisper(y, d)
  permutest(beta)

  beta<-betadisper(y, f)
  permutest(beta)

  dev.off()
  #PERMANOVA
  #ORDINATION
  pdf(file=paste0(Output,"ordination.pdf",sep=""))
  nmds<-ordinate(this1,method="NMDS",distance="bray")
  plot_ordination(physeq=this1,ordination=nmds,color="System.loc", shape="Soil_Zone")
  dev.off()
  #add more ordination types...

}
