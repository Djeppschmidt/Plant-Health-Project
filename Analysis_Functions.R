#analysis functions


GenStat<-function(Input, model, Output){
  require(phyloseq)
  this1<-readRDS(Input) #import
  #alpha diversity
  pdf(file=paste0(Output,DivPlot.pdf,sep="/"), append=TRUE)

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
  sink(paste0(Output,DivStats.txt,sep="/"), append=TRUE)

  
summary(aov(model))

dev.off()
  #PERMANOVA
  #ORDINATION


}
