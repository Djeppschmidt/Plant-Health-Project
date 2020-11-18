# variance test
# example input: variance.test(fsp.RRcorn, sample_data(fsp.RRcorn)$Sampling_date)
# I tried to get this to work by doing the parsing within the function, but wouldn't read properly... always got an error. Works this way now...

variance.test<-function(ps, group){
  require(phyloseq)
  require(vegan)
  dist<-vegdist(otu_table(ps), method="euclidean")
  a<-betadisper(dist, group)
  b<-anova(a)
  out<-list("var"=a, "anova"=b)
}

extract.varPre<-function(x){mean(x$var$distance[x$var$group=="pre"])}
extract.varPost<-function(x){mean(x$var$distance[x$var$group=="post"])}


Var.summ<-function(x){
  pre<-lapply(x, extract.varPre)
  post<-lapply(x, extract.varPost)
  df<-data.frame(pre,post)
  colnames(df)<-c("pre", "post")
  rownames(df)<-names(x)
  df
}

# scratch space ####

Var.summ(list(var))
extract.varPre(var)
extract.varPost(var)
lapply(x, extract.varPre)
lapply(x, extract.varPost)
extract.varPre<-function(x){mean(x$var$distance[x$var$group=="pre"])}
extract.varPost<-function(x){mean(x$var$distance[x$var$group=="post"])}
mean(var$var$distances[var$var$group=="pre"])
mean(var$var$distances[var$var$group=="post"])


var<-variance.test(fsp.RRcorn, sample_data(fsp.RRcorn)$Sampling_date)

dist<-vegdist(otu_table(fsp.RRcorn), method="euclidean")
a<-betadisper(dist, sample_data(fsp.RRcorn)$Sampling_date)
a
?vegdist
?betadisper
length(sample_data(fsp.RRcorn)$Sampling_date)
fsp.RRcorn

var$anova$`Pr(>F)`

extract.var<-function(x, t){mean(x$var$group==t)}


mean(var$var$distances[var$var$group=="pre"])
mean(var$var$distances[var$var$group=="post"])


variance.test<-function(ps, group){
  require(phyloseq)
  require(vegan)
  dist<-vegdist(otu_table(ps), method="euclidean")
  a<-betadisper(dist, sample_data(ps)$group)
  b<-anova(a)
  out<-list("var"=a, "anova"=b)
}
