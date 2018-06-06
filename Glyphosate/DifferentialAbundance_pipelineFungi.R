# fungal pipeline
library(ggplot2)
library(data.table)
library(phyloseq)
library(gplots)
library(VennDiagram)

require(devtools)
#install_version("vegan", version = "2.4-6", repos = "http://cran.us.r-project.org")
library(vegan)
library(DESeq2)


# import global dataset
global.fungi<-readRDS("/Users/maullabmacbookpro/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/Raw/global_fungi_phyloseq.rds")
global.fungi<-subset_taxa(global.fungi, Kingdom=="k__Fungi")

global.fungi.rare<-readRDS("/Users/maullabmacbookpro/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/Rarefied/fungi_RAREcomDat.rds")
global.fungi.rare<-subset_taxa(global.fungi.rare, Kingdom=="k__Fungi")
# check to make sure taxa names are decent
taxa_names(global.fungi)

# subset

fungi.fsp.RRsoy.raw<-subset_samples(global.fungi, Location == "Beltsville" & genotype == "RR" & crop == "soy")
fungi.fsp.RRcorn.raw<-subset_samples(global.fungi, Location == "Beltsville" & genotype == "RR" & crop == "corn")
fungi.sv.RRsoy.raw<-subset_samples(global.fungi, Location == "Stoneville" & genotype == "RR" & crop == "soy")
fungi.sv.RRcorn.raw<-subset_samples(global.fungi, Location == "Stoneville" & genotype == "RR" & crop == "corn")

fungi.fsp.RRsoy.rare<-subset_samples(global.fungi.rare, Location == "Beltsville" & genotype == "RR" & crop == "soy")
fungi.fsp.RRcorn.rare<-subset_samples(global.fungi.rare, Location == "Beltsville" & genotype == "RR" & crop == "corn")
fungi.sv.RRsoy.rare<-subset_samples(global.fungi.rare, Location == "Stoneville" & genotype == "RR" & crop == "soy")
fungi.sv.RRcorn.rare<-subset_samples(global.fungi.rare, Location == "Stoneville" & genotype == "RR" & crop == "corn")

# vst transformations ####


fun.fsp.corn.Deseq.res<-diffDESeq(fungi.fsp.RRcorn.raw)
fun.fsp.soy.Deseq.res<-diffDESeq(fungi.fsp.RRsoy.raw)
fun.SV.corn.Deseq.res<-diffDESeq(fungi.sv.RRcorn.raw)
fun.SV.soy.Deseq.res<-diffDESeq(fungi.sv.RRsoy.raw)

fun.fsp.corn.DeseqR.res<-diffDESeq(fungi.fsp.RRcorn.rare)#not run bc need to add correct metadata
fun.fsp.soy.DeseqR.res<-diffDESeq(fungi.fsp.RRsoy.rare)#not run
fun.SV.corn.DeseqR.res<-diffDESeq(fungi.sv.RRcorn.rare)#not run
fun.SV.soy.DeseqR.res<-diffDESeq(fungi.sv.RRsoy.rare)#not run

# explore output... looking for spp order ####
head(rownames(fun.SV.corn.Deseq.res))
head(rownames(tax_table(fun.sv.RRcorn.Deseq.raw)))
identical(rownames(fun.SV.corn.Deseq.res),rownames(tax_table(fun.sv.RRcorn.Deseq.raw)))


# permanova test ####

fun.fsp.Soy.fusarium<-makefusarium(fun.fsp.RRsoy.Deseq.raw)
anova.fsp.RRsoy.F<-permanova(fun.fsp.RRsoy.Deseq.raw, n=2.5, fun.fsp.Soy.fusarium)

fun.fsp.corn.fusarium<-makefusarium(fun.fsp.RRcorn.Deseq.raw)
anova.fsp.RRcorn.F<-permanova(fun.fsp.RRcorn.Deseq.raw, n=2.5, fun.fsp.corn.fusarium)

fun.sv.Soy.fusarium<-makefusarium(fun.sv.RRsoy.Deseq.raw)
anova.sv.RRsoy.F<-permanova2(fun.sv.RRsoy.Deseq.raw, n=2.5, fun.sv.Soy.fusarium)

fun.sv.corn.fusarium<-makefusarium(fun.sv.RRcorn.Deseq.raw)
anova.sv.RRcorn.F<-permanova2(fun.sv.RRcorn.Deseq.raw, n=2.5, fun.sv.corn.fusarium)



anova.fsp.RRsoy.FR<-permanova(fun.fsp.RRsoy.Deseq.rare, n=2.5, fun.fsp.Soy.fusarium)
anova.fsp.RRcorn.FR<-permanova(fun.fsp.RRcorn.Deseq.rare, n=2.5, fun.fsp.corn.fusarium)
anova.sv.RRsoy.FR<-permanova2(fun.sv.RRsoy.Deseq.rare, n=2.5, fun.sv.Soy.fusarium)
anova.sv.RRcorn.FR<-permanova2(fun.sv.RRcorn.Deseq.rare, n=2.5, fun.sv.corn.fusarium)


# add p-value to tax table! ####
fun.fsp.RRcorn.Deseq.raw<-extractP(fungi.fsp.RRcorn.raw, fun.fsp.corn.Deseq.res, anova.fsp.RRcorn.F)
fun.fsp.RRsoy.Deseq.raw<-extractP(fungi.fsp.RRsoy.raw, fun.fsp.soy.Deseq.res, anova.fsp.RRsoy.F)
fun.sv.RRcorn.Deseq.raw<-extractP(fungi.sv.RRcorn.raw,fun.SV.corn.Deseq.res, anova.sv.RRcorn.F)
fun.sv.RRsoy.Deseq.raw<-extractP(fungi.sv.RRsoy.raw, fun.SV.soy.Deseq.res, anova.sv.RRsoy.F)



fun.fsp.RRcorn.Deseq.rare<-extractP(fun.fsp.RRcorn.Deseq.rare, fun.fsp.corn.DeseqR.res)#not run
fun.fsp.RRsoy.Deseq.rare<-extractP(fun.fsp.RRsoy.Deseq.rare, fun.fsp.soy.DeseqR.res)#not run
fun.sv.RRcorn.Deseq.rare<-extractP(fun.sv.RRcorn.Deseq.rare,fun.SV.corn.DeseqR.res)#not run
fun.sv.RRsoy.Deseq.rare<-extractP(fun.sv.RRsoy.Deseq.rare, fun.SV.soy.DeseqR.res)#not run

#check:
rank_names(tax_table(fun.sv.RRsoy.Deseq.raw))
# transform sequence table to vst, remove zeros



# save new RDS for PERMANOVA pipeline input! ####
#saveRDS(fun.fsp.RRsoy.Deseq.raw, file="/Users/maullabmacbookpro/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/Raw/beltsville_soy/Beltsville_soy_padj.RDS")

#saveRDS(fun.fsp.RRcorn.Deseq.raw, file="/Users/maullabmacbookpro/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/Raw/beltsville_corn/Beltsville_corn_padj.RDS")

#saveRDS(fun.SV.RRsoy.Deseq.raw, file="/Users/maullabmacbookpro/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/Raw/stoneville_soy/Stoneville_soy_padj.RDS")

#saveRDS(fun.SV.RRcorn.Deseq.raw, file="/Users/maullabmacbookpro/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/Raw/stoneville_corn/Stoneville_corn_padj.RDS")

?saveRDS
getwd()
library(ggplot2)
# plotsss ####

p1<-plotCoeffF(fun.sv.RRsoy.Deseq.raw, title="Stoneville Soy Fusarium", xlab="Taxon", ylab="Coefficient Value")#test
p2<-plotCoeffF(fun.sv.RRcorn.Deseq.raw, title="Fungi Stoneville Corn Fusarium", xlab="Taxon", ylab="Coefficient Value")#test
p3<-plotCoeffF(fun.fsp.RRsoy.Deseq.raw, title="Fungi Beltsville Soy Fusarium", xlab="Taxon", ylab="Coefficient Value")
p4<-plotCoeffF(fun.fsp.RRcorn.Deseq.raw, title="Fungi Beltsville Corn Fusarium", xlab="Taxon", ylab="Coefficient Value")

multiplot(p1, p2, p3, p4, cols=2)
# scratch space ####
head(sample_data(fun.sv.RRcorn.Deseq.raw))
table(sample_data(fun.sv.RRcorn.Deseq.raw)$System.loc)
