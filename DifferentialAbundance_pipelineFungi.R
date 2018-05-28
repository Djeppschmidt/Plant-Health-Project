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

global.fungi.rare<-readRDS("")

# check to make sure taxa names are decent
taxa_names(global.fungi)

# subset

fungi.fsp.RRsoy.raw<-subset_samples(global.fungi, Location == "Beltsville" & genotype == "RR" & crop == "soy")
fungi.fsp.RRcorn.raw<-subset_samples(global.fungi, Location == "Beltsville" & genotype == "RR" & crop == "corn")
fungi.sv.RRsoy.raw<-subset_samples(global.fungi, Location == "Stoneville" & genotype == "RR" & crop == "soy")
fungi.sv.RRcorn.raw<-subset_samples(global.fungi, Location == "Stoneville" & genotype == "RR" & crop == "corn")

ungi.fsp.RRsoy.rare<-subset_samples(global.fungi.rare, Location == "Beltsville" & genotype == "RR" & crop == "soy")
fungi.fsp.RRcorn.rare<-subset_samples(global.fungi.rare, Location == "Beltsville" & genotype == "RR" & crop == "corn")
fungi.sv.RRsoy.rare<-subset_samples(global.fungi.rare, Location == "Stoneville" & genotype == "RR" & crop == "soy")
fungi.sv.RRcorn.rare<-subset_samples(global.fungi.rare, Location == "Stoneville" & genotype == "RR" & crop == "corn")

# vst transformations ####


fun.fsp.corn.Deseq.res<-diffDESeq(fungi.fsp.RRcorn.raw)
fun.fsp.soy.Deseq.res<-diffDESeq(fungi.fsp.RRsoy.raw)
fun.SV.corn.Deseq.res<-diffDESeq(fungi.sv.RRcorn.raw)
fun.SV.soy.Deseq.res<-diffDESeq(fungi.sv.RRsoy.raw)

fun.fsp.corn.DeseqR.res<-diffDESeq(fungi.fsp.RRcorn.rare)
fun.fsp.soy.DeseqR.res<-diffDESeq(fungi.fsp.RRsoy.rare)
fun.SV.corn.DeseqR.res<-diffDESeq(fungi.sv.RRcorn.rare)
fun.SV.soy.DeseqR.res<-diffDESeq(fungi.sv.RRsoy.rare)

# explore output... looking for spp order ####
head(rownames(fun.SV.corn.Deseq.res))
head(rownames(tax_table(fun.sv.RRcorn.Deseq.raw)))
identical(rownames(fun.SV.corn.Deseq.res),rownames(tax_table(fun.sv.RRcorn.Deseq.raw)))



# add p-value to tax table! ####
fun.fsp.RRcorn.Deseq.raw<-extractP(fun.fsp.RRcorn.Deseq.raw, fun.fsp.corn.Deseq.res)
fun.fsp.RRsoy.Deseq.raw<-extractP(fun.fsp.RRsoy.Deseq.raw, fun.fsp.soy.Deseq.res) #use this to test permanova!/color
fun.sv.RRcorn.Deseq.raw<-extractP(fun.sv.RRcorn.Deseq.raw,fun.SV.corn.Deseq.res)
fun.sv.RRsoy.Deseq.raw<-extractP(fun.sv.RRsoy.Deseq.raw, fun.SV.soy.Deseq.res) #rerun deseq var stab part


fun.fsp.RRcorn.Deseq.rare<-extractP(fun.fsp.RRcorn.Deseq.rare, fun.fsp.corn.DeseqR.res)
fun.fsp.RRsoy.Deseq.rare<-extractP(fun.fsp.RRsoy.Deseq.rare, fun.fsp.soy.DeseqR.res)
fun.sv.RRcorn.Deseq.rare<-extractP(fun.sv.RRcorn.Deseq.rare,fun.SV.corn.DeseqR.res)
fun.sv.RRsoy.Deseq.rare<-extractP(fun.sv.RRsoy.Deseq.rare, fun.SV.soy.DeseqR.res)

#check:
rank_names(tax_table(fun.fsp.RRsoy.Deseq.raw))
# transform sequence table to vst, remove zeros



# save new RDS for PERMANOVA pipeline input! ####
saveRDS(fun.fsp.RRsoy.Deseq.raw, file="/Users/maullabmacbookpro/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/Raw/beltsville_soy/Beltsville_soy_padj.RDS")

saveRDS(fun.fsp.RRcorn.Deseq.raw, file="/Users/maullabmacbookpro/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/Raw/beltsville_corn/Beltsville_corn_padj.RDS")

saveRDS(fun.SV.RRsoy.Deseq.raw, file="/Users/maullabmacbookpro/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/Raw/stoneville_soy/Stoneville_soy_padj.RDS")

saveRDS(fun.SV.RRcorn.Deseq.raw, file="/Users/maullabmacbookpro/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/Raw/stoneville_corn/Stoneville_corn_padj.RDS")

?saveRDS
getwd()
