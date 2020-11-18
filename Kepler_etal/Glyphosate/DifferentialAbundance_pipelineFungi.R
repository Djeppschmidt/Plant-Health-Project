# fungal pipeline

# vst transformations ####

fun.sv.RRcorn.Deseq.raw<-readRDS("~/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/Raw/Stoneville_corn/Stoneville_corn.rds")

fun.sv.RRsoy.Deseq.raw<-readRDS("~/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/Raw/Stoneville_soy/Stoneville_soy.rds")

fun.fsp.RRcorn.Deseq.raw<-readRDS("~/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/Raw/beltsville_corn/Beltsville_corn.rds")

fun.fsp.RRsoy.Deseq.raw<-readRDS("~/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/Raw/beltsville_soy/Beltsville_soy.rds")


fun.fsp.corn.Deseq.res<-diffDESeq(fun.fsp.RRcorn.Deseq.raw)
fun.fsp.soy.Deseq.res<-diffDESeq(fun.fsp.RRsoy.Deseq.raw)
fun.SV.corn.Deseq.res<-diffDESeq(fun.sv.RRcorn.Deseq.raw)
fun.SV.soy.Deseq.res<-diffDESeq(fun.sv.RRsoy.Deseq.raw)

# explore output... looking for spp order ####
head(rownames(fun.SV.corn.Deseq.res))
head(rownames(tax_table(fun.sv.RRcorn.Deseq.raw)))
identical(rownames(fun.SV.corn.Deseq.res),rownames(tax_table(fun.sv.RRcorn.Deseq.raw)))



# add p-value to tax table! ####
fun.fsp.RRcorn.Deseq.raw<-extractP(fun.fsp.RRcorn.Deseq.raw, fun.fsp.corn.Deseq.res)
fun.fsp.RRsoy.Deseq.raw<-extractP(fun.fsp.RRsoy.Deseq.raw, fun.fsp.soy.Deseq.res) #use this to test permanova!/color
fun.sv.RRcorn.Deseq.raw<-extractP(fun.sv.RRcorn.Deseq.raw,fun.SV.corn.Deseq.res)
fun.sv.RRsoy.Deseq.raw<-extractP(fun.sv.RRsoy.Deseq.raw, fun.SV.soy.Deseq.res) #rerun deseq var stab part


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
