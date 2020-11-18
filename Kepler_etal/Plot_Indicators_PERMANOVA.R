# ggplot2 of indicator spp ####

# Bacteria ####

# fsp Corn ####
p<-ggplot(anova.fsp.RRcorn$sppInt, aes(x=reorder(factor.indics., -Glyphosphate_Treatment1.Sampling_date1), Glyphosphate_Treatment1.Sampling_date1))
p+geom_bar(stat="identity")+theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))+coord_flip() + #color by p-value cutoff from DESeq2



# fsp Soy ####

p<-ggplot(anova.fsp.RRsoy$sppInt, aes(x=reorder(factor.indics., -Glyphosphate_Treatment1.Sampling_date1), Glyphosphate_Treatment1.Sampling_date1))
p+geom_bar(stat="identity")+theme_classic() + coord_flip() + theme(axis.text.x=element_text(angle=90, hjust=1))

head(anova.fsp.RRsoy.F$sppInt)

# sv Corn ####
p<-ggplot(anova.sv.RRcorn$sppInt, aes(x=reorder(factor.indics., -Glyphosphate_Treatment1.Sampling_date1), Glyphosphate_Treatment1.Sampling_date1))
p+geom_bar(stat="identity") + theme_classic() + coord_flip() + theme(axis.text.x=element_text(angle=90, hjust=1))

# sv Soy ####

p<-ggplot(anova.sv.RRsoy$sppInt, aes(x=reorder(factor.indics., -Glyphosphate_Treatment1.Sampling_date1), Glyphosphate_Treatment1.Sampling_date1))
p+geom_bar(stat="identity")+theme_classic() + coord_flip() + theme(axis.text.x=element_text(angle=90, hjust=1))


# Fungi ####


# fsp Soy ####

#for variance stabilized dataset; doing raw permanova first
#fun.fsp.vst<-readRDS("/Users/maullabmacbookpro/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/VST/beltsville_soy/Beltsville_soy_vst.rds")
#dim(fun.fsp.vst)
#fun.fsp.vst<-otu_table(fun.fsp.vst, taxa_are_rows = TRUE)

#otu_table(fun.fsp.RRsoy.Deseq.vst)<-fun.fsp.vst

#rank_names(tax_table(fun.fsp.RRsoy.Deseq.vst))

#as.data.frame(as.matrix(tax_table(fun.fsp.RRsoy.Deseq.vst)))$Genus

#fun.fsp.RRsoy.Deseq.vst<-extractP(fun.fsp.RRsoy.Deseq.vst, fun.fsp.soy.Deseq.res)

fun.fsp.Soy.fusarium<-makefusarium(fun.fsp.RRsoy.Deseq.raw)

anova.fsp.RRsoy<-permanova(fun.fsp.RRsoy.Deseq.raw, n=2.5, fun.fsp.Soy.fusarium)

seqs<-taxa_names(fun.fsp.RRsoy.Deseq.raw)
name<-paste0("Seq", seq(ntaxa(fun.fsp.RRsoy.Deseq.raw)))
nameKey<-data.frame(name, seqs)

taxa_names(fun.fsp.RRsoy.Deseq.raw)<-name

head(anova.fsp.RRsoy$sppInt)
test<-subset_taxa(fun.sv.RRsoy.Deseq.raw, Genus=="g__Fusarium")
test2<-as.data.frame(tax_table(test))$DESeq_padj
test3<-psmelt(test)
head(test3)

p<-ggplot(anova.fsp.RRsoy.F$sppInt, aes(x=reorder(factor.indics., -Glyphosphate_Treatment1.Sampling_date1), Glyphosphate_Treatment1.Sampling_date1))
p+geom_bar(stat="identity")+theme_classic() + coord_flip() + theme(axis.text.x=element_text(angle=90, hjust=1)) + scale_fill_manual("legend", values=c("sig"))

t<-plot_bar(subset_taxa(fun.sv.RRsoy.Deseq.raw, Genus=="g__Fusarium"), fill="DESeq_padj")

