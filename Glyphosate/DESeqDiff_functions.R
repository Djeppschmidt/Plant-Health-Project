#Fungal differential abundance test ->-> DESeq2

# differential abundance function

diffDESeq<-function(ps){
  L1 <-phyloseq_to_deseq2(ps, ~ group + Sampling_date + group:Sampling_date)
  
  diagdds = DESeq(L1, test="Wald", fitType="parametric")
  
  res = results(diagdds)
  
 return(res) 
}


# extract values, put in tax table ####
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

table(tax_table(fun.sv.RRcorn.Deseq.raw))
rank_names(fun.sv.RRcorn.Deseq.raw)


# test space ####
NAME<-"Input.file.rds"



fun.sv.RRcorn.Deseq.raw<-readRDS("~/Desktop/PhD/PHP/StatPlot_Outputs/Fungi/Raw/beltsville_corn/Beltsville_corn.rds")

fun.SV.corn.Deseq.res<-diffDESeq(fun.sv.RRcorn.Deseq.raw)

test<-extractP(fun.fsp.RRsoy.Deseq.raw,fun.fsp.soy.Deseq.res) # worked!


identical(rownames(tax_table(fun.fsp.RRsoy.Deseq.raw)), rownames(fun.fsp.soy.Deseq.res))
head(rownames(fun.SV.soy.Deseq.res))
rank_names(tax_table(test)) ### exellent! test worked, added column of p-values!
