### Preprocess and make phyloseq object ###

#start interactive session:
#srun --pty -p short -t 05:00:00 -n 10 -N 1 /bin/bash -l
getwd() #make sure in correct working directory

library(phyloseq)

sampleDat<-as.data.frame(as.matrix(read.csv("Master_SampleData_Table_PHPv3.csv")))

Head(sampleDat) #check import

com.Dat<-readRDS("allseq_nochim.rds")
com.Dat #check import

tax<-tax_table(com.Dat) # isolate tax table for later ... check syntax!

c.factor<-com.Dat$Correction #name QPCR correction factor to be applied to sequence counts

com.Dat<-as.data.frame(as.matrix(otu_table(com.Dat))) #prepare sequence counts for transformation

adj.comDat<-com.Dat*c.factor #do transformation

head(adj.comDat) #visually chech transformation

adj.comDat<-round(adj.comDat) #get rid of decimals/ convert to count data again

adj.comDat<-otu_table(adj.comDat, taxa_are_rows=#F?) #change based on output from above!!!

sampleDat<-sample_data(sampleDat)


PSadj.comDat<-phyloseq(adj.comDat, sampleDat)

saveRDS(, "QPCR_Corrected_totalPS.rds")
