### Preprocess and make phyloseq object ###
#upload supplemental data
#scp /Users/maullabmacbookpro/Documents/GitHub/Plant-Health-Project/Master_SampleData_Table_PHPv3_2.csv dietrich.schmidt@scinet-login.bioteam.net:/project/php-bacarc/new_preproc/sequence_variants2/Master_SampleData_Table_PHPv3_2.csv

#start interactive session:
#srun --pty -p short -t 05:00:00 -n 10 -N 1 /bin/bash -l
#module load r
#R
getwd() #make sure in correct working directory

library(phyloseq)
### CHECK FIRST IMPORT STEP!! SOMETHING WRONG WHEN MULTIPLYING!!! ### <- this is because urbana is missing sample names!! ugh!
##I think this ^ is resolved...
sampleDat<-as.data.frame(read.csv("Master_SampleData_Table_PHPv3_2.csv", stringsAsFactors=FALSE))
c.factor<-as.numeric(sampleDat$Correction)
c.factor[c.factor==NA]<-1
rownames(sampleDat)<-sampleDat$X.SampleID
sampleDat<-sample_data(sampleDat)

head(sampleDat) #check import

com.Dat<-readRDS("allseq_nochim.rds")
dim(com.Dat) #check import

tax<-readRDS("tax_training.rds") # import tax assingments
dim(tax) #check import
tax<-tax_table(tax)
 #name QPCR correction factor to be applied to sequence counts

com.Dat<-as.data.frame(com.Dat) #prepare sequence counts for transformation
com.Dat[1:3,1:3]

sampleDat[sampleDat$Correction==NA]<-1

#adj.comDat<-com.Dat*sampleDat$Correction #do transformation/ test case not working
adj.comDat<-com.Dat*c.factor

head(adj.comDat[1:3,1:3]) #visually chech transformation

adj.comDat<-round(adj.comDat) #get rid of decimals/ convert to count data again

adj.comDat<-otu_table(adj.comDat, taxa_are_rows=F) #change based on output from above!!!

sampleDat<-sample_data(sampleDat)


PSadj.comDat<-merge_phyloseq(adj.comDat, sampleDat, tax)

saveRDS(PSadj.comDat, "QPCR_Corrected_totalPS.rds")
