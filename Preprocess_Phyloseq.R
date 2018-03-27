### Preprocess and make phyloseq object ###
#upload supplemental data
#scp /Users/maullabmacbookpro/Documents/GitHub/Plant-Health-Project/Master_SampleData_Table_PHPv3_2.csv dietrich.schmidt@scinet-login.bioteam.net:/project/php-bacarc/new_preproc/sequence_variants2/Master_SampleData_Table_PHPv3_2.csv

#start interactive session:
#srun --pty -p short -t 05:00:00 -n 10 -N 1 /bin/bash -l
#module load r
#R
getwd() #make sure in correct working directory

set.seed(397652)
library(phyloseq)
### CHECK FIRST IMPORT STEP!! SOMETHING WRONG WHEN MULTIPLYING!!! ### <- this is because urbana is missing sample names!! ugh!
##I think this ^ is resolved...
sampleDat<-as.data.frame(read.csv("Master_SampleData_Table_PHPv3_2.csv", stringsAsFactors=FALSE))
#c.factor<-as.numeric(sampleDat$Correction) # needs to happen after rarefaction
#c.factor[c.factor==NA]<-1
rownames(sampleDat)<-sampleDat$X.SampleID
sampleDat<-sample_data(sampleDat)

head(sampleDat) #check import

com.Dat<-readRDS("allseq_nochim.rds")
dim(com.Dat) #check import
com.Dat<-otu_table(com.Dat)
#first rarefy the community data!!

tax<-readRDS("tax_training.rds") # import tax assingments
dim(tax) #check import
tax<-tax_table(tax)
 #name QPCR correction factor to be applied to sequence counts

ps<-merge_phyloseq(com.Dat, sampleDat, tax)

#run this after making phyloseq object
com.rare<-rarefy_even_depth(ps, sample.size = 20000, rngseed = FALSE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

saveRDS(com.rare, "rarefied_com.rds")


#cd /project/php-bacarc/new_preproc/sequence_variants2/
#com.rare<-readRDS("rarefied_com.rds")
 #name QPCR correction factor to be applied to sequence counts


rare.OTU<-as.data.frame(as.matrix(otu_table(com.rare))) #prepare sequence counts for transformation
rare.OTU[1:3,1:3]

c.factor<-as.numeric(sample_data(com.rare)$Correction)
c.factor[c.factor==NA]<-1
#adj.comDat<-com.Dat*sampleDat$Correction #do transformation/ test case not working
adj.rareCom<-rare.OTU*c.factor

head(adj.rareCom[1:3,1:3]) #visually check transformation

adj.rareCom<-round(adj.rareCom) #get rid of decimals/ convert to count data again

adj.rareCom<-otu_table(adj.rareCom, taxa_are_rows=F) #change based on output from above!!!

#sampleDat<-sample_data(sampleDat)

PSadj.comDat<-com.rare
otu_table(PSadj.comDat)<-adj.rareCom
otu_table(PSadj.comDat)[1:3,1:3]#check that OTU table was added correctly

saveRDS(PSadj.comDat, "QPCR_Corrected_totalPS.rds")


#split into bac vs arch
#com<-readRDS("QPCR_Corrected_totalPS.rds")
Arc<-subset_taxa(com, Kingdom=="Archaea")
Bac<-subset_taxa(com, Kingdom=="Bacteria")

saveRDS(Arc, "ArchaealComTest.rds")
saveRDS(Bac, "BacterialComTest.rds")

#split by site location
