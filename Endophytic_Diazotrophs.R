###
###
### Testing frankia
###
###

bacarc<-readRDS("~/Desktop/PhD/PHP/bacarc_RAWcomDat.rds")
#load a rarefied version!

Q.rare<-readRDS("~/Desktop/PhD/PHP/bacarc_Qrare.rds")



?order
?sort
genera<-unique(as.data.frame(as.matrix(tax_table(bacarc)))$Genus)

order(as.character(genera))
sort(genera)

?subset_taxa
#Burkholderia-Paraburkholderia
Burk<-subset_taxa(Q.rare, Genus=="Burkholderia-Paraburkholderia")
Burk.R<-estimate_richness(Burk, measures="Observed")
plot_richness(Burk)
?estimate_richness
?plot_richness
sample_names(Burk)
plot(get_taxa(Burk, "128.S10"))
head(tax_table(Burk))
plot_bar(Burk, x="System.loc")

Burk.B<-subset_samples(Q.rare, Location=="Beltsville")
plot_bar(Burk.B, x="System.loc")
plot_bar(Burk.B, x="crop")
#Klebsiella 
Kleb<-subset_taxa(Q.rare, Genus=="Klebsiella")
plot_bar(Kleb, x="System.loc")
plot_bar(Kleb, x="crop")

#Enterobacter
Ent<-subset_taxa(Q.rare, Genus=="Enterobacter")
plot_bar(Ent, x="System.loc")
plot_bar(Ent, x="crop")


#Azospirillum
Azo<-subset_taxa(Q.rare, Genus=="Azospirillum")
plot_bar(Azo, x="System.loc")
plot_bar(Azo, x="crop")

#Bacillus
Baci<-subset_taxa(Q.rare, Genus=="Bacillus")
plot_bar(Baci, x="System.loc")
plot_bar(Baci, x="crop")

#Derxia
Derx<-subset_taxa(Q.rare, Genus=="Derxia")
plot_bar(Derx, x="System.loc")
plot_bar(Derx, x="crop")

#Erwinia doesn't exist in our dataset
#Erw<-subset_taxa(Q.rare, Genus=="Erwinia")

#Herbispirilum
Herb<-subset_taxa(Q.rare, Genus=="Herbaspirillum")
plot_bar(Herb, x="System.loc")
plot_bar(Herb, x="crop")

#pseudomonas
Pse<-subset_taxa(Q.rare, Genus=="Pseudomonas")
plot_bar(Pse, x="System.loc")
plot_bar(Pse, x="crop")
