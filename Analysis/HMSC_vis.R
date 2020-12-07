# visualization and data analysis for PHP HMSC
# To do on local machine:
# extract + associations
# extract - associations
# filter for large effects in each
# filter for bacterial-fungal interactions [examine seperately unique interactions]
# Construct + / - association networks
# Identify taxa that are central (either bacteria or fungi)
# identify hubs
# Explore each of the hubs and central taxa individually
# (should the networks include taxa that only strongly associate with others of same type? / compare levels between networks? are there more fungi vs bacteria than within each group?)
# for local computer
# Make positive and negative association networks
# collect descriptive stats on networks ( can get from model 1)
# Run stats on difference among networks in farming systems
setwd()


# to do list:
# construct taxa metadata structure             [_]
# separate bac, fun and fun-bac interactions    [_]
# separate + / - interactions                   [_]
# get metrics for each network                  [_]
# make network using igraph                     [_]

# correlation structure: use logical 0 = false; 1=True thus: T+T=2 T+F=1 F+F=0 (or just code directly as logical)
# need info for network Q2: ID of taxa (fungi vs bacteria), Effect of factor on taxa; Delta interaction
# association network for bacteria-bacteria ([,1]=b & [,2]=b)*mean
# 
# 

# imaginary run 1 without effect of farming system 
# each of these are stored as matrices
# need to make a matching matrix of b-b and b-f and f-f interactions for each ...
# test case:
# 4X4 matrix

mat<-matrix(nrow=4, ncol=4)
rownames(mat)<-c("one", "two", "three", "four")
colnames(mat)<-c("one", "two", "three", "four")
grouping<-c(1,1,2,2)
mat[grouping==1, grouping==1]<-FALSE
mat[grouping==2, grouping==1]<-TRUE
mat[grouping==1, grouping==2]<-TRUE
mat[grouping==2, grouping==2]<-FALSE
mat

# test case 2:

identical(colnames(t.ct$LogPoi$OmegaCor.lp[[1]]$mean), taxa_names(CT.101c)) # true
identical(taxa_names(CT.101c), rownames(tax_table(CT.101c)))# true
# therefore order is preserved ...

function(x){
  net<-x$LogPoi$OmegaCor.lp[[1]]$mean
  BBnet<-net
  BBnet<-BBnet[BBnet,]
}

# make BB table

# make FF table

# make B-F table
tax_table(CT.101c)[,"Kingdom"]

# now make an index of taxon kingdom
s<-tax_table(CT.101c)[,1]
s[s=="Bacteria"]<-0
s
CT$LogPoi$OmegaCor.lp[[3]]$support[1,1] # support (significance)
CT$LogPoi$OmegaCor.lp[[3]]$mean # mean (effect size)

# imaginary run 2 with effect of farming system
CT2$LogPoi$OmegaCor.lp[[3]]$support # support (significance)
CT2$LogPoi$OmegaCor.lp[[3]]$mean # mean (effect size)

if(names(CT2))

View(as.data.frame(as.matrix(((t.ct$LogPoi$OmegaCor.lp[[2]]$support>0.95)+(t.ct$LogPoi$OmegaCor.lp[[2]]$support<(1-0.95))>0)*t.ct$LogPoi$OmegaCor.lp[[2]]$mean)))
# scratch
# probit model (update to reflect the actual analysis?)
RH.Ydat<-as.data.frame(as.matrix(otu_table(RH)))# dim(b.Ydat)
RH.XDat1<-as.data.frame(as.matrix(sample_data(RH)), stringsAsFactors = TRUE)
RH.XDat1<-b.XDat1[,c()] # this needs to be updated!!
rownames(b.XDat1)<-c(1:nrow(b.XDat1))
RH.XDat1$Sample<-as.factor(c(1:nrow(b.XDat1)))
RH.XFormula1= ~Glyphosphate_Treatment + genotype + System.loc + Loc_plot_ID + crop + pH + OM... # changed since last time
RH.XFormula2= ~Glyphosphate_Treatment + genotype +System.loc + Loc_plot_ID + crop + pH + OM...+ NLFA.AM.Fungi
studyDesign = data.frame("Sample"=b.XDat1$Sample,"Loc_plot_ID"=b.XDat1$Loc_plot_ID, "Sampling_date"=b.XDat1$Sampling_date)
rL1 <- HmscRandomLevel(units=studyDesign$Sample)
rL2 <- HmscRandomLevel(units=studyDesign$Loc_plot_ID) # set random level to loc_plot_ID
rL3 <- HmscRandomLevel(units=studyDesign$Sampling_date)
rL4 <- HmscRandomLevel(units=studyDesign$)# second random level is sampling 
# final model should have pH and organic matter as random levels
# first do probit model ####
b.m1<-Hmsc(Y=as.matrix(t(b.Ydat)), XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="probit")
b.m1<-sampleMcmc(b.m.lp, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)

OmegaCor.lp=computeAssociations(bf.lp)
probit3$OmegaCor.probit=computeAssociations(probit1$bfm1)
probit3$preds<-computePredictedValues(probit1$bfmq)
probit3$MF<-evaluateModelFit(probit1$bgm1, predY=probit1$preds)
probit3$VP<-computeVariancePartitioning(probit1$bfm1, group=c(), groupnames=c())
probit3$postBeta<-getPostEstimate(probit1$bfm1, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal
probit3$GradientNT.AM<-constructGradient(bf.m1, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "NT-MD")))
probit3$GradientCT.AM<-constructGradient(bf.lp, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "CT-MD")))
probit3$GradientO3.AM<-constructGradient(bf.lp, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "Org_3")))
probit13$GradientO6.AM<-constructGradient(bf.lp, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "Org_6")))
probit3$predAMGradientNT<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)# rename things!!
probit3$predAMGradientCT<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)
probit3$predAMGradientO3<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)
probit3$predAMGradientO6<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)

Y2<-as.matrix(t(b.Ydat))
Y2[Y2==0]<-NA
# lognormal poisson:
b.m1<-Hmsc(Y=as.matrix(t(b.Ydat)), XData=b.XDat, XFormula=b.XFormula, studyDesign=studyDesign, ranLevels=list("Sample"=rL1, "Loc_plot_ID"=rL2, "Sampling_date"=rL3), distr="lognormal poisson")
b.m1<-sampleMcmc(b.m.lp, thin=2, samples=1000, transient=500, nChains=2, nParallel=2, verbose=100)

OmegaCor.lp=computeAssociations(bf.lp)
probit1$OmegaCor.probit=computeAssociations(probit1$bfm1)
probit1$preds<-computePredictedValues(probit1$bfmq)
probit1$MF<-evaluateModelFit(probit1$bgm1, predY=probit1$preds)
probit1$VP<-computeVariancePartitioning(probit1$bfm1, group=c(), groupnames=c())
probit1$postBeta<-getPostEstimate(probit1$bfm1, parName="Beta") #beta is species abundance ; gamma is traits; rho is phylogenetic signal
probit3$GradientNT.AM<-constructGradient(bf.lp, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "NT-MD")))
probit3$GradientCT.AM<-constructGradient(bf.lp, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "CT-MD")))
probit3$GradientO3.AM<-constructGradient(bf.lp, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "Org_3")))
probit13$GradientO6.AM<-constructGradient(bf.lp, focalVariable="NLFA.AM.Fungi", non.focalVariables=list("OM..."=list(1),"pH"=list(1), "System.loc"=list(3, "Org_6")))
probit3$predAMGradientNT<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)# rename things!!
probit3$predAMGradientCT<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)
probit3$predAMGradientO3<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)
probit3$predAMGradientO6<-predict(probit3$GradientNT.AM, XData=probit3$GradientNT.AM, studyDesign= , ranLevels=, expected=TRUE)