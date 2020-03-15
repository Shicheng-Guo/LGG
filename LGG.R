BiocManager::install("ChAMP") 
BiocManager::install("doParallel") 
BiocManager::install("benchmarkme") 
BiocManager::install("DMRcate") 
BiocManager::install("minfiData")
BiocManager::install("missMethyl")
BiocManager::install("minfiData")

library("ChAMP")
library("ggplot2")
require("minfi")
library("knitr")
library("limma")
library("minfi")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("IlluminaHumanMethylation450kmanifest")
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library("IlluminaHumanMethylationEPICmanifest")
library("RColorBrewer")
library("missMethyl")
library("minfiData")
library("Gviz")
library("DMRcate")
library("stringr")
library("ChAMP")
library("doParallel")

Dir="C:/Users/Schrodi Lab/Documents/GitHub/LGG/extdata"
setwd("C:/Users/Schrodi Lab/Documents/GitHub/LGG/extdata")
set.seed(11)

targets <- read.metharray.sheet(Dir)
RGSet <- read.metharray.exp(targets = targets)
phenoData <- pData(RGSet)
manifest <- getManifest(RGSet)
head(getProbeInfo(manifest))
myNormalRGSet<-preprocessFunnorm(RGSet, nPCs=4, sex = NULL, bgCorr = TRUE,dyeCorr = TRUE, keepCN = TRUE, ratioConvert = TRUE,verbose = TRUE)
predictedSex <- getSex(myNormalRGSet, cutoff = -2)$predictedSex
predictedSex==targets$Gender

myLoad <- champ.load(Dir,filterBeads=TRUE,arraytype="450k")
champ.QC(beta = myLoad$beta,resultsDir="../CHAMP_Raw_QCimages/")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450k",cores=1)
champ.QC(beta = myNorm,pheno=myLoad$pd$Sample_Group,resultsDir="./CHAMP_Norm_QCimages/")

myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide"))
##########################################################################
# don't use all the cores which will easily be killed by system
detectCores()
seed=sample(seq(1,10000,by=1),1)
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450k",cores=1)
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$pureG3,compare.group=c("Case","Control"),arraytype="450k")
write.table(myDMP,file=paste("AtrialFibrillation.",seed,".24Case4Control.myDMP.txt",sep=""),col.names = NA,row.names = T,quote=F,sep="\t")


myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group,arraytype="450k")
write.table(myDMP,file=paste("AtrialFibrillation.",seed,".CaseControl.myDMP.txt",sep=""),col.names = NA,row.names = T,quote=F,sep="\t")

myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Young_Old,arraytype="450k",adjPVal = 0.1)
write.table(myDMP,file="AtrialFibrillation.YoungOld.myDMP.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Pmr_Enrollment,arraytype="450k",adjPVal = 0.1)
write.table(myDMP,file="AtrialFibrillation.Enrollment.myDMP.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,method="Bumphunter",arraytype="450k",minProbes=2,cores=6,maxGap=1000)
write.table(myDMR,file="AtrialFibrillation.CaseControl.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Young_Old,method="Bumphunter",arraytype="450k",minProbes=2,cores=6,maxGap=1000)
write.table(myDMR,file="AtrialFibrillation.YoungOld.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Pmr_Enrollment,method="Bumphunter",arraytype="450k",minProbes=2,cores=6,maxGap=1000)
write.table(myDMR,file="AtrialFibrillation.Enrollment.myDMR.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="450k",cores=6)
write.table(myBlock$Block,file="AtrialFibrillation.CaseControl.myblock.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="450k",cores=6)
write.table(myBlock$Block,file="AtrialFibrillation.YoungOld.myblock.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="450k",cores=6)
write.table(myBlock$Block,file="AtrialFibrillation.Enrollment.myblock.txt",col.names = NA,row.names = T,quote=F,sep="\t")

myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="450k",adjPval=0.05, method="fisher")
myebayGSEA <- champ.ebGSEA(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="450k")
write.table(myebayGSEA,file="AtrialFibrillation.myebayGSEA.txt",col.names = NA,row.names = T,quote=F,sep="\t")
myEpiMod <- champ.EpiMod(beta=myNorm,pheno=myLoad$pd$Sample_Group)
myCNA <- champ.CNA(intensity=myLoad$intensity,pheno=myLoad$pd$Sample_Group)
myRefBase1 <- champ.refbase(beta=myNorm,arraytype="450k")
myRefBase2 <- champ.refbase(beta=myNorm,arraytype="450K")
