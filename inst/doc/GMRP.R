### R code from vignette source 'GMRP.Rnw'

###################################################
### code chunk number 1: <style-Sweave
###################################################
BiocStyle::latex()

#library("knitr")
#opts_chunk$set(tidy=FALSE,dev="pdf",fig.show="hide",
 #              fig.width=4,fig.height=4.5,
  #             dpi=300,# increase dpi to avoid pixelised pngs
  #             message=FALSE)
               
###################################################
### code chunk number 2: GMRP.Rnw:40-43
###################################################
set.seed(102)
options(width = 90)


###################################################
### code chunk number 3: GMRP.Rnw:45-47
###################################################
library(GMRP)

###################################################
### code chunk number 4: fmerge.Rnw:56-72
###################################################
data1 <- matrix(NA, 20, 4)
data2 <- matrix(NA, 30, 7)
SNPID1 <- paste("rs", seq(1:20), sep="")
SNPID2 <- paste("rs", seq(1:30), sep="")
data1[,1:4] <- c(round(runif(20), 4), round(runif(20), 4), round(runif(20), 4), round(runif(20), 4))
data2[,1:4] <- c(round(runif(30), 4),round(runif(30), 4), round(runif(30), 4), round(runif(30), 4))
data2[,5:7] <- c(round(seq(30)*runif(30), 4), round(seq(30)*runif(30), 4), seq(30))
data1 <- cbind(SNPID1, as.data.frame(data1))
data2 <- cbind(SNPID2, as.data.frame(data2))
dim(data1)
dim(data2)
colnames(data1) <- c("SNP", "var1", "var2", "var3", "var4")
colnames(data2) <- c("SNP", "var1", "var2", "var3", "var4", "V1", "V2", "V3")
data1<-DataFrame(data1)
data2<-DataFrame(data2)
data12 <- fmerge(fl1=data1, fl2=data2, ID1="SNP", ID2="SNP", A=".dat1", B=".dat2", method="No")

###################################################
### code chunk number 5: disease data load:90-93
###################################################
data(cad.data)
#cad <- DataFrame(cad.data)
cad<-cad.data
head(cad)
###################################################
### code chunk number 6:lipid data load: 95-99
###################################################
data(lpd.data)
#lpd <- DataFrame(lpd.data)
lpd<-lpd.data
head(lpd)

###################################################
### code chunk number 7: choose SNPs: 127-152
###################################################

###################################################
# Step1: calculate pvj: 129-136
###################################################

pvalue.LDL <- lpd$P.value.LDL
pvalue.HDL <-lpd$P.value.HDL
pvalue.TG <- lpd$P.value.TG
pvalue.TC <- lpd$P.value.TC
pv <- cbind(pvalue.LDL, pvalue.HDL, pvalue.TG, pvalue.TC)
pvj <- apply(pv, 1, min)

###################################################
# Step2: retrieve causal variables from data: 139-145
###################################################

beta.LDL <- lpd$beta.LDL
beta.HDL <- lpd$beta.HDL
beta.TG <- lpd$beta.TG
beta.TC <- lpd$beta.TC
beta <- cbind(beta.LDL, beta.HDL, beta.TG, beta.TC)

###################################################
#     Step3: construct a matrix for allele1: 148-154
###################################################

a1.LDL <- lpd$A1.LDL
a1.HDL <- lpd$A1.HDL
a1.TG <- lpd$A1.TG
a1.TC <- lpd$A1.TC
alle1 <- cbind(a1.LDL, a1.HDL, a1.TG, a1.TC)

###################################################
# Step4:  calculate pcj:157-165
###################################################

N.LDL <- lpd$N.LDL
N.HDL <- lpd$N.HDL
N.TG <- lpd$N.TG
N.TC <- lpd$N.TC
ss <- cbind(N.LDL, N.HDL, N.TG, N.TC)
sm <- apply(ss,1,sum)
pcj <- round(sm/max(sm), 6)

################################################### 
# Step5: Construct matrix for frequency of allele1:168-174
###################################################

freq.LDL<-lpd$Freq.A1.1000G.EUR.LDL
freq.HDL<-lpd$Freq.A1.1000G.EUR.HDL
freq.TG<-lpd$Freq.A1.1000G.EUR.TG
freq.TC<-lpd$Freq.A1.1000G.EUR.TC
freq<-cbind(freq.LDL,freq.HDL,freq.TG,freq.TC)

###################################################
# Step6: construct matrix for sd of each causal variable: 177-183 
###################################################

sd.LDL <- rep(37.42, length(pvj))
sd.HDL <- rep(14.87, length(pvj))
sd.TG <-rep(92.73, length(pvj))
sd.TC <- rep(42.74, length(pvj))
sd <- cbind(sd.LDL, sd.HDL, sd.TG, sd.TC)

###################################################
# Step7: SNPID and position: 186-189
###################################################
<<Step7, keep.source=TRUE, eval=FALSE>>=
hg19 <- lpd$SNP_hg19.HDL
rsid <- lpd$rsid.HDL

###################################################
# Step8: separate chromosome number and SNP position: 192-194
###################################################
chr<-chrp(hg=hg19)

###################################################
# Step9: get new data: 197-201
###################################################

newdata<-cbind(freq,beta,sd,pvj,ss,pcj)
newdata<-cbind(chr,rsid,alle1,as.data.frame(newdata))
dim(newdata)


###################################################
# Step10: retrieve data from cad and calculate pdj:204-215
###################################################

hg18.d <- cad$chr_pos_b36
SNP.d <- cad$SNP #SNPID
a1.d<- tolower(cad$reference_allele)
freq.d <- cad$ref_allele_frequency
pvalue.d <- cad$pvalue
beta.d <- cad$log_odds
N.case <- cad$N_case
N.ctr <- cad$N_control
N.d <- N.case+N.ctr
freq.case <- N.case/N.d

###################################################
 # Step11: combine these cad variables into new data sheet:218-222
 ###################################################

newcad <- cbind(freq.d, beta.d, N.case, N.ctr, freq.case)
newcad <- cbind(hg18.d, SNP.d, a1.d, as.data.frame(newcad))
dim(newcad)

################################################### 
#Step12: give name vector of causal variables: 225-227
###################################################
varname <-c("CAD", "LDL", "HDL", "TG", "TC")

###################################################
### Step 13: choose SNPs and make beta table
###################################################
mybeta <- mktable(cdata=newdata, ddata=newcad, rt="beta", 
varname=varname, LG=1, Pv=0.00000005, Pc=0.979, Pd=0.979)
dim(mybeta)
beta <- mybeta[,4:8]   #  standard beta table for path analysis
snp <- mybeta[,1:3]   #  snp data for annotation analysis
beta<-DataFrame(beta)
head(beta)

###################################################
### load beta data: 256-264
###################################################
data(beta.data)
beta.data<-DataFrame(beta.data)
CAD <- beta.data$cad
LDL <- beta.data$ldl
HDL <- beta.data$hdl
TG <- beta.data$tg
TC <- beta.data$tc

###################################################
### two way scatter: 266-280
###################################################
par(mfrow=c(2, 2), mar=c(5.1, 4.1, 4.1, 2.1), oma=c(0, 0, 0, 0))
plot(LDL,CAD, pch=19, col="blue", xlab="beta of SNPs on LDL", 
ylab="beta of SNP on CAD", cex.lab=1.5, cex.axis=1.5, cex.main=2)
abline(lm(CAD~LDL), col="red", lwd=2)
plot(HDL, CAD, pch=19,col="blue", xlab="beta of SNPs on HDL", 
ylab="beta of SNP on CAD", cex.lab=1.5, cex.axis=1.5, cex.main=2)
abline(lm(CAD~HDL), col="red", lwd=2)
plot(TG, CAD, pch=19, col="blue", xlab="beta of SNPs on TG", 
ylab="beta of SNP on CAD",cex.lab=1.5, cex.axis=1.5, cex.main=2)
abline(lm(CAD~TG), col="red", lwd=2)
plot(TC,CAD, pch=19, col="blue", xlab="beta of SNPs on TC", 
ylab="beta of SNP on CAD", cex.lab=1.5, cex.axis=1.5, cex.main=2)
abline(lm(CAD~TC), col="red", lwd=2)

###################################################
### MR and Path Analysis: 296-300
###################################################
data(beta.data)
mybeta <- DataFrame(beta.data)
mod <- CAD~LDL+HDL+TG+TC
pathvalue <- path(betav=mybeta, model=mod, outcome="CAD")

###################################################
### Create Path Diagram: 305-312
###################################################

mypath <- matrix(NA,3,4)
mypath[1,] <- c(1.000000, -0.066678, 0.420036, 0.764638)
mypath[2,] <- c(-0.066678, 1.000000, -0.559718, 0.496831)
mypath[3,] <- c(0.420036, -0.559718, 1.000000, 0.414346)
colnames(mypath) <- c("LDL", "HDL", "TG", "path")
mypath<-as.data.frame(mypath)
mypath

###################################################
### load package: diagram: 324-326
###################################################
library(diagram)

###################################################
# make path diagram
###################################################
pathdiagram(pathdata=mypath, disease="CAD", R2=0.988243, range=c(1:3))

pathD<-matrix(NA,4,5)
pathD[1,] <- c(1,	-0.070161, 0.399038, 0.907127, 1.210474)
pathD[2,] <- c(-0.070161,	1, -0.552106, 0.212201, 0.147933)
pathD[3,] <- c(0.399038, -0.552106, 1, 0.44100, 0.64229)
pathD[4,] <- c(0.907127, 0.212201, 0.441007, 1, -1.035677)
colnames(pathD) <- c("LDL", "HDL", "TG", "TC", "path")
pathD<-as.data.frame(pathD)
pathD

pathdiagram2(pathD=pathD,pathO=mypath,rangeD=c(1:4),rangeO=c(1:3),
disease="CAD", R2D=0.536535,R2O=0.988243)

###################################################
# load SNP358 data
###################################################

data(SNP358.data)
SNP358 <- DataFrame(SNP358.data)
head(SNP358)

###################################################
#  load package library(graphics)
###################################################

library(graphics)

###################################################
# chromosome SNP position histogram
##########################################
snpPositAnnot(SNPdata=SNP358,SNP_hg19="chr",main="A")


###################################################
# SNP function annotation analysis
###################################################

###################################################
# load SNP annotation data: 431-434
###################################################

data(SNP368annot.data)
SNP368<-DataFrame(SNP368annot.data)
SNP368[1:10, ]

###################################################
# load package plotrix
###################################################

###################################################
# make 3D Pie
#############################################
ucscannot(UCSCannot=SNP368,SNPn=368)

ucscannot(UCSCannot=SNP368,SNPn=368,A=3,B=2,C=1.3,method=2)


sessionInfo()


