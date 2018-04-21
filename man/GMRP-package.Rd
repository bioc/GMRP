\name{GMRP-package}
\alias{GMRP-package}
\alias{GMRP}
\docType{package}
\title{
\bold{GWAS}-based Mendelian Randomization Path Analysis
}
\description{
\bold{GMRP} is used to perform Mendelian randomization analysis of causal variables on disease of study using SNP beta data from \bold{GWAS} or \bold{GWAS} meta analysis and furthermore execute path analysis of causal variables onto the disease.
}
\details{
\bold{GMRP} can perform analyses of Mendelian randomization (MR),correlation, path of causal variables onto disease of interest and SNP annotation summarization analysis. MR includes SNP selection with given criteria and regression analysis of causal variables on the disease to generate beta values of causal variables on the disease. Using the beta values, \bold{GMRP} performs correlation and path analyses to construct  path diagrams of causal variables to the disease. \bold{GMRP} consists
of 8 functions: \emph{chrp, fmerge, mktable,  pathdiagram2, pathdiagram, path, snpPositAnnot, ucscannot} and 5 datasets: \code{beta.data, cad.data, lpd.data, SNP358.data and SNP368_annot.data}. Function \emph{chrp} is used to separate string vector hg19 into two numeric vectors: chromosome number and SNP position on chromosomes. Function \emph{fmerge} is used to merge two datasets into one dataset. Function \emph{mktable} performs SNP selection and creates a standard beta table for
function path to do path analyses. Function \emph{pathdiagram} is used to create a path diagram of causal variables onto the disease or onto outcome. Function \emph{pathdiagram2} can merge two-level pathdiagrams into one nested pathdiagram where inner path diagram is a path diagram of causal variables contributing onto outcome and the outside path diagram is a diagram of path of causal variables including outcome onto the disease. The five datasets provide examples for running these
functions. \code{lpd.data} and \code{cad.data} provide an example to create a standard beta dataset for \emph{path} function to do path analysis and SNP data for SNP annotation analysis by performing \emph{mktable} and \emph{fmerge}. beta.data are a standard beta dataset for path analysis. \code{SNP358.data} provide an example for function \emph{snpPositAnnot} to do SNP position annotation analysis and \code{SNP368_annot.data} are for function \emph{ucscannot} to do SNP function annotation analysis. It is specially emphasized that except for that making standard beta table using \emph{mktable} must be done in Unix/Linux system,GMRP can be performed in Windows or Mac OS. This is because GWAS datasets usually are very huge but standard beta table is small. If users'Unix/Linux system has X11 or the other graphics system, then user should perform GMRP in Unix/Linux system, otherwise, user should transfer a standard beta table to a local computer and run GMRP in it.      

}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
\cr
Maintainer: Yuan-De Tan
}
\references{
Do, R. et al. 2013. Common variants associated with plasma triglycerides and risk for coronary artery disease. Nat Genet 45: 1345-1352.
\cr
Sheehan, N.A. et al. 2008. Mendelian randomisation and causal inference in observational epidemiology. PLoS Med 5: e177.
\cr 
Sheehan, N.A.,et al. 2010. Mendelian randomisation: a tool for assessing causality in observational epidemiology. Methods Mol Biol 713: 153-166.\\Wright, S. 1921. Correlation and causation. J. Agricultural Research 20: 557-585.
\cr
Wright, S. 1934. The method of path coefficients Annals of Mathematical Statistics 5 (3): 161-215 5: 161-215.

}
\keyword{ package }
\seealso{
  \code{\link{path}}, \code{\link{mktable}}, \code{\link{pathdiagram}},\code{\link{pathdiagram2}}, \code{\link[diagram]{plotmat}}, \code{\link[diagram]{plotweb}}
  }

\examples{

data(beta.data)
mybeta<-DataFrame(beta.data)
CAD<-beta.data$cad
LDL<-beta.data$ldl
HDL<-beta.data$hdl
TG<-beta.data$tg
TC<-beta.data$tc
#par(mfrow=c(2,2))
plot(LDL,CAD,pch=19,col="blue",xlab="beta of SNPs on LDL",ylab="beta of SNP on CAD", 
 main="A",cex.lab=1.5,cex.axis=1.5,cex.main=2)
abline(lm(CAD~LDL),col="red",lwd=2)
plot(HDL,CAD,pch=19,col="darkgreen",xlab="beta of SNPs on HDL",ylab="beta of SNP on 
 CAD", main="B",cex.lab=1.5,cex.axis=1.5,cex.main=2)
abline(lm(CAD~HDL),col="red",lwd=2)
plot(TG,CAD,pch=19,col=colors()[96],xlab="beta of SNPs on TG",ylab="beta of SNP on
 CAD",main="C",cex.lab=1.5,cex.axis=1.5,cex.main=2)
abline(lm(CAD~TG),col="red",lwd=2)
plot(TC,CAD,pch=19,col=colors()[123],xlab="beta of SNPs on TC",ylab="beta of SNP on
 CAD",main="D",cex.lab=1.5,cex.axis=1.5,cex.main=2)
abline(lm(CAD~TC),col="red",lwd=2)

mod<-cad~ldl+hdl+tg+tc
pathvalue<-path(betav=mybeta,model=mod,outcome="cad")

mypath<-matrix(NA,3,4)
mypath[1,]<-c(1.000000,-0.066678, 0.420036,0.764638)
mypath[2,]<-c(-0.066678,1.000000,-0.559718,0.496831)
mypath[3,]<-c(0.420036,-0.559718,1.000000,0.414346)
colnames(mypath)<-c("ldl","hdl","tg","path")
mypath<-DataFrame(mypath)
#mypath
#DataFrame with 3 rows and 4 columns
#        ldl       hdl        tg      path
#  <numeric> <numeric> <numeric> <numeric>
#1  1.000000 -0.066678  0.420036  0.764638
#2 -0.066678  1.000000 -0.559718  0.496831
#3  0.420036 -0.559718  1.000000  0.414346

#> pathdiagram(pathdata=mypath,disease="cad",R2=0.988243,range=c(1:3))
#Loading required package: shape
#Error in pathcad$path : $ operator is invalid for atomic vectors
mypath<-as.data.frame(mypath)
pathdiagram(pathdata=mypath,disease="cad",R2=0.988243,range=c(1:3))

}
