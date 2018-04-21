\name{mktable}
\alias{mktable}
\title{
Selection of SNPs and Creation of A Standard Table for Mendelian Randomization and Path Analyses
}
\description{
mktable is used to choose SNPs with \code{LG, Pv, Pc} and \code{Pd} and create a standard \verb{SNP} beta table for \code{Mendelian randomization} and \code{path analysis}, see details.
}
\usage{
mktable(cdata, ddata,rt, varname, LG, Pv, Pc, Pd)
}
\arguments{
  \item{cdata}{
causal variable \bold{GWAS} data or \bold{GWAS} meta-analysed data containing \verb{SNP} \code{ID}, \verb{SNP} position, chromosome, allele, allelic frequency, beta value, \code{sd}, sample size, etc.
}
  \item{ddata}{
disease \bold{GWAS} data or \bold{GWAS} meta-analysed data containing \verb{SNP} \code{ID}, \verb{SNP} position, chromosome, allele, allelic frequency, beta value, \code{sd}, sample size, etc.
}
  \item{rt}{
a string that specifies type of returning table. It has two options: \code{rt="beta"} returns beta table or \code{rt="path"} returns \verb{SNP} direct path coefficient table. Default is "beta".
}
  \item{varname}{
a required string set that lists names of undefined causal variables for Mendelian randomization and path analyses. The first name is disease name. Here an example given is \code{varname <-c("CAD","LDL","HD","TG","TC")}.
}
  \item{LG}{
a numeric parameter. \code{LG} is a given minimum interval distnce between \verb{SNP}s and used to choose \verb{SNP}s with. Default \code{LG=1}
}
  \item{Pv}{
a numeric parameter. \code{Pv} is a given maximum p-value that is used to choose \verb{SNP}s. Default Pv=5e-8
}
  \item{Pc}{
a numeric parameter. \code{Pc} is a given proportion of sample size to maximum sample size in causal variable data and used to choose \verb{SNP}s. Default \code{Pc=0.979}
}
  \item{Pd}{
a numeric parameter. Pd is a given proportion of sample size to the maximum sample size in disease data and used to choose \verb{SNP}s. Default \code{Pd =0.979}.
}
}
\details{
The standard \bold{GWAS} cdata set should have the format with following columns: chrn, posit, rsid, \code{a1.x1, a1.x2}, \dots, \code{a1.xn}, \code{freq.x1, freq.x2}, \dots, \code{freq.xn}, \code{beta.x1, beta.x2}, \dots, \code{beta.xn}, \code{sd.x1, sd.x2}, \dots, \code{sd.xn}, \code{pvj}, \code{N.x1, N.x2}, \dots, \code{N.xn}, \code{pcj}. The standard \bold{GWAS} \emph{ddata} set should have\code{hg.d}, \code{SNP.d},\code{a1.d}, \code{freq.d}, \code{beta.d}, \code{N.case},\code{N.ctr},\code{freq.case} where \code{x1, x2}, \dots, \code{xn} are causal variables. See example.
\describe{
\item{beta}{ is a numeric vector that is a column of beta values for regression of SNPs on variable vector \code{X={x1, x2, \dots, xn}}.
}	
\item{freq}{ is a numeric vector that is a column of frequencies of allele 1 with respect to variable vector \code{X={x1, x2, \dots, xn}}.
}
\item{sd}{is a numeric vector that is a column of standard deviations of variable \code{x1,x2}, \dots, \code{xn} specific to \verb{SNP}. Note that here sd is not beta standard deviation. If sd is not specifical to \verb{SNP}s, then sd.xi has the same value for all SNPs in variable \code{i}.
	}
\item{d}{denotes disease.}
\item{N}{is sample size.}
\item{freq.case}{is frequency of disease.}
\item{chrn}{is a numeric vector for chromosome #.}
\item{posit}{is a numeric vector for \verb{SNP} positions on chromosome #. Some time, \verb{chrn} and posit are combined into string vector: \code{hg19/hg18}.
	}
\item{pvj}{is defined as p-value, \code{pcj} and \code{pdj} as proportions of sample size for \verb{SNP} \code{j} to the maximum sample size in the causal variable data and in disease data, respectively.
  }
 }
}
\value{
Return a standard \verb{SNP} beta or \verb{SNP} path table containing \code{m} \verb{SNP}s chosen with \code{LG, Pv, Pc and Pd} and \code{n} variables and disease for \code{Mendelian randomization} and \code{path analysis}.
}
\references{
Do, R. et al. 2013. Common variants associated with plasma triglycerides and risk for coronary artery disease. \emph{Nat Genet} \bold{45}: 1345-1352.
\cr
Sheehan, N.A. et al. 2008. Mendelian randomisation and causal inference in observational epidemiology. \emph{PLoS Med} \bold{5}: e177.
\cr
Sheehan, N.A.,et al. 2010. Mendelian randomisation: a tool for assessing causality in observational epidemiology. \emph{Methods Mol Biol} \bold{713}: 153-166.
\cr
Willer, C.J. Schmidt, E.M. Sengupta, S. Peloso, G.M. Gustafsson, S. Kanoni, S. Ganna, A. Chen, J.,Buchkovich, M.L. Mora, S. et al (2013) Discovery and refinement of loci associated with lipid levels. \emph{Nat Genet} \bold{45}: 1274-1283.
}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}
\note{
The order of column variables must be \verb{chrn} \verb{posit} \verb{rsid} \code{a1.x1} \dots \code{a1.xn} \code{freq.x1} \dots \code{freq.xn} \code{beta.x1} \dots \code{beta.x1} \dots \code{beta.xn} \code{sd.x1} \dots \code{sd.xn} \dots otherwise, mktable would have error. see example.
}

\seealso{
     \code{\link{path}}
}
\examples{
data(lpd.data)
#lpd<-DataFrame(lpd.data)
lpd<-lpd.data
data(cad.data)
#cad<-DataFrame(cad.data)
cad<-cad.data
# step 1: calculate pvj
pvalue.LDL<-lpd$P.value.LDL
pvalue.HDL<-lpd$P.value.HDL
pvalue.TG<-lpd$P.value.TG
pvalue.TC<-lpd$P.value.TC
pv<-cbind(pvalue.LDL,pvalue.HDL,pvalue.TG,pvalue.TC)
pvj<-apply(pv,1,min)

#step 2: construct beta table of undefined causal variables:
beta.LDL<-lpd$beta.LDL
beta.HDL<-lpd$beta.HDL
beta.TG<-lpd$beta.TG
beta.TC<-lpd$beta.TC
beta<-cbind(beta.LDL,beta.HDL,beta.TG,beta.TC)

#step 3: construct a matrix for allele 1 in each undefined causal variable:
a1.LDL<-lpd$A1.LDL
a1.HDL<-lpd$A1.HDL
a1.TG<-lpd$A1.TG
a1.TC<-lpd$A1.TC
alle1<-cbind(a1.LDL,a1.HDL,a1.TG,a1.TC)

#step 4: calculate sample sizes of causal variables and calculate pcj
N.LDL<-lpd$N.LDL
N.HDL<-lpd$N.HDL
N.TG<-lpd$N.TG
N.TC<-lpd$N.TC
ss<-cbind(N.LDL,N.HDL,N.TG,N.TC)
sm<-apply(ss,1,sum)
pcj<-sm/max(sm)

#step 5: construct a matrix for frequency of allele1 in each undefined causal variable in 1000G.EUR
freq.LDL<-lpd$Freq.A1.1000G.EUR.LDL
freq.HDL<-lpd$Freq.A1.1000G.EUR.HDL
freq.TG<-lpd$Freq.A1.1000G.EUR.TG
freq.TC<-lpd$Freq.A1.1000G.EUR.TC
freq<-cbind(freq.LDL,freq.HDL,freq.TG,freq.TC)

#step 6: construct matrix for sd of each causal variable (here sd is not specific to SNPj)
# the sd values were averaged over 63 studies see reference Willer et al(2013) 
sd.LDL<-rep(37.42,length(pvj))
sd.HDL<-rep(14.87,length(pvj))
sd.TG<-rep(92.73,length(pvj))
sd.TC<-rep(42.74,length(pvj))
sd<-cbind(sd.LDL,sd.HDL,sd.TG,sd.TC)

#step 7: retriev SNP ID and position:
hg19<-lpd$SNP_hg19.HDL
rsid<-lpd$rsid.HDL

#step 8: invoke chrp to separate chromosome number and SNP position:
chr<-chrp(hg=hg19)

#step 9: get new data of causal variables:
newdata<-cbind(freq,beta,sd,pvj,ss,pcj)
newdata<-cbind(chr,rsid,alle1,as.data.frame(newdata))
dim(newdata)
#[1] 120165     25

#step 10: retrieve cad data from cad and calculate pdj and frequency of cad in population
hg18.d<-cad$chr_pos_b36
SNP.d<-cad$SNP #SNPID
a1.d<-tolower(cad$reference_allele)
freq.d<-cad$ref_allele_frequency
pvalue.d<-cad$pvalue
beta.d<-cad$log_odds
N.case<-cad$N_case
N.ctr<-cad$N_control
N.d<-N.case+N.ctr
freq.case<-N.case/N.d


#step 11: get new cad data:
newcad<-cbind(freq.d,beta.d,N.case,N.ctr,freq.case)
newcad<-cbind(hg18.d,SNP.d,a1.d,as.data.frame(newcad))
dim(newcad)

#step 12: give variable list
varname<-c("CAD","LDL","HDL","TG","TC")
#step 3: create beta table with function mktable 
mybeta<-mktable(cdata=newdata,ddata=newcad,rt="beta",varname=varname,LG=1, Pv=0.00000005,
Pc=0.979,Pd=0.979)

beta<-mybeta[,4:8] # save beta for path analysis
snp<-mybeta[,1:3] # save snp for annotation analysis
beta<-DataFrame(beta)
}
\keyword{Selection of SNPs}
\keyword{Mendelian Randomiztion}
