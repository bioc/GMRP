\name{fmerge}
\alias{fmerge}
\title{
Merge Two \bold{GWAS} Result Data Sheets
}
\description{
\emph{fmerg} can be used to merge two \bold{GWAS} result data sheets with the same key \verb{ID(SNP ID)} into one data sheet.
}
\usage{
fmerge(fl1, fl2,ID1, ID2, A, B, method)
}
\arguments{
  \item{fl1}{
R object: data file 1
}
  \item{fl2}{
R object: data file 2
}
  \item{ID1}{
key id (SNP ID such as rsid) in file 1
}
  \item{ID2}{
key id (SNP ID such as rsid) in file 2
}
  \item{A}{
postfix for file 1: A=".W1". W1 may be any identifier in file 1. Default is A="".
}
  \item{B}{
postfix for file 2: B=".W2". W2 may be any identifier in file2. Default is B="".
}
  \item{method}{
method for merging. See details.
}
}
\details{
fl1 and fl2 are two \bold{GWAS} result data files from different studies or with different risk variables. They contain \verb{SNPID}, \verb{hg18}, \verb{hg19}(positions), beta values, allele, frequency, and so on. The method has four options: method="No","NO" or "no" means that all data with unmatch \verb{SNP}s are not saved in the merged file; method="All","ALL" or "all" lets fmerge save all the data with unmatched \verb{SNP}s from two files but they are not paired one-by-one. This is
different from R \emph{merge} function. method="file1" will save the data with unmatched \verb{SNP}s only from file 1 in the merged file and method="file2" allows function \emph{fmerge} to save the data with unmatched \verb{SNP}s from file2 in the merged file.

}
\value{
Return a joined data sheet.
}

\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}
\note{
Function fmerg can also be applied to the other types of data.
}

\seealso{
 \code{\link[base]{merge}}
}
\examples{
data1<-matrix(NA,20,4)
data2<-matrix(NA,30,7)
SNPID1<-paste("rs",seq(1:20),sep="")
SNPID2<-paste("rs",seq(1:30),sep="")
data1[,1:4]<-c(round(runif(20),4),round(runif(20),4),round(runif(20),4),round(runif(20),4))
data2[,1:4]<-c(round(runif(30),4),round(runif(30),4),round(runif(30),4),round(runif(30),4))
data2[,5:7]<-c(round(seq(30)*runif(30),4),round(seq(30)*runif(30),4),seq(30))
data1<-cbind(SNPID1,as.data.frame(data1))
data2<-cbind(SNPID2,as.data.frame(data2))
dim(data1)
dim(data2)
colnames(data1)<-c("SNP","var1","var2","var3","var4")
colnames(data2)<-c("SNP","var1","var2","var3","var4","V1","V2","V3")
data12<-fmerge(fl1=data1,fl2=data2,ID1="SNP",ID2="SNP",A=".dat1",B=".dat2",method="No")
#data12[1:3,]
#  SNP.dat1 var1.dat1 var2.dat1 var3.dat1 var4.dat1 SNP.dat2 var1.dat2 var2.dat2
#1      rs1    0.9152    0.9853    0.9879    0.9677      rs1    0.5041    0.5734
#2     rs10    0.3357     0.116    0.3408    0.1867     rs10    0.9147    0.9294
#3     rs11    0.8004    0.8856    0.2236    0.4642     rs11    0.9262    0.5831
#  var3.dat2 var4.dat2     V1     V2 V3
#1    0.4933    0.6766 0.1864 0.6836  1
#2    0.4104    0.1327 3.2192 1.4166 10
#3    0.8541    0.6228 1.1803 1.9044 11

}
\keyword{data}
\keyword{merge}
