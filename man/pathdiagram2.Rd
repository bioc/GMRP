\name{pathdiagram2}
\alias{pathdiagram2}
\title{
Two-level Nested Pathdiagram
}
\description{
This function is used to create two-level nested pathdiagram to represent causal effects of risk factors on outcome and on the disease of study. The nested path is a child path, which is related to outcome and the outside path is parent path with respect to disease.
}
\usage{
pathdiagram2(pathD, pathO, rangeD, rangeO, disease, R2D, R2O)
}
\arguments{
  \item{pathD}{
R object that is disease path result data consisting of correlation matrix of causal variables to be identified in Mendelian randomization analysis and path coefficient vector of these variables directly causing the disease of study.
}
  \item{pathO}{
R object that is outcome path result data consisting of correlation matrix of causal variables and path coefficient vector of these variables directly contributing to outcome. This outcome variable may be or not be one of risk factors or causal variables in disease path data. These variables are the same with those in \emph{pathD}.
}
  \item{rangeD}{
specifies column range for correlation coefficient matrix  in \emph{pathD}, for example, \code{rangeD=c(2:4)} means the correlation coefficient begins with column 2 and end at column 4.
}
  \item{rangeO}{
specifies column range for correlation coefficient matrix in \emph{pathO}, see example in \verb{rangeD}.
}
  \item{disease}{
a string that specifies disease name. If disease name is long or has multiple words, then we suggest an abbreviated name, for example, "coronary artery disease" can be shortened as \verb{CAD}.
}
  \item{R2D}{
a required numeric parameter and its value is R-square obtained from path analysis of disease data.
}
  \item{R2O}{
a required numeric parameter and its value is R-square obtained from path analysis of outcome data.
}
}
\details{
Two path datasets must contain correlation matrix of variables detected to be risk factor of disease and a vector of direct path coefficients obtained from path analysis of beta data of \verb{SNP}s on causal variables and disease. Columns must have shortened variable names and path word (see examples). \emph{pathdiagram2} requires two path data have the same causal variable names and the same name order. The outcome in the outcome path data must be the last variable in the correlation matrix in disease path data (see examples). Otherwise, \emph{pathdiagram2} would give an error message.
}
\value{
Null. Function \emph{pathdiagram2} creates a nested two-level path diagram labeled with color.
}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}

}

\seealso{
  \code{\link{pathdiagram}}, \code{\link[diagram]{plotmat}}, \code{\link[diagram]{plotweb}}
}
\examples{
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

#In this pathdiagram, the outcome is TC
pathD<-matrix(NA,4,5)
pathD[1,]<-c(1,	-0.070161,0.399038,0.907127,1.210474)
pathD[2,]<-c(-0.070161,	1,-0.552106,0.212201,0.147933)
pathD[3,]<-c(0.399038,-0.552106,1,0.44100,0.64229)
pathD[4,]<-c(0.907127	,0.212201,0.441007,1,-1.035677)
colnames(pathD)<-c("ldl","hdl","tg","tc","path")

#pathD
#        LDL       HDL        TG        TC      path
#1  1.000000 -0.070161  0.399038  0.907127  1.210474
#2 -0.070161  1.000000 -0.552106  0.212201  0.147933
#3  0.399038 -0.552106  1.000000  0.441000  0.642290
#4  0.907127  0.212201  0.441007  1.000000 -1.035677

pathD<-as.data.frame(pathD)
## tc is outcome in my path
pathdiagram2(pathD=pathD,pathO=mypath,rangeD=c(1:4),rangeO=c(1:3),disease="CAD",
R2D=0.536535,R2O=0.988243)

}

\keyword{path }
\keyword{diagram}
