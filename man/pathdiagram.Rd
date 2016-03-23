\name{pathdiagram}
\alias{pathdiagram}
\title{
Path Diagram
}
\description{
Create a directed acyclic diagram to represent causal effects of risk factors on the disease of study.
}
\usage{
pathdiagram(pathdata, disease, R2, range)
}
\arguments{
  \item{pathdata}{
R object that is dataset consisting of correlation matrix of risk factors and a numeric vector of direct path coefficients.
}
  \item{disease}{
a string that specifies outcome name or disease name. If disease name is long or has multiple words, then we suggest an abbreviated name, for example, coronary artery disease may be shortened as \code{CAD}.
}
  \item{R2}{
a numeric parameter, is R-square obtained from path analysis.
}
  \item{range}{
range of specified columns for correlation matrix. For example, \code{range = c(2:4)} means the correlation coefficient begins with column 2 and end at column 4.
}
}
\details{
The \emph{pathdata} contains correlation matrix of risk factors and a numeric vector of direct path coefficients obtained from path analysis of beta data of \code{SNP}s on risk factors and disease. Columns must have risk factor names and path.
}
\value{
NULL. pathdiagram will create one-level path diagram labeled with colors.
}

\author{
Yuan-De Tan
\email{tanyuande@gmail.com}

\\Dajiang Liu
}

\seealso{
    \code{\link{pathdiagram2}}, \code{\link[diagram]{plotmat}}, \code{\link[diagram]{plotweb}}
}

\examples{

mypath<-matrix(NA,3,4)
mypath[1,]<-c(1.000000,-0.066678, 0.420036,0.764638)
mypath[2,]<-c(-0.066678,1.000000,-0.559718,0.496831)
mypath[3,]<-c(0.420036,-0.559718,1.000000,0.414346)
colnames(mypath)<-c("ldl","hdl","tg","path")

#mypath
#        ldl       hdl        tg      path
#1  1.000000 -0.066678  0.420036  0.764638
#2 -0.066678  1.000000 -0.559718  0.496831
#3  0.420036 -0.559718  1.000000  0.414346

mypath<-as.data.frame(mypath)
pathdiagram(pathdata=mypath,disease="cad",R2=0.988243,range=c(1:3))
}

\keyword{path }
\keyword{diagram}
