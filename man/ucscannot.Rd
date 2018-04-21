\name{ucscannot}
\alias{ucscannot}
\title{
Functional Annotation of \emph{SNP}s Chosen
}
\description{
This function is used to give proportion of \emph{SNP}s derived from functional elements of genes.
}
\usage{
ucscannot(UCSCannot,SNPn,A=3,B=1.9,C=1.3,method=1)
}
\arguments{
  \item{UCSCannot}{
annotation data obtained by performing \verb{SNP} tools.
}
\item{SNPn}{numeric parameter for number of \verb{SNP}s contained in \verb{UCSCannot}}
\item{A}{numeric parameter for title size, default=2.5}
\item{B}{numeric parameter for label size, default=1.5}
\item{C}{numeric parameter for \verb{labelrad} distance,default=0.1}
\item{method}{numeric parameter for choosing figure output methods. It has two options: method=1 has no legend but color and pie components are labeled with gene elements, method=2 has legend over pie. The default = 1.	
	 }
}
\details{
SNPs chosen by performing mktable should be copied to \bold{Batch Query Box} in \verb{SNP} annotation tool. After setting parameters and running by clicking run button, \emph{SNP} annotation data will be obtained after running for a while. Consequence sheet of \bold{UCSC} should be copied to excel sheet,"Predicted function" column name is changed to "function_unit" name and save it as \verb{csv} format. These parametric defaults are used as graph image for publication, user can expand image to the maximum size and copy it to powerpoint that will give ideal effect. User also can use R package\code{link{VariantAnnotation}} to get SNP annotation result but the result must be constructed a table with function_unit column listing gene elements and Symbol column listing genes, otherwise, ucscannot will get an error.  
}
\value{
Create a color \emph{pie3D} diagram and return a set of numeric values: proportions of code region, intron,3'and 5' UTRs and upstream and downstream etc.
}
\references{
\url{http://snp-nexus.org/index.html}
}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}

}
\note{
This function just need  data of "Predicted function" and "symbol", so the other column data in UCSCannot do not impact the results of analysis.
}
\seealso{
 \code{\link{mktable}}, \code{\link[plotrix]{pie3D}},\code{link{VariantAnnotion}}
}
\examples{
data(SNP368annot.data)
SNP368<-DataFrame(SNP368annot.data)
ucscannot(UCSCannot=SNP368,SNPn=368,A=1.5,B=1,C=1.3)
ucscannot(UCSCannot=SNP368,SNPn=368,A=1.5,B=1,C=1.3,method=2)
}
\keyword{SNP}
\keyword{annotation}
