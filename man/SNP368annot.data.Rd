\name{SNP368annot.data}
\alias{SNP368annot.data}
\docType{data}
\title{
Annotation data of \verb{368SNPs}
}
\description{
The annotation data of 368SNPs are used to construct SNP distribution in gene elements (coding region, introne,UTR, etc). The data contain 12 vectors or variables but only Symbol and function_unit are used by ucscannot.R to build SNP distribution in gene elements.
}
\usage{data("SNP368annot.data")}
\format{
  A data frame with 1053 observations on the following 6 variables.
  \describe{
    \item{\code{SNP}}{a string vector}
    \item{\code{Allele}}{a string vector} 
    \item{\code{Strand}}{a numeric vector}
    \item{\code{Symbol}}{a string vector}
    \item{\code{Gene}}{a string vector}
    \item{\code{function_unit}}{a string vector}
  }
}
\details{
\code{SNP368annot.data} were obtained by performing mktable with \code{PV=5x10e-08},\code{Pc=Pd=0.979} on lpd.data and cad.data and SNP tools. \code{SNP368annot.data} provides an practical example for constructing distribution of SNPs in gene elements. Note that \verb{function_units} are gene elements.
}
\value{
A dataset with 1053 rows and 6 columns for results of SNP annotation analysis. See format above.	
}
\source{
http://csg.sph.umich.edu//abecasis/public/lipids2013/

}
\references{
http://snp-nexus.org/test/snpnexus_19427/
}
\examples{
data(SNP368annot.data)
}
\keyword{datasets}
