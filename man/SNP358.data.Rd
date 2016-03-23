\name{SNP358.data}
\alias{SNP358.data}
\docType{data}
\title{
Data of 358 SNPs
}
\description{
\code{SNP358.data} were obtained from \bold{GWAS} Meta-analyzed datasets of lipoprotein cholesterols and coronary artery disease. The data contain three numeric vectors (columns): SNPID(\verb{rsid}), chromosome number (\verb{chr}) and SNP position on chromosome (posit).
}
\usage{data("SNP358.data")}
\format{
  A data frame with 358 observations on the following 3 variables.
  \describe{
    \item{\code{rsid}}{a character vector}
    \item{\code{chr}}{a numeric vector}
    \item{\code{posit}}{a numeric vector}
  }
}
\details{
These 358 SNPs were chosen by using mktable with Pv = code{5x10e-08},\code{Pc = Pd =0.979} from code{lpd.data and cad.data}. They provide a data example for how to perform annotation analysis of SNP positions.
}
\value{
A set of data with 358 rows(SNPs) and 3 columns(SNP ID, chromosome # and SNP position on chromosomes).
}

\examples{
data(SNP358.data)
## maybe str(SNP358.data) ; plot(SNP358.data) ...
}
\keyword{datasets}
