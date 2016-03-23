\name{lpd.data}
\alias{lpd.data}
\docType{data}
\title{
\bold{GWAS} Meta-analyzed Data of Lipoprotein Cholesterols
}
\description{
lpd.data are standard \bold{GWAS} Meta-analyzed dataset of lipoprotein cholesterols. It was constructed by merging four datasets: \code{Mc_HDL.txt,Mc_LDL.txt, Mc_TC.txt and Mc_TG.txt}.
}
\usage{data("lpd.data")}
\format{
  A data frame with 1609 observations on the following 40 variables.
  \describe{
    \item{\code{SNP_hg18.HDL}}{a character vector}
    \item{\code{SNP_hg19.HDL}}{a character vector}
    \item{\code{rsid.HDL}}{a character vector}
    \item{\code{A1.HDL}}{a character vector}
    \item{\code{A2.HDL}}{a character vector}
    \item{\code{beta.HDL}}{a numeric vector}
    \item{\code{se.HDL}}{a numeric vector}
    \item{\code{N.HDL}}{a numeric vector}
    \item{\code{P.value.HDL}}{a numeric vector}
    \item{\code{Freq.A1.1000G.EUR.HDL}}{a numeric vector}
    \item{\code{SNP_hg18.LDL}}{a character vector}
    \item{\code{SNP_hg19.LDL}}{a character vector}
    \item{\code{rsid.LDL}}{a character vector}
    \item{\code{A1.LDL}}{a character vector}
    \item{\code{A2.LDL}}{a character vector}
    \item{\code{beta.LDL}}{a numeric vector}
    \item{\code{se.LDL}}{a numeric vector}
    \item{\code{N.LDL}}{a numeric vector}
    \item{\code{P.value.LDL}}{a numeric vector}
    \item{\code{Freq.A1.1000G.EUR.LDL}}{a numeric vector}
    \item{\code{SNP_hg18.TG}}{a character vector}
    \item{\code{SNP_hg19.TG}}{a character vector}
    \item{\code{rsid.TG}}{a character vector}
    \item{\code{A1.TG}}{a character vector}
    \item{\code{A2.TG}}{a character vector}
    \item{\code{beta.TG}}{a numeric vector}
    \item{\code{se.TG}}{a numeric vector}
    \item{\code{N.TG}}{a numeric vector}
    \item{\code{P.value.TG}}{a numeric vector}
    \item{\code{Freq.A1.1000G.EUR.TG}}{a numeric vector}
    \item{\code{SNP_hg18.TC}}{a character vector}
    \item{\code{SNP_hg19.TC}}{a character vector}
    \item{\code{rsid.TC}}{a character vector}
    \item{\code{A1.TC}}{a character vector}
    \item{\code{A2.TC}}{a character vector}
    \item{\code{beta.TC}}{a numeric vector}
    \item{\code{se.TC}}{a numeric vector}
    \item{\code{N.TC}}{a numeric vector}
    \item{\code{P.value.TC}}{a numeric vector}
    \item{\code{Freq.A1.1000G.EUR.TC}}{a numeric vector}
  }
}
\details{
These four \bold{GWAS} Meta-analyzed data of lipoprotein cholesterols were downloaded from the following website.
}
\value{
A data sheet consisting of 1609 rows (SNPs) and 40 columns(character vectors such as SNPID and allele, numeric vector such as allele frequency, beta coefficient. See data format above). 		
}
\source{
http://csg.sph.umich.edu//abecasis/public/lipids2013/
}
\references{
Willer CJ et al. Discovery and refinement of loci associated with lipid levels. \emph{Nat. Genet.} 2013. doi:10.1038/ng.2797.
}
\examples{
data(lpd.data)
## maybe str(lpd.data) ; plot(lpd.data) ...
}
\keyword{datasets}
