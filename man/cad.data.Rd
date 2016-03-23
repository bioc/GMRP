\name{cad.data}
\alias{cad.data}
\docType{data}
\title{
bold{GWAS} Meta-analyzed Data of Coronary Artery Disease
}
\description{
\code{cad.data} are a matrix dataset consisting of 12 variables such as \verb{SNPID}, \verb{SNP} position on chromosomes, allele and alternative allele, allelic frequencies and 1069 \verb{SNP}s.
}
\usage{data("cad.data")}
\format{
  A data frame with 1609 observations on the following 12 variables.
  \describe{
    \item{\code{SNP}}{a character vector}
    \item{\code{chr_pos_b36}}{a character vector}
    \item{\code{reference_allele}}{a character vector}
    \item{\code{other_allele}}{a character vector}
    \item{\code{ref_allele_frequency}}{a numeric vector}
    \item{\code{pvalue}}{a numeric vector}
    \item{\code{het_pvalue}}{a numeric vector}
    \item{\code{log_odds}}{a numeric vector}
    \item{\code{log_odds_se}}{a numeric vector}
    \item{\code{N_case}}{a numeric vector}
    \item{\code{N_control}}{a numeric vector}
    \item{\code{model}}{a character vector}
  }
}
\details{
\code{cad.data}, also called \verb{CARDIoGRAM} \bold{GWAS}, are a meta-analyzed GWAS data from \bold{GWAS} studies of European descent imputed to \bold{HapMap2} involving 22,233 cases and 64,762 controls. The data were downloaded from the following website.
}
\value{
A data sheet consisting of 1609 rows(SNPs) and 12 columns(character vectors such as \verb{SNPID}  and allele, numeric vector such as allele frequency and beta coefficient. See data format above). 	
}
\source{
http://www.cardiogramplusc4d.org/downloads/
}
\references{
Schunkert H, Konig IR, Kathiresan S, Reilly MP, Assimes TL, Holm H et al. Large-scale association analysis identifies 13 new susceptibility loci for coronary artery disease. \bold{Nat Genet.} 2011 \bold{43}: 333-338.
}
\examples{
data(cad.data)
## maybe str(cad.data) ; plot(cad.data) ...
}
\keyword{datasets}
