\name{beta.data}
\alias{beta.data}
\docType{data}
\title{
Beta Data Of SNP Regressed on Causal Variables and Disease 
}
\description{
Beta data are a matrix dataset consisting of 5 columns: cad, ldl, hdl, tg, and tc with 368 rows.
}
\usage{data("beta.data")}
\format{
  A data frame with 368 observations on the following 5 variables.
  \describe{
    \item{\code{cad}}{a numeric vector}
    \item{\code{ldl}}{a numeric vector}
    \item{\code{hdl}}{a numeric vector}
    \item{\code{tg}}{a numeric vector}
    \item{\code{tc}}{a numeric vector}
  }
}
\details{
Beta data are a matrix consisting of regression coefficients of 368 SNPs on cad, ldl, hdl, tg, tc where \verb{cad} is coronary artery disease, \verb{ldl} is low-density lipoprotein cholesterol, \verb{hdl}, high-density lipoprotein cholesterol, \verb{tg},triglycerides and \verb{tc}, total cholesterol in plasma. These 368 SNPs were obtained by using \verb{mktable} from \bold{GWAS} meta-analyzed data.
}
\value{
A set of real regression coefficients of 368 SNPs on disease and causal variables.	
}
\source{
http://csg.sph.umich.edu//abecasis/public/lipids2013/

\\http://www.cardiogramplusc4d.org/downloads/

}
\references{
Willer CJ et al. Discovery and refinement of loci associated with lipid levels. Nat. Genet. 2013. doi:10.1038/ng.2797.

\\Schunkert, H. et al. 2011. Large-scale association analysis identifies 13 new susceptibility loci for coronary artery disease. Nat Genet 43: 333-338.[online]

\\Schunkert H, Konig IR, Kathiresan S, Reilly MP, Assimes TL, Holm H et al. Large-scale association analysis identifies 13 new susceptibility loci for coronary artery disease. Nat Genet. 2011 43: 333-338.
}
\examples{
data(beta.data)
## maybe str(beta.data) ; plot(beta.data) ...
}
\keyword{datasets}
