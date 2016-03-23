\name{path}
\alias{path}
\title{
Path Analysis
}
\description{
\emph{path} is used to perform path analysis of multiple causal or risk variables on an outcome or disease of study.
}
\usage{
path(betav,model,outcome)
}
\arguments{
   \item{betav}{
a matrix numeric data with p rows and q columns in which the first column must be outcome and other columns are risk variables.
}
   \item{model}{
a linear model for multivariate linear regression analysis. The model must be given in \code{R Console}. For example, \code{mymodel<-CAD~LDL+HDL+TG+TC}.
}
  \item{outcome}{
a string that is required to give outcome name or disease name. For example, \code{outcome="CAD"}.
}

}
\details{
\emph{path} is originally planed to perform causal analysis of risk variables on disease of study based on the results of the Mendelian randomization analysis of \verb{SNP}s on these risk variables and disease. In the \bold{GMRP} package, the betav is a matrix of beta coefficients of linear regression analyses of chosen \verb{SNP}s on the risk (or causal) variables and disease or outcome. The beta values are equivalently quantitative values, so this path function can also be used to analyze direct and indirect contributions of quantitative traits to economic traits.
}
\value{
Return three matrices: beta coefficients of regressions of risk variables on outcome or disease, correlation matrix and path matrix and also return director path coefficients and R-square.
}
\references{
Do, R. et al. 2013 Common variants associated with plasma triglycerides and risk for coronary artery disease. \emph{Nat Genet} \bold{45}: 1345-1352.

\\Sheehan, N.A. et al. 2008 Mendelian randomisation and causal inference in observational epidemiology. PLoS Med 5: e177.

\\Sheehan, N.A.,et al. 2010 Mendelian randomisation: a tool for assessing causality in observational epidemiology. \emph{Methods Mol Biol} \bold{713}: 153-166.

\\Wright, S. 1921 Correlation and causation. \emph{J.Agricultural Research} \bold{20}: 557-585.

\\Wright, S. 1934 The method of path coefficients. \emph{Annals of Mathematical Statistics} \bold{5} (3): 161-215.
}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}

\\Dajiang Liu
}
\note{
betav may also be a matrix of SNP path coefficients onto risk variables and disease.
}
\seealso{
\code{link[base]{lm}},\code{link[stats]{cor}}	
}
\examples{
data(beta.data)	
mybeta<-DataFrame(beta.data)
#mybeta<-as.data.frame(beta.data)
mod<-cad~ldl+hdl+tg+tc
pathvalue<-path(betav=mybeta,model=mod,outcome="cad")

}
\keyword{path}
\keyword{structural equation model}
