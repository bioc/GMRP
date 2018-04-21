\name{snpPositAnnot}
\alias{snpPositAnnot}
\title{
\emph{SNP} Position Annotation
}
\description{
This function is used to perform position annotation analysis of \verb{SNP}s chosen from \verb{GWAS}.
}
\usage{
snpPositAnnot(SNPdata, SNP_hg19="chr", main)
}
\arguments{
  \item{SNPdata}{
\verb{SNPdata} may be hg19 that is a string vector(chr##.########) or two numeric vectors (chromosome number and SNP position).
}
  \item{SNP_hg19}{
a string parameter. It may be "hg19" or "chr". If SNP_hg19="hg19",then \verb{SNPdata} contain a string vector of hg19 or if \verb{SNP_hg19}="chr", then \verb{SNPdata} consist of at lest two columns: \verb{chr} and \verb{posit}. \verb{chr} is chromosome number and posit is \verb{SNP} physical position on chromosomes. If data sheet has chromosome X, then character "X" should be changed to 23 in chr vector or chr23.######## in hg19 vector.
}
   \item{main}{
a string which is title of graph. If no title is given, then main="".
  }
}
\value{
Return a set of numbers of SNPs between which interval length > \bold{LG} on 23 chromosomes. This function also creates a histogram for averaged distances between SNPs and SNP numbers on chromosomes.
}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}
\note{
This function can also be applied to hg18 data with \code{SNP_hg19}="hg18". 
}

\seealso{
 \code{\link[graphics]{barplot}}, \code{\link[graphics]{text}}, \code{\link{chrp}}
}
\examples{
data(SNP358.data)
SNP358<-DataFrame(SNP358.data)
snpPositAnnot(SNPdata=SNP358,SNP_hg19="chr",main="A")
}
\keyword{graphics}
\keyword{SNP position}
