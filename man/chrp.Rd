\name{chrp}
\alias{chrp}
\title{
Split hg19
}
\description{
Split string hg19 into two numeric columns: chr and posit.
}
\usage{
chrp(hg)
}

\arguments{
  \item{hg}{
character vector represented with code{chr##.#######} where chr## is chromosome number and \code{.#######} is SNP site (bp).
}
}
\details{
chrp can convert \code{chr##.#######} into two numeric columns: chr(chromosome number) and posit(SNP position)
}
\value{
Return two numeric vectors: chrosomome number and SNP position
}
\note{
If there is chrX.######### in the data sheet, then user should change chrX.######## to chr23.#######.	
}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}

\\Dajiang Liu
}
\note{
hg may also be hg18. User can also use packages GenomicRanges to retrieve chromosome # and SNP position.
}

\seealso{
 \code{\link{mktable}},\code{link{GenomicRanges}[Granges]},\code{link{GenomicRanges}[IRangs]},\code{link{GenomicRanges}[DataFrame]}
}
\examples{
data(lpd.data)
lpd<-lpd.data
hg19<-lpd$SNP_hg19.HDL
chr<-chrp(hg=hg19)
}
\keyword{chromosome}
\keyword{SNP}
