\name{peel_off}
\alias{peel_off}
\title{The iterative peel-off procedure to extract the independent peak regions.}
\description{
This function implements the peel-off algorithm to extract the independent regions having the minimum q-value (lower than the given threshold) within each chromosome. The function returns for each aberration and for each chromosome the list of aberrant regions. This function uses as support the function \code{search_peaks_in_region} that extract the primary peaks.
}
\author{Sandro Morganella et al.

Maintainer: S. Morganella <morganellaalx@gmail.com>
}
\note{
This function uses the R package \code{qvalue} available at the bioconductor repository.\cr
To install the \code{qvalue} package, start R and enter:\cr
source("http://bioconductor.org/biocLite.R")\cr
biocLite("qvalue")
}
