\name{synthCNV_Matrix}
\alias{synthCNV_Matrix}
\docType{data}
\title{Example of aberration descriptor data matrix}
\description{
This data matrix simulates 10 observations (samples) for 24 chromosomes (of 1000 probes) and for 3 kinds of aberrations 0, 1 and 2 with the following aberrant regions:\cr

Chromosome 1: aberrant region of kind 0 from 301 to 700 in 100\% of samples\cr
Chromosome 2: aberrant region of kind 0 from 301 to 700 in 80\% of samples\cr
Chromosome 3: aberrant region of kind 0 from 301 to 700 in 60\% of samples\cr
Chromosome 4: aberrant region of kind 0 from 301 to 700 in 40\% of samples\cr
Chromosome 5: aberrant region of kind 0 from 301 to 700 in 20\% of samples\cr
\cr\cr
Chromosome 10: aberrant region of kind 1 from 1 to 700 in 100\% of samples\cr
Chromosome 11: aberrant region of kind 1 from 1 to 700 in 80\% of samples\cr
Chromosome 12: aberrant region of kind 1 from 1 to 700 in 60\% of samples\cr
Chromosome 13: aberrant region of kind 1 from 1 to 700 in 40\% of samples\cr
Chromosome 14: aberrant region of kind 1 from 1 to 700 in 20\% of samples\cr
\cr\cr
Chromosome 20: aberrant region of kind 2 from 801 to 1000 in 100\% of samples\cr
Chromosome 21: aberrant region of kind 2 from 801 to 1000 in 80\% of samples\cr
Chromosome 22: aberrant region of kind 2 from 801 to 1000 in 60\% of samples\cr
Chromosome 23: aberrant region of kind 2 from 801 to 1000 in 40\% of samples\cr
Chromosome 24: aberrant region of kind 2 from 801 to 1000 in 20\% of samples\cr

}

\author{Sandro Morganella et al.

Maintainer: S. Morganella <morganellaalx@gmail.com>
}
\usage{data(synthCNV_Matrix)}

\format{
	This data matrix is organized as a matrix having a row for each observed aberrant region. Each aberrant region is described by the following columns:\cr
	Sample Name - Chromosome - Start - End - Num of Markers - CN\cr
	"Sample Name" indicates the name of the sample. "Chromosome", "Start", "End", "Num of Markers" and "CN" indicate for each aberrant region the respective chromosome, the start and the end position (in bp) the number of markers contained within the region and the found aberrations. Note that "CN" represents the estimated copy number for the segmented region and it must be an integer in the range 0..(K-1) where K is the number of the considered aberrations. In this data matrix three different aberration kinds are considered: 0, 1 and 2.
}
\examples{
data(synthCNV_Matrix)

}
\keyword{datasets}
