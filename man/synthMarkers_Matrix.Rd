\name{synthMarkers_Matrix}
\alias{synthMarkers_Matrix}
\docType{data}
\title{Example of marker descriptor data matrix}
\description{
	This matrix simulates the marker descriptors of 24 chromosomes (of 1000 probes) where the start and the end position of each probe coincide. 
}

\usage{data(synthMarkers_Matrix)}
\author{Sandro Morganella et al.

Maintainer: S. Morganella <morganellaalx@gmail.com>
}
\format{
	This data matrix is organized as a matrix having a row for each measured probe. Each probe is described by the following columns:\cr
	Probe Name - Chromosome - Start\cr

  	Where: "Probe Name" is the name of the observed probe; "Chromosome" is the chromosome where the probe is located and "Start" is the start point (in bp) of the probe. Note that the matrix can have also a column ("End") specifying the position in which the probe ends, this column is optional and if it is missed than start and the end positions will coincide.\cr
	\code{synthMarkers_Matrix} reports the marker descriptor for 24 chromosomes of 1000 probes (where the sex chromosomes X and Y are indicated with 23 and 24 respectively.
}
\examples{
data(synthMarkers_Matrix)

}
\keyword{datasets}
