\name{write_significant_regions}
\alias{write_significant_regions}
\title{This function writes into a tab delimited file the significant aberrant regions.}
\description{
 This function writes into a tab delimited file the significant aberrant regions having the following format:\cr
 Chromosome - Aberration Kind - Region Start [bp] - Region End [bp] - Region Size [bp] - q-value\cr
	Chromosome: The chromosome where the aberration is found\cr
	Aberration Kind: The kind of the aberration found for the region\cr
	Region Start [bp]: The bp starting point of the aberrant region\cr
	Region End [bp]: The bp ending point of the aberrant region\cr
	Region Size [bp]: The size of the region in terms of bp\cr
	q-value: The q-value assessed for the aberrant region\cr}

\value{
	This function return a matrix containing all significant aberrant regions.
}
\author{Sandro Morganella et al.

Maintainer: S. Morganella <morganellaalx@gmail.com>
}
