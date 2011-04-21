\name{crc}
\alias{crc}
\docType{data}
\title{Real aCGH Dataset of Colorectal Cancer (CRC)}
\description{
Dataset of CRC published by Venkatachalam et al. (Identification of candidate predisposing copy number variants in familial and early-onset 
colorectal cancer patients. Int. J. Cancer, 2010). The dataset contains 30 samples that were hybridized on SNP 250k Affymetrix GeneChip arrays. Raw data are available in GEO with identifier GSE13429.
}

\author{Sandro Morganella et al.

Maintainer: S. Morganella <morganellaalx@gmail.com>
}
\usage{data(crc)}

\format{
	Data were preprocessed by PennCNV tool and discretized by Vega algorithm (Morganella et al. VEGA: variational segmentation for copy number detection. Bioinformatics, 2010) which is available as an R/Bioconductor package at the url (http://www.bioconductor.org/packages/devel/bioc/html/Vega.html). Data are organized as a matrix having a row for each observed aberrant region. Each aberrant region is described by the following columns:\cr
	Sample Name - Chromosome - Start - End - Num of Markers - CN\cr
}
\examples{
data(crc)

}
\keyword{datasets}
