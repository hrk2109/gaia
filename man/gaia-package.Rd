\name{gaia}
\alias{gaia}
\alias{gaia}
\docType{package}
\title{
GAIA: An R package for genomic analysis of significant chromosomal aberrations.
}
\description{
GAIA (Genomic Analysis of Important Aberrations) allows to assess the statistical significance of chromosomal aberrations. A permutation test is used to compute the probability distribution of the normal case (no significant aberrations are present in the data) so that we can estimate the statistical significance of the observed data. In order to correct for multiple hypothesis testing the False Discovery Rate approach proposed by Storey et al. (2004) is used. Finally an iterative "peel-off" procedure is used to identify the most significant independent regions.

GAIA is described in Morganella et al. (2011).

}
\details{
\tabular{ll}{
Package: \tab gaia\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2010-09-13\cr
License: \tab GNU GPL\cr
LazyLoad: \tab yes\cr
}

}
\author{Sandro Morganella et al.

Maintainer: S. Morganella <morganellaalx@gmail.com>
}
\references{
Morganella S. et al. (2011). Finding recurrent copy number alterations preserving within-sample homogeneity. Bioinformatics. DOI: 10.1093/bioinformatics/btr488.

Storey JD. et al. (2004). Strong control, conservative point estimation, and simultaneous conservative consistency of false discovery rates: A unified approach. Journal of the Royal Statistical Society. 66:187-205.
}
\examples{
# Load the matrix containing the informations about the markers
data(synthMarkers_Matrix)

# Use the function load_markers to obtain the marker descriptor data object
markers_obj <- load_markers(synthMarkers_Matrix)

# Load the matrix containing the informations about the aberrant regions
data(synthCNV_Matrix)

# Use the function load_cnv to obtain the aberrant region descriptor data object
cnv_obj <- load_cnv(synthCNV_Matrix, markers_obj, 10)

# run GAIA algorithm and save the results within the file "results.txt"
runGAIA(cnv_obj, markers_obj, "results.txt")
}
