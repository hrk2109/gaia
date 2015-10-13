\name{load_cnv}
\alias{load_cnv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{This function create the object containing all needed informations about the aberrant regions.}
\description{
This function loads the informations about the aberrant regions contained within the matrix passed as argument. It creates for all chromosomes and all kind of aberration (e.g. loss and gain) a matrix of dimension NxM (N observed samples and M observed probes).
}
\usage{
load_cnv(segmentation_matrix, markers_list, num_of_samples)
}
\arguments{
  \item{segmentation_matrix}{A matrix containing the aberrant regions where each row in the file reports the information of an aberrant region. In particular the matrix has the following column:\cr
  	Sample Name - Chromosome - Start - End - Num of Markers - CN\cr
	"Sample Name" indicates the name of the sample. "Chromosome", "Start", "End", "Num of Markers" and "CN" indicate for each aberrant region the respective chromosome, the start and the end position the number of markers contained within the region and the found aberrations. Note that "CN" represents the estimated copy number for the segmented region and it must be an integer in the range 0..(K-1) where K is the number of the considered aberrations Therefore if we are considering only loss, LOH, gain in the file passed as argument the only possible kind of aberrations is 0, 1 and 2.
}
  \item{markers_list}{The marker descriptor object obtained by the function \code{load_markers}.}
  \item{num_of_samples}{The number of analyzed samples.}
}
\value{ This function returns a list having the following structure:
	\item{CNV_matrix_list[[i]][[j]]}{contains the informations for the j-th chromosome on the i-th aberration. This element is a matrix of dimension NxM (N observed samples and M observed probes).}
	An example of the data produced by this function can be found in \code{synthCNV}

}
%\references{ ~put references to the literature/web site here ~ }
\author{Sandro Morganella et al.

Maintainer: S. Morganella <morganellaalx@gmail.com>
}

%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
# Load the matrix containing the informations about the markers
data(synthMarkers_Matrix)

# Use the function load_markers to obtain the marker descriptor data object
marks <- load_markers(synthMarkers_Matrix)

# Load the matrix containing the informations about the aberrant regions
data(synthCNV_Matrix)

# Use the function load_cnv to obtain the aberrant region descriptor data object
cnv <- load_cnv(synthCNV_Matrix, marks, 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
