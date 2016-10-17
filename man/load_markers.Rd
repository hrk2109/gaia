\name{load_markers}
\alias{load_markers}
\title{This function create the marker descriptor object containing all needed marker informations.}
\description{
This function loads the markers contained within the matrix passed as argument and creates for all chromosomes an ordered vector containing the position of each marker.These vectors are loaded within a list.

}
\usage{
load_markers(marker_matrix)
}

\arguments{
  \item{marker_matrix}{ contains the marker descriptions as a matrix with the following structure:\cr
  Probe Name - Chromosome - Start - End\cr
  Note that the End position column is optional (in this case start and the end positions coincide). The sex chromosomes X and Y must be indicated with 23 and 24 respectively.}
}
\value{
 This function returns a list having the following structure:
 \item{chromosome_marker_list[[i]]}{is a matrix of dimension 2xN (N is the number of observed probes for the i-th chromosome) the first and the second row contains the start and the end position of each marker of the i-th chromosome respectively.\cr\cr}
An example of the data produced by this function can be found in \code{synthMarkers}.
}
%\references{}
\author{Sandro Morganella et al.

Maintainer: S. Morganella <morganellaalx@gmail.com>
}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
# Load the matrix containing the informations about the markers
data(synthMarkers_Matrix)

# Use the function load_markers to obtain the marker descriptor data object
marks <- load_markers(synthMarkers_Matrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
