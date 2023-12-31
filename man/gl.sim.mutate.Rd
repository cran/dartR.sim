% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.sim.mutate.r
\name{gl.sim.mutate}
\alias{gl.sim.mutate}
\title{Simulates mutations within a genlight object}
\usage{
gl.sim.mutate(x, mut.rate = 1e-06)
}
\arguments{
\item{x}{Name of the genlight object containing the SNP data [required].}

\item{mut.rate}{Constant mutation rate over nInd*nLoc*2 possible locations
[default 1e-6]}
}
\value{
Returns a genlight object with the applied mutations
}
\description{
This script is intended to be used within the simulation framework of dartR. 
It adds the ability to add a constant mutation rate across all loci. Only 
works currently for biallelic data sets (SNPs). Mutation rate is checking for 
all alleles position and mutations at loci with missing values are ignored 
and in principle 'double mutations' at the same loci can occur, but should be 
rare.
}
\examples{
b2 <- gl.sim.mutate(bandicoot.gl,mut.rate=1e-4 )
#check the mutations that have occurred
table(as.matrix(bandicoot.gl), as.matrix(b2))
}
\author{
Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
}
