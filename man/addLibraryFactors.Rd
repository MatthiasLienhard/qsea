\name{addLibraryFactors}
\alias{addLibraryFactors}


\title{Estimate effective library size}
\description{Normalization factors for effective library size are computed 
using the trimmed mean of m-values approach (TMM).}
\usage{
addLibraryFactors(qs, factors,\dots)
}
\arguments{
\item{qs}{The qseaSet object}
\item{factors}{In case normalization factors have been pre-computed by the user,
they can be passed with this parameter. 
In this case QSEA adds this factors to the qseaSet object and does not compute 
normalization factors.}
\item{\dots}{Further parameters used for the TMM normalization (see details)}
}
\details{
The user can specify the TMM normalization by setting the following 
additional parameters, which are passed to the internal functions.
\\trimA [default: c(.5,.99)] lower and upper quantiles for trimming of A values
\\trimM [default: c(.1,.9)] lower and upper quantiles for trimming of M values
\\doWeighting [default: TRUE] computes a weighted TMM
\\ref [default: 1] the index of the reference sample
\\plot [default: FALSE] if set to TRUE, MvsA plots depicting the TMM 
normalization are created. 
\\nReg [default: 500000] Number of regions to be analyzed for normalization. 
Regions are drawn uniformly over the whole genome. 
}


\value{This function returns the qseaSet object, containing effective 
library size normalization factors.}

\author{
Mathias Lienhard
}

\seealso{edgeR::calcNormFactors}
\examples{
qs=getExampleQseaSet(expSamplingDepth=500*10^(1:5), repl=5)
#in this case, the first sample has only view reads, so it is important to set
#the reference sample
qs=addLibraryFactors(qs, plot=TRUE, ref="Sim5N")

}

