\name{plotEnrichmentProfile}
\alias{plotEnrichmentProfile}
\alias{plotEPmatrix}

\title{Plotting functions for enrichment profiles}
\description{Plots the estimated sequence pattern dependent enrichment 
profile for one or several samples as a matrix of plots}
\usage{
plotEnrichmentProfile(qs,sample, sPoints=seq(0,30,1), 
    fitPar=list(lty=2, col="green"),cfPar=list(lty=1), densityPar, meanPar,... )
plotEPmatrix(qs, sa=getSampleNames(qs),nrow=ceiling(sqrt(length(sa))), 
    ncol=ceiling(length(sa)/nrow), ...)
}
\arguments{
\item{qs}{The qseaSet object}
\item{sample}{The index of the sample for which the 
enrichment profile should be depicted}
\item{sPoints}{The values at which the enrichment profile 
function is evaluated}
\item{fitPar}{List of parameters for depiction of the fitted 
enrichment profile function (see details)}
\item{cfPar}{List of parameters for depiction of the empirical 
enrichment profile (see details)}
\item{densityPar}{List of parameters for depiction high density 
scatterplot of coverage and pattern density (see details)}
\item{meanPar}{List of parameters for depiction of the mean 
coverage per pattern density bin (see details)}
\item{sa}{vector of samples to be depicted in matrix plot}
\item{nrow}{number of rows in matrix plot}
\item{ncol}{number of columns in matrix plot}
\item{\dots}{Further graphical parameters may also be supplied}

}
\details{
Parameter lists for lines in the plot (e.g. fitPar, cfPar and meanPar) are 
passed to graphics::lines(), densityPar are passed to graphics::smoothScatter() 
function. 
}

\value{
plotEnrichmentProfile returns the coordinates of the enrichment profile. 
plotEPmatrix returns enrichment profile coordinates for all depicted samples. 
}
\author{
Mathias Lienhard
}

\seealso{addEnrichmentParameters}
\examples{
#create example object with different sequencing depth
qs=getExampleQseaSet(expSamplingDepth=50*10^(1:4), repl=4)
#enrichment profile for one sample
plotEnrichmentProfile(qs, "Sim4T")
#enrichment profile for all samples
plotEPmatrix(qs)
}

