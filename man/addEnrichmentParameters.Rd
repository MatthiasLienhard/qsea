\name{addEnrichmentParameters}
\alias{addEnrichmentParameters}


\title{Enrichment analysis}
\description{This function analyses the dependency of enrichment on a 
sequence pattern, based on a subset of windows for which the signal is known.}
\usage{
addEnrichmentParameters(qs, enrichmentPattern, signal, windowIdx, 
min_wd=5,bins=seq(.5,40.5,1))
}
\arguments{
\item{qs}{The qseaSet object}
\item{enrichmentPattern}{The name of the pattern, on which the enrichment 
depends on (usually CpG for methylation analysis). 
This name must correspond to the name specified in addPatternDensity()}
\item{windowIdx}{vector of window indices, for which "true" 
values are known (or can be estimated)}
\item{signal}{Matrix containing the known (or estimated) values for 
all samples and all specified windows, as a numeric matrix. 
These values are expected to be between 0 and 1.}
\item{bins}{For the enrichment analysis, windows are binned according 
to pattern density. This parameter specifies the bins. }
\item{min_wd}{minimal number of windows per bin to be considered}
}
\value{The function returns the qseaSet object, extended by the parameters 
of the enrichment profiles for all samples}

\author{
Mathias Lienhard
}

\seealso{plotEnrichmentProfile, addPatternDensity}

\examples{
qs=getExampleQseaSet(enrichmentAnalysis=FALSE)
#this procedure assumes that regions with low CpG density is 80% methylated
#on average, and regions within CpG islands are 25% methylated on average. 
wd=which(getRegions(qs)$CpG_density>1 &
    getRegions(qs)$CpG_density<15)
signal=(15-getRegions(qs)$CpG_density[wd])*.55/15+.25
signal=matrix(signal,nrow=length(signal),ncol=length(getSampleNames(qs)))
qs=addEnrichmentParameters(qs, enrichmentPattern="CpG", 
    windowIdx=wd, signal=signal)
}

