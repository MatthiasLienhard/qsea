\name{addOffset}
\alias{addOffset}


\title{Estimate background reads}
\description{This function sets the background reads offset parameters for the 
qseaSet object, either by estimating offset reads, or by setting user provided 
values. }
\usage{
addOffset(qs,enrichmentPattern , maxPatternDensity=0.01,offset)

}
\arguments{
\item{qs}{the qseaSet object}
\item{enrichmentPattern}{name of the enrichment pattern, as specified in 
addPatternDensity}
\item{maxPatternDensity}{Maximum pattern density, at which the window is 
treated as pattern free.}
\item{offset}{This parameter alternatively allows to specify the amount of 
background reads for each sample manually. In this case, please provide 
average background reads for CNV free windows in rpkm scale.}
}
\value{The function returns the qseaSet object, extended by the estimated 
amount of background reads for all samples}

\author{
Mathias Lienhard
}

\seealso{addPatternDensity, getOffset}

\examples{
#simulate data with varing background fractions
qs=getExampleQseaSet(expSamplingDepth=5e4, repl=5,bgfraction=seq(0,.8,.2))
#estimate the background in simulated data
addOffset(qs, "CpG", maxPatternDensity=0.7)
#return the background on different scales
getOffset(qs, scale="fraction") #estimated fraction of total reads
getOffset(qs, scale="rpw") #average background reads per CNV free window

}

