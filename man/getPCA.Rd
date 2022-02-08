\name{getPCA}
\alias{getPCA}

\title{Principle Component Analysis (PCA) in QSea}
\description{The getPCA() function performs a
Principle Component Analysis (PCA) of the coverage profiles from a qsea object
for exploratory data analysis.
}
\usage{
getPCA(qs, chr= getChrNames(qs),ROIs, minRowSum=20, keep ,
    norm_method=normMethod(logRPM =
    c("log", "library_size", "cnv", "preference", "psC10")), topVar=1000,
    samples=getSampleNames(qs))
}
\arguments{
\item{qs}{DIPSset (mandatory)}
\item{chr}{chromosomes to consider}
\item{ROIs}{If specified, only windows overlapping ROIs are considered. }
\item{minRowSum}{minimal number of total read counts per window
        over all samples}
\item{keep}{windows to consider}
\item{norm_method}{name of predefined normalization (e.g. "beta"),
        or user defined normalization by calling normMethod() function }
\item{topVar}{only the top variable windows are considered}
\item{samples}{names of samples to be considered}
\item{minEnrichment}{
    for transformation to absolute methylation level,
    you can specify the minimal number of expected reads for a fully
    methylated window. This avoids inaccurate estimates, due to low enrichment.
}

}
\details{
The principle component analysis is calculated using the
        singular value decomposition (svd).
}
\value{
    getPCA() returns a list object, containing the svd and
    information on the selected windows.
}


\author{
Mathias Lienhard
}
\seealso{plotPCA}
\examples{
qs=getExampleQseaSet( repl=5)
pca=getPCA(qs, norm_method="beta")
colors=c(rep("red", 5), rep("green", 5))
plotPCA(pca, bgColor=colors)
#plotPCAfactors is more interesting, if ROIs have been specified in getPCA
plotPCAfactors(pca)

}

