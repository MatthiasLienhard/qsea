\name{plotCNV}
\alias{plotCNV}


\title{Plots a Heatmap-like Overview of the CNVs}
\description{
This function plots the Copy Number Variations (CNVs) 
of the samples in a heatmap like representation. 
Amplified regions are depicted in red, whereas 
deletions are depicted green, and CNV free regions blue. 
The samples are ordered by an hierarchical clustering. 
}
\usage{plotCNV(qs, dist = c("euclid", "cor")[1], clust_method = "complete",
            chr = getChrNames(qs), samples =getSampleNames(qs),
            cex = 1, labels = c(TRUE, TRUE,TRUE, TRUE), naColor = "darkgrey",
            indicateLogFC = TRUE )
}
\arguments{
\item{qs}{a qseaSet object (mandatory)}
\item{dist}{distance measure for clustering. dQuote{euclid}ian 
or dQuote{cor}relation based (1-cor)}
\item{clust_method}{method to be passed to hclust}
\item{chr}{vector of chromosomes to be depicted}
\item{samples}{samples for which CNVs are depicted}
\item{cex}{font size of labels}
\item{labels}{Boolean vector of length four (bottom, left, top, right), 
specifying the sides of the map to be labeled}
\item{naColor}{Color for regions without CNV information}
\item{indicateLogFC}{indicate the CNV logFC values in the legend}
}
\value{
This function returns the pairwise distances of the CNV profiles, 
on which the clustering is based on. 
}
\examples{
qs=getExampleQseaSet()
plotCNV(qs, labels=c(FALSE, TRUE, TRUE, FALSE))
}
\author{
Mathias Lienhard
}

