\name{plotCoverage}
\alias{plotCoverage}


\title{Plots a genome-browser-like image of a region}
\description{This function plots the normalized coverage of specified samples 
in a specified region, together with annotations, 
in a genome-browser-like fashion}
\usage{plotCoverage(qs,test_results, chr, start, end, samples,samples2, 
norm_method="nrpkm", yoffset, xlab="Position", 
ylab="MeDIP seq", col="black", main, reorder="non", indicate_reorder=TRUE, 
distfun=dist, clustmethod="complete", scale=TRUE, steps=TRUE, space=0.05, 
baselines=TRUE, scale_val, scale_unit=NULL, logFC_pc=.1, cex=1, smooth_width, 
smooth_function=mean, regions, regions_lwd=1, regions_col, 
regions_offset, regions_names, regions_dash=0.1)
}
\arguments{
\item{qs}{a qseaSet object}
\item{chr}{the chromosome of the region to be depicted}
\item{start}{the start position of the region to be depicted}
\item{end}{the end position of the region to be depicted}
\item{samples}{the indices of the samples to be depicted}
\item{samples2}{if specified, used to calculated logFC (samples/samples2) 
profiles, must be of same length as samples}
\item{logFC_pc}{if samples2 is specified and logFC are calculated, 
this parameter specifies the pseudocount to avoid division by zero}
\item{norm_method}{a vector of normalization methods to be combined}
\item{yoffset}{horizontal offset, used to adjust the space between the profiles}
\item{xlab}{title for the x axis}
\item{ylab}{title for the y axis}
\item{main}{an overall title for the plot}
\item{col}{color vector for the samples (is recycled)}
\item{reorder}{indicate whether, and if yes how, the samples are reordered. 
Valid values are "non", "clust", "max", "minP", or a genomic position 
within the range that is depicted }
\item{test_results}{a qseaGLM object, used to find the region with minimal 
p value (only if reorder="minP")}
\item{indicate_reorder}{indicate the window that has been used 
for reordering by an arrow.}
\item{distfun}{if reorder="clust": for hierarchical clustering for reordering}
\item{clustmethod}{if reorder="clust": 
for hierarchical clustering for reordering}
\item{scale}{if set TRUE, print a bar scale}
\item{scale_val}{length of the bar scale}
\item{scale_unit}{unit of the bar scale}
\item{steps}{plot the coverage as step function (steps=TRUE), or as lines}
\item{space}{fraction of the plot set aside for sample names etc.}
\item{baselines}{depict the baselines (zero) of the coverage profiles}
\item{cex}{font size}
\item{smooth_width}{number of windows to be considered 
for sliding window smoothing}
\item{smooth_function}{function to be applied on the sliding 
windows for smoothing}
\item{regions}{named list of GenomicRanges objects, containing 
annotation (eg exons) to be depicted below the coverage profiles}
\item{regions_lwd}{vector of line width for the }
\item{regions_col}{vector of colors for the regions}
\item{regions_offset}{offset value, defining the space between the regions}
\item{regions_names}{vector of column names, 
that store the names of the regions}
\item{regions_dash}{vector, specifying the length of the end 
dashes of the regions}

}

\value{list containing a table containing the plotted coverage values, 
the position that has been used for ordering, and the image coordinates}


\author{
Mathias Lienhard
}
\examples{
qs=getExampleQseaSet(repl=5)
colors=c(rep("red", 5), rep("green", 5))
plot(1)
plotCoverage(qs,samples=getSampleNames(qs), 
    chr="chr1", start=1960001, end=1970001,col=colors, 
    norm_method="beta", yoffset=1,space=.2, reorder=1964500)
plotCoverage(qs,samples=getSampleNames(qs), 
    chr="chr1", start=1960001, end=1970001,col=colors, 
    norm_method="beta", yoffset=1,space=.2, reorder="clust")

}

