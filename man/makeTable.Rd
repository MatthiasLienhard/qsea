\name{makeTable}
\alias{makeTable}


\title{Create a Results Table}
\description{
    This function creates a table from the qsea objects 
    qseaSet and qseaTvN_glm
}
\usage{
makeTable(qs,glm,norm_methods="counts",samples,groupMeans, keep, ROIs,
    annotation, minPvalSummarize, CNV=FALSE, verbose=TRUE, minEnrichment=3)
}
\arguments{
\item{qs}{a qseaSet object (mandatory)}
\item{glm}{a list of one or more qseaGLM objects (optional)}
\item{norm_methods}{
    ether a character vector of pre-defined normalization combinations, 
    or a list defining normalization combinations. 
    This affects both individual and mean values. 
}
\item{samples}{
    The indices of the samples for which individual values 
    are to be written out in the specified order
}
\item{groupMeans}{a named list of indices vectors, 
    defining groups for which mean values are to be written out}
\item{keep}{a vector of indices of the windows 
    that are considered (as created by isSignificant)}
\item{ROIs}{A GRanges object, containing regions of interest (ROIs). 
    Only windows overlapping ROIs are considered.}
\item{annotation}{
    a named list of GRange objects, 
    containing annotations (e.g. genes, CpG islands, ...) 
    that are added to the table. 
}
\item{minPvalSummarize}{
    If ROIs are given, you can specify a QseaTvN_glm object. 
    For each ROI the window with the most significant differential 
    coverage is written out
}
\item{CNV}{
    If set TRUE, the CNV logFC for the 
    samples specified by samples are written out. 
}
\item{verbose}{ verbosity level}
\item{minEnrichment}{
    for transformation to absolute methylation level, 
    you can specify the minimal number of expected reads for a fully 
    methylated window. This avoids inaccurate estimates, due to low enrichment.
}
}
\details{
    Note that, if overlapping ROIs are specified, windows might 
    emerge in the table several times. 
}
\value{
    A result table containing the specified normalized values 
    for the selected windows and samples/groups
}


\author{
Mathias Lienhard
}
\seealso{isSignificant}
\examples{
#create example set
qs=getExampleQseaSet()
design=model.matrix(~group, getSampleTable(qs))
TvN_glm=fitNBglm(qs, design, norm_method="beta")
TvN_glm=addContrast(qs,TvN_glm, coef=2, name="TvN")
sig=isSignificant(TvN_glm, fdr_th=0.01)

##Table containing all significant windows
tab1=makeTable(qs=qs, glm=TvN_glm, 
    keep=sig, samples=getSampleNames(qs))
##additional CNV logFC for the selected samples
tab2=makeTable(qs=qs, glm=TvN_glm, 
    keep=sig, samples=getSampleNames(qs), CNV=TRUE)
##explicit selection of normalization: 
##counts (i.e. no normalization, only counts)
tab3=makeTable(qs=qs, glm=TvN_glm, keep=sig,
    samples=getSampleNames(qs), norm_method="counts")

##counts AND %methylation values for individual samples and group means
tab4=makeTable(qs=qs, glm=TvN_glm, keep=sig,
    samples=getSampleNames(qs), groupMeans=getSampleGroups(qs),
    norm_method=c("counts", "beta"))

}

