\name{isSignificant}
\alias{isSignificant}

\title{Finds Significant Regions}
\description{This function looks for regions, where the test statistic 
    is below the defined thresholds}
\usage{
isSignificant(glm, contrast = NULL, fdr_th = NULL, pval_th = NULL,
            absLogFC_th = NULL, direction = "both")
}
\arguments{
\item{glm}{A qseaGLM object (mandatory)}
\item{contrast}{name of contrast to be used}
\item{fdr_th}{a threshold for the false discovery rate}
\item{pval_th}{a p value threshold}
\item{absLogFC_th}{the threshold for the absolute value of logFC}
\item{direction}{direction of change: either "both", "loss", or "gain"}
}
\details{
If a threshold is NULL, it is ignored.

For the direction parameter, the following synonyms are valid:

"loss" == "less" == "hypo"

"gain" == "more" == "hyper"

}
\value{A vector with indices of significant windows, which can be passed 
    to keep parameter of makeTable() function}


\author{
Mathias Lienhard
}
\seealso{makeTable}
\examples{
qs=getExampleQseaSet()
design=model.matrix(~group, getSampleTable(qs))
TvN_glm=fitNBglm(qs, design, norm_method="beta")
TvN_glm=addContrast(qs,TvN_glm, coef=2, name="TvN")
sig=isSignificant(TvN_glm, fdr_th=0.01)
}

