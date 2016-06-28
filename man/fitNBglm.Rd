\name{fitNBglm}
\alias{fitNBglm}

\title{Fit GLM for each window}
\description{This function fits a negative binomial GLM for each genomic 
window, according to the design matrix.}
\usage{
fitNBglm(qs,design,link="log",keep, disp_method="region_wise",
norm_method="rpkm",init_disp=0.5 ,verbose=TRUE, minRowSum=10, pseudocount=1 )
}
\arguments{
\item{qs}{a qseaSet object}
\item{design}{the design matrix for the GLMs}
\item{link}{name of the link function. Currently, only the canonical 
dQuote{log} link function is implemented. }
\item{keep}{indices of windows to be included in the analysis.}
\item{disp_method}{method to estimate dispersion parameters. 
Allowed values are dQuote{region_wise} for independent window wise estimates, 
dQuote{common} for a single estimate for all windows, 
dQuote{cutAtQuantiles} for window wise estimates trimmed at the 25\% and 75\% 
quantiles, 
or dQuote{initial} for using the dispersion parameters provided with the 
init_disp parameter. }
\item{norm_method}{normalization method, as defined by normMethod() function}
\item{init_disp}{initial estimate for dispersion parameter. Either a single 
parameter for all regions, or a vector with window wise parameters.}
\item{verbose}{more messages that document the process}
\item{minRowSum}{filter out windows with less than minRowSum reads over all 
samples}
\item{pseudocount}{this value is added to the read counts, to avoid log(0)}
}
\value{This function returns a qseaGLM object, containing the fitted 
coefficients of the GLMs.}

\author{
Mathias Lienhard
}

\seealso{addContrast()}

\examples{
#tbd
qs=getExampleQseaSet()
design=model.matrix(~group, getSampleTable(qs))
TvN_glm=fitNBglm(qs, design, norm_method="beta")
}

