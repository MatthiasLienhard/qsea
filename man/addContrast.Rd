\name{addContrast}
\alias{addContrast}

\title{fit GLMs to reduced model and test for significance}
\description{This function fits negative binomial GLMs to reduced models 
defined either by the "contrast" parameter, or by one or several model 
coefficients (specified by "coef" parameter) set to zero.
Subsequently, a likelihood ratio test is applied, to identify windows 
significantly dependent on the tested coefficient.
}
\usage{
addContrast(qs,glm,contrast,coef,name,verbose=TRUE )

}
\arguments{
\item{qs}{a qseaSet object}
\item{glm}{a qseaGLM object}
\item{contrast}{numeric vector specifying a contrast of the model coefficients. 
This contrast can for example be defined using limma::makeContrasts()}
\item{coef}{alternatively defines the contrast by coefficient(s) of the model 
tested to be equal to zero.}
\item{name}{short descriptive name for the contrast (as "TvsN"), used for 
examples in columns of result tables}
\item{verbose}{more messages that document the process}
}
\value{This function returns the qseaGLM object, extended by the fitted 
coefficients of the reduced GLMs, as well as the test statistics. 
Note that one qseaGLM object can contain several contrasts.}

\author{
Mathias Lienhard
}

\seealso{limma::makeContrasts(), fitNBglm(), isSignificant()}

\examples{
qs=getExampleQseaSet()
design=model.matrix(~group, getSampleTable(qs))
TvN_glm=fitNBglm(qs, design, norm_method="beta")
TvN_glm=addContrast(qs,TvN_glm, coef=2, name="TvN")
}

