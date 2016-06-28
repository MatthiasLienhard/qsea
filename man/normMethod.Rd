\name{normMethod}
\alias{normMethod}

\title{Definition of normalization procedure}
\description{This function allows to define normalization 
methods by specifying components.}
\usage{
normMethod(methods, \dots)
}
\arguments{
\item{methods}{names of predefined normalization methods 
( for a list of predefined methods, see details)}
\item{\dots}{sets of normalization components, that can be 
combined to user defined normalization methods}
}
\value{a list object, containing the components for the 
specified normalization procedure}
\details{
Predefined normalization methods:

\dQuote{counts}: no normalization, simply raw count values

\dQuote{reads}: same as counts

\dQuote{rpm}: reads per million mappable reads

\dQuote{nrpm}: CNV normalized reads per million mappable reads

\dQuote{beta}: transformation to \% methylation, posterior mean point estimator

\dQuote{logitbeta}: logit transformed beta values

\dQuote{betaLB}: 2.5% quantile of beta posterior distribution, 
    lower bound for the point estimator

\dQuote{betaUB}: 97.5% quantile of beta posterior distribution, 
    upper bound for the point estimator



Allowd components for user defined normalization methods:

\dQuote{library_size}: scale by effective library size

\dQuote{region_length}: scale by window size

\dQuote{preference}: scale by positional sequencing preference

\dQuote{cnv}: scale by CNV ratio

\dQuote{enrichment}: use enrichment profiles for transformation 
to absolute methylation level

\dQuote{qXY}: quantile estimator for transformation to absolute methylation 
level. XY must be replaced by the quantile 
(see example with self defined lower and upper bound)

\dQuote{offset}: consider background reads


WARNING: not all combinations are allowed (eg qXY requires enrichment) 
and not all allowed combinations are meaningful. 
Inexperienced users should stick to predefined normalization methods.
}
\author{
Mathias Lienhard
}

\seealso{makeTable}

\examples{
#simply raw counts
nm=normMethod("counts")
#beta-values (% methylation) including lower and upper bounds
nm=normMethod(c("beta", "betaLB", "betaUB"))
#self defined lower and upper bound: 10\% and 90\% quantile
nm=normMethod("beta", 
betaLB_10=c("enrichment", "cnv", "library_size",
    "region_length", "preference","q10", "offset"),
betaUB_90=c("enrichment", "cnv", "library_size",
    "region_length", "preference","q90", "offset")
)
}

