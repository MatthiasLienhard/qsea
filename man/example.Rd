\name{getExampleQseaSet}
\alias{getExampleQseaSet}

\title{Simulation of MeDIP seq QSEA set}
\description{Creates a example qseaSet object by sampling reads 
    for simulated Tumor and Normal samples. Number of replicates, 
    sequencing depth and fraction of background reads can be specified.}
\usage{
getExampleQseaSet(CpG=TRUE,CNV=TRUE,repl=2, 
    doSampling=TRUE,enrichmentAnalysis=TRUE, expSamplingDepth=50000, 
    bgfraction=.1)
}
\arguments{
\item{CpG}{if TRUE CpG density is added to the object}
\item{CNV}{if TRUE CNV are emulated for the tumor samples}
\item{repl}{number of replicates for tumor and normal samples}
\item{doSampling}{if TRUE, read counts are sampled and added to the object}
\item{enrichmentAnalysis}{if TRUE, parameters for enrichment profiles are added}
\item{expSamplingDepth}{expected value of sequencing depth}
\item{bgfraction}{fraction of background reads}
}
\details{
The function creates an example and test qseaSet object for an toy example 
genome (one chromosome, 50kb) with 500 bases windows. 
}
\value{
The qseaSet object
}


\author{
Mathias Lienhard
}

\seealso{createQseaSet()}

\examples{
qs=getExampleQseaSet()
}

