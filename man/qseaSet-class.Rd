\Rdversion{1.1}
\docType{class}
\name{qseaSet-class}
\alias{qseaSet-class}
\alias{qseaSet}
\alias{show,qseaSet-method}
\alias{getSampleTable,qseaSet-method}
\alias{getOffset,qseaSet-method}
\alias{getWindowSize,qseaSet-method}
\alias{hasCNV,qseaSet-method}
\alias{getSampleNames,qseaSet-method}
\alias{getSampleGroups,qseaSet-method}
\alias{getChrNames,qseaSet-method}
\alias{getParameters,qseaSet-method}
\alias{getRegions,qseaSet-method}
\alias{getLibSize,qseaSet-method}
\alias{getNormFactors,qseaSet-method}
\alias{getCNV,qseaSet-method}
\alias{getCounts,qseaSet-method}
\alias{getPloidity,qseaSet-method}
\alias{setPloidity,qseaSet-method}

\alias{getWindowSize}
\alias{getSampleTable}
\alias{getSampleGroups}
\alias{getSampleNames}
\alias{hasCNV}
\alias{getCNV}
\alias{getOffset}
\alias{getChrNames}
\alias{getCounts}
\alias{getLibSize}
\alias{getNormFactors}
\alias{getParameters}
\alias{getRegions}
\alias{getPloidity}
\alias{setPloidity}

\title{qseaSet class and its methods}

\description{The qseaSet class is used in qsea to store information about the 
coverage, the dependent organism, the chromosomes included in the input file, 
the length of the included chromosomes (automatically loaded), the number of 
regions, and optionally CNV information.
}

\section{Slots}{
    \describe{
    \item{\code{sampleTable}:}{Object of class \code{"data.frame"}: 
    the sample table}
    \item{\code{count_matrix}:}{Object of class \code{"matrix"}: 
    matrix containing the coverage for all samples}
    \item{\code{ploidity}:}{Object of class \code{"matrix"}: 
    matrix containing the ploidity for all chromosomes and all samples}
    \item{\code{regions}:}{Object of class \code{"GenomicRanges"}: 
    the genomic regions for the coverage matrix}
    \item{\code{parameters}:}{Object of class \code{"list"}: 
    the parameter list used to create this object}
    \item{\code{cnv}:}{Object of class \code{"GenomicRanges"}: 
    CNV ranges and logFCs}
    \item{\code{enrichment }:}{Object of class \code{"list"}:
    parameters of the sequence pattern enrichment analysis}
    \item{\code{libraries }:}{Object of class \code{"matrix"}:
    parameters of the sequencing libraries}
    }
}

\section{Methods}{
    \describe{
    \item{getSampleTable}{\code{signature(object = "qseaSet")}: 
    extracts the sample table of a qsea set}
    \item{getSampleNames}{\code{signature(object = "qseaSet")}: 
    extracts the sample names of a qsea set}
    \item{getSampleGroups}{\code{signature(object = "qseaSet")}: 
    extracts the sample groups of a qsea set}
    \item{getChrNames}{\code{signature(object = "qseaSet")}: 
    returns the analysed chromosomes}
    \item{getCounts}{\code{signature(object = "qseaSet")}: 
    extracts the count matrix a qsea set}
    \item{getRegions}{\code{signature(object = "qseaSet")}: 
    extracts the regions object of a qsea set}
    \item{getParameters}{\code{signature(object = "qseaSet")}: 
    extracts the parameter list of a qsea set}
    \item{getLibSize}{\code{signature(object = "qseaSet")}: 
    extracts the library size (eg the total number of read counts per sample)}
    \item{getNormFactors}{\code{signature(object = "qseaSet")}: 
    extracts the list with the different normalization factors}
    \item{hasCNV}{\code{signature(object = "qseaSet")}: 
    TRUE if CNV information is present, FALSE otherwise}
    \item{getCNV}{\code{signature(object = "qseaSet")}: 
    extracts the CNV regions and logFCs}
    \item{getOffset}{\code{signature(object = "qseaSet")}: 
    extracts offset of rpkm scaled background reads}
    \item{getWindowSize}{\code{signature(object = "qseaSet")}: 
    returns the window size of the object}
    \item{getPloidity}{\code{signature(object = "qseaSet")}: 
    returns the ploidity matrix of the object}
    \item{setPloidity}{\code{signature(object = "qseaSet", ploidityMatrix)}: 
    sets the ploidity matrix, and resets CNV}
    }
}
\author{Matthias Lienhard}
\examples{
showClass("qseaSet")
}
\keyword{classes}

