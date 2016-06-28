\Rdversion{1.1}
\docType{class}
\name{qseaPCA-class}
\alias{qseaPCA-class}
\alias{qseaPCA}
\alias{show,qseaPCA-method}
\alias{getSampleNames,qseaPCA-method}
\title{qseaPCA class and its methods}

\description{
    The qseaPCA class is used in qsea to store results of the principle 
    component analysis.
}

\section{Slots}{
    \describe{
    \item{\code{svd}:}{singular value decomposition}
    \item{\code{sample_names }:}{names of the samples}
    \item{\code{factor_names}:}{names of the genomic windows involved}
    }
}

\author{Matthias Lienhard}
\examples{
showClass("qseaPCA")
}
\keyword{classes}

