\Rdversion{1.1}
\docType{class}
\name{qseaGLM-class}
\alias{qseaGLM-class}
\alias{qseaGLM}
\alias{show,qseaGLM-method}
\alias{getParameters,qseaGLM-method}
\alias{getSampleNames,qseaGLM-method}

\title{qseaGLM class and its methods}

\description{
    The qseaGLM class is used in qsea to store fitted coefficients of the GLM.
}

\section{Slots}{
    \describe{
    \item{\code{fullModelDesign}:}{design matrix of full model }
    \item{\code{fullModel }:}{list containing parameters and fitted coefficients
of full model}
    \item{\code{parameters}:}{list of parameters used to create the object}
    \item{\code{contrast}:}{list of lists containing parameters and the fitted 
model coeficients of the reduced models }
    \item{\code{windows}:}{vector of window indices, for which GLMs have been 
fitted}
    }
}

\author{Matthias Lienhard}
\examples{
showClass("qseaGLM")
}
\keyword{classes}




