\docType{methods}
\name{plotPCA}
\alias{plotPCA}
\alias{plotPCAfactors}
\alias{plotPCA,qseaPCA-method}
\alias{plotPCAfactors,qseaPCA-method}

\title{Plots for Principle Component Analysis (PCA) in QSEA}
\description{
The principle components can be depicted using the plotting 
methods plotPCA and plotPCAfactors}

\usage{
\S4method{plotPCA}{qseaPCA}(object,plotComponents=c(1,2), fgColor="black",
    bgColor = "white", legend, plotLabels=TRUE, radius=5, labelOffset=.5,
    labelPos=1, labelAdj, labelColor="black", cex=1, \dots)

\S4method{plotPCAfactors}{qseaPCA}(object,plotComponents=c(1,2), 
    fgColor="black",bgColor = "white", plotTopLabels=100, labelsOfInterest, 
    radius=1, labelOffset=.5,labelPos=1,labelColor="black", cex=1, \dots)
}
\arguments{
\item{object}{the qseaPCA object, resulting from the getPCA function}
\item{plotComponents}{vector of the two components of the PCA}
\item{fgColor}{vector of foreground colors for the circles}
\item{bgColor}{vector of background colors for the circles}
\item{legend}{add a legend to the plot}
\item{plotLabels}{if set TRUE, the labels of 
        the samples are written in the plot}
\item{radius}{defines the size of the plotted circles}
\item{labelOffset}{defines the offset of the labels to the circles}
\item{labelPos}{specify position of the labels 
        in the plot (see graphics::text)}
\item{labelAdj}{alternative way to specify position of the 
        labls in the plot (see graphics::text)}
\item{labelColor}{a vector of colors for the labels}
\item{cex}{font size of the labels}
\item{plotTopLabels}{labels of factors with strongest contribution 
        to plotted components are shown}
\item{labelsOfInterest}{vector of factor names that are highlighted 
        and labeled in the plot}
\item{\dots}{further graphical parameters}

}
\value{
The functions return a list with the coordinates of the depicted components
}

\author{
Mathias Lienhard
}
\seealso{plotPCA}
\examples{ 
qs=getExampleQseaSet( repl=5)
pca=getPCA(qs, norm_method="beta")
colors=c(rep("red", 5), rep("green", 5))

plotPCA(pca, bgColor=colors)
#plotPCAfactors is more interesting, if ROIs have been specified in getPCA
plotPCAfactors(pca)

}

