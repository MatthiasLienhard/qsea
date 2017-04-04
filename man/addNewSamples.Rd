\name{addNewSamples}
\alias{addNewSamples}

\title{
Extends an exisiting qseaSet by new samples
}
\description{
This function allows the qseaSet to be extended by new samples, 
provided in the sample table.
}
\usage{
addNewSamples(qs, sampleTable, force=FALSE, parallel=FALSE)
}
\arguments{
    \item{qs}{The qseaSet object to be extended}
    \item{sampleTable}{
        data.frame, describing the samples. Must be in same format as 
        getSampleTable(qs)}
    \item{force}{force adding of new samples, even if existing CNV or 
        enrichment information requires recomputation}
    \item{parallel}{parallel processing of alignment files}
}

\value{
An object of class qseaSet, including the new samples. 
}


\author{
Mathias Lienhard
}

\examples{
library("BSgenome.Hsapiens.UCSC.hg19")
data(samplesNSCLC, package="MEDIPSData")
path=system.file("extdata", package="MEDIPSData")
samples_NSCLC$file_name=paste0(path,"/",samples_NSCLC$file_name )
originalQseaSet=createQseaSet(sampleTable=samples_NSCLC[1:4,], 
        BSgenome="BSgenome.Hsapiens.UCSC.hg19", chr.select="chr22", 
        window_size=500)
originalQseaSet=addCoverage(originalQseaSet, uniquePos=TRUE, paired=TRUE)
qseaSet=addNewSamples(originalQseaSet, samples_NSCLC)
}
