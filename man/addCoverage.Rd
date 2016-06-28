\name{addCoverage}
\alias{addCoverage}


\title{Import sequencing data}
\description{This function imports the alignment files (in sam/bam format) and 
counts the reads per genomic window or directly imports 
coverage files (in wiggle/bigwiggle format)}
\usage{
addCoverage(qs, fragment_length, uniquePos=TRUE, minMapQual=1, paired=FALSE)
}
\arguments{
\item{qs}{qseaSet object, e.g. produced by the createQseaSet() function}
\item{fragment_length}{
For single end data, provide the expected fragment length}
\item{paired}{If set to TRUE, data is considered to be paired end sequencing, 
and the actual fragments size is used. }
\item{uniquePos}{If set to TRUE, fragments with same position and orientation 
are considered to be PCR duplicates and replaced by one representative. }
\item{minMapQual}{The minimal mapping quality for reads to be considered. 
Note that the definition of mapping quality depends on the alignment tool.}
}
\details{The coverage is imported from the files specified in the file_name 
column of the sample table, provided for the createQseaSet() function. 
In case of alignment files, the reads are counted for the window at the center 
of the sequencing fragment. For single end data, 
Filetypes is detected automatically from the file suffix.}
\value{The function returns the qseaSet object, extended by the number of 
reads per window for all samples}

\author{
Mathias Lienhard
}

\seealso{crateQseaSet}

\examples{
library("BSgenome.Hsapiens.UCSC.hg19")

bam_hESCs_1 = system.file("extdata", 
    "hESCs.MeDIP.Rep1.chr22.bam", package="MEDIPSData")
bam_hESCs_2 = system.file("extdata", 
    "hESCs.MeDIP.Rep2.chr22.bam", package="MEDIPSData")
sample_table=data.frame(sample_name=paste0("hESCs_", 1:2), 
    file_name=c(bam_hESCs_1,bam_hESCs_2), 
    group=rep("hESC",2), stringsAsFactors=FALSE)
qseaSet=createQseaSet(sampleTable=sample_table, 
        BSgenome="BSgenome.Hsapiens.UCSC.hg19", 
        chr.select="chr22", 
        window_size=500)
qseaSet=addCoverage(qseaSet, fragment_length=300)
}
