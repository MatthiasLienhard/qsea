\name{addPatternDensity}
\alias{addPatternDensity}


\title{Infer sequence pattern density values and add to qseaSet object}
\description{This function estimates the average occurrences of a sequence 
pattern (such as CpG dinucleotides) within the overlapping sequencing 
fragments for each genomic window}
\usage{
addPatternDensity(qs, pattern,name, fragment_length, fragment_sd,
patternDensity, fixed=TRUE, masks=c("AGAPS","AMB", "RM", "TRF")[1:2])
}

\arguments{
\item{qs}{a qseaSet object}
\item{pattern}{actual sequence of the pattern (e.g. \dQuote{CG})},
\item{name}{a name for the sequence pattern(e.g. \dQuote{CpG})},
\item{fragment_length}{the average fragment length to be assumed for pattern 
density estimation. If omitted, this parameter is taken from the 
qseaSet object.}
\item{fragment_sd}{the standard deviation of fragment length to be assumed for 
pattern density estimation. If omitted, this parameter is taken from the 
qseaSet object.}
\item{patternDensity}{this parameter alternatively allows to specify the 
pattern density manually. In this case, please provide a numerical vector, 
containing a value (greater than 0) for each genomic window.}
\item{fixed}{if FALSE, an IUPAC ambiguity code in the pattern can match 
any letter in the reference genome that is associated with
        the code, and vice versa.}
\item{masks}{names of the BSgenome masks to be active.}

}
\value{The function returns the qseaSet object, extended by the pattern 
density for all genomic windows}

\author{
Mathias Lienhard
}

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
qseaSet=addPatternDensity(qseaSet, "CG", name="CpG", fragment_length=300)


}

