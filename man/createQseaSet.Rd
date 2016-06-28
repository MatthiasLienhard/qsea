\name{createQseaSet}
\alias{createQseaSet}

\title{
Prepares a qseaSet Object
}
\description{
    This method prepares the qseaSet object, and prepares genome wide bins. 
Coverage and normalization parameters are added in succeeding functions. 
}
\usage{
createQseaSet(sampleTable,BSgenome, chr.select,Regions, window_size=250 )
}
\arguments{
    \item{BSgenome}{name of BSgenome package}
    \item{Regions}{GRanges object. If specified, 
        only selected regions are processed}
    \item{chr.select}{If specified, only selected chromosomes are processed}
    \item{sampleTable}{
data.frame, containing at least 3 columns: 

the sample names (sample_name), 
paths to alignment or coverage file in 
sam/bam/wiggle/bigwig format (file_name), 
and one or more test condition(s) (group). 

Optionally it may contain a column with alignment or coverage files for CNV 
analysis, and further information in the samples that are of interest for the 
analysis.
}
    \item{window_size}{size for the genome wide bins in base pairs}
}

\value{
An object of class qseaSet, containing the sample and genome information.
}


\author{
Mathias Lienhard
}

\examples{
library("BSgenome.Hsapiens.UCSC.hg19")
bam_hESCs_1 = system.file("extdata", "hESCs.MeDIP.Rep1.chr22.bam", 
    package="MEDIPSData")
bam_hESCs_2 = system.file("extdata", "hESCs.MeDIP.Rep2.chr22.bam", 
    package="MEDIPSData")
samplesTable=data.frame(sample_name=paste0("hESCs_", 1:2), 
file_name=c(bam_hESCs_1,bam_hESCs_2), 
group=rep("hESC",2),stringsAsFactors=FALSE)
qs=createQseaSet(samplesTable, BSgenome="BSgenome.Hsapiens.UCSC.hg19", 
    chr.select="chr22", window_size=500)

}
