\name{addCNV}
\alias{addCNV}

\title{estimate CNV information and add to qseaSet object}
\description{This function adds information on Copy Number Variation (CNV) to 
the qseaSet object, which is used for normalization. 
Sample wise CNV information can either be provided, or estimated from input or 
enrichment sequencing data, by incorporating functions of the HMMcopy package. 
}
\usage{
addCNV(qs,file_name, window_size=1000000, paired=FALSE, fragment_length,cnv, 
mu=log2(c(1/2, 2/3, 1, 3/2,2,3)), normal_idx, plot_dir, MeDIP=FALSE, 
parallel=FALSE)
}
\arguments{
\item{qs}{the qseaSet object}
\item{cnv}{pre-computed CNV information for each sample. If provided, 
    the following parameters are ignored}
\item{file_name}{column name of the sample table for the sequencing files, 
    from which CNV information are computed}
\item{window_size}{window size for CNV analysis}
\item{paired}{are files in file_name column paired end}
\item{fragment_length}{for single end sequencing, 
    provide the average fragment length}
\item{mu}{a priori CNV levels of different states, parameter passed to HMMcopy }
\item{normal_idx}{index of samples which are assumed to be CNV free. 
    The median of these samples serves as "normal" CNV reference level, 
    and CNV are computed relative to this reference level. By default, 
    QSEA looks for samples with "normal" or "control" in its name.}
\item{plot_dir}{If provided, detail CNV plots for each chromosome and 
    each sample are created in the provided directory}
\item{MeDIP}{If set TRUE, QSEA assumes that provided files are methylation 
    enriched sequencing data. In this case, only fragments without CpG 
    dinucleotides are considered. This option allows QSEA to infer CNV 
    information from MeDIP or MDB seq experiments directly}
\item{parallel}{Switch for parallel computing, using BiocParallel}
}
\value{
The qseaSet object, extended by the CNV information
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

#this is an example for computing CNVs from MeDIP data. A very limited example
#however, since the samples do not contain CNVs.
qseaSet=addCNV(qseaSet, fragment_length=300, file_name="file_name", MeDIP=TRUE,
window_size=1000000)


}
\author{
Mathias Lienhard
}

\seealso{HMMsegment}

