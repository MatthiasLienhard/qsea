\name{addSeqPref}
\alias{addSeqPref}

\title{Add sequence preference to qseaSet object}
\description{This function allows to add window specific sequencing preference, 
that can be used by the normalization procedure. This preference can be defined 
by the user, or estimated from sequencing of input libraries.}
\usage{
addSeqPref(qs, seqPref,file_name, fragment_length, paired=FALSE, 
uniquePos=TRUE, alpha=0.05, pseudocount=5, cut=3)
}
\arguments{
\item{qs}{a qseaSet object}
\item{seqPref}{A vector with predefined sequencing preference for each window. 
Values are interpreted as log2 ratios relative to normal/average sequencing 
preference. }
\item{file_name}{alternatively, the sequencing preference can be estimated 
from input sequencing. In this case, provide the column of the sample table 
that contains the file names for input sequencing alignment or coverage files.}
\item{fragment_length}{for single end data, provide the expected 
fragment length}
\item{paired}{if set to TRUE, data is considered to be paired end sequencing, 
and the actual fragments size is used. }
\item{uniquePos}{if set to TRUE, fragments with same position and orientation 
are considered to be PCR duplicates and replaced by one representative. }
\item{alpha}{currently ignored}
\item{pseudocount}{this value is added to the coverage of each window, 
to smooth the estimates.}
\item{cut}{absolute log2 value threshold for windows to be excluded from 
later analysis due to extreme preference values.}
}
\value{the function returns the qseaSet object, extended by the sequencing 
preference for all genomic windows.}

\author{
Mathias Lienhard
}

