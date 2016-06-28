\name{regionStats}
\alias{regionStats}

\title{Counts the Windows in Regions of Interest}
\description{This function takes a list of window indices and a list of ROIs
 and counts the number of overlapping windows}
\usage{
regionStats(qs, subsets = list(covered = which(rowSums(getCounts(qs)) >= 20)),
 ROIs = list(), minoverlap = 1, maxgap = 0)
}
\arguments{
\item{qs}{A qsea Set object}
\item{subsets}{A list of window indices}
\item{ROIs}{A list of Regions of Interest}
\item{minoverlap}{Passed to findOverlaps}
\item{maxgap}{Passed to findOverlaps}

}
\value{a matrix, containing the total number of windows overlapping the ROIs
 and the numbers of windows from the subset list overlapping ROIs}

\author{
Mathias Lienhard
}

\seealso{findOverlaps}

\examples{
qs=getExampleQseaSet()
#as an example, we analyze the fraction of reads covered by at least 10
#or at least 20 reads, for bins of CpG density
ROIs=list()
regs=getRegions(qs)
cpg=getRegions(qs)$CpG_density
bins=seq(0,30,5)
for(i in 1:(length(bins)-1)){
    n=paste0(bins[i],"-",bins[i+1]," CpGs")
    ROIs[[n]]=regs[which(cpg>=bins[i] & cpg < bins[i+1])]
}
subsets = list(
    ">10" = which(rowSums(getCounts(qs)) >= 10),
    ">20" = which(rowSums(getCounts(qs)) >= 20))
coverage_stats=regionStats(qs, subsets, ROIs)
coverage_stats_rel=coverage_stats[,-1]/coverage_stats[,1]
x=barplot(t(coverage_stats_rel)*100,ylab="fraction of windows[\%]",
    beside=TRUE, legend=TRUE, las=2, args.legend=list(x="topleft"), 
    main="Covered Windows")


}

