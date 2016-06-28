getCoverage<-function(Regions=NULL,file_name=NULL, fragment_length=NULL,
    minMapQual=1, uniquePos=TRUE, paired=FALSE, CpGrange=NULL, CpGpos=NULL){

    if(!paired && is.null(fragment_length) )
        stop("for single end reads, please provide an estimate ",
            "for the fragment length")
    chr.select=as.vector(na.omit(mixedsort(levels(seqnames(Regions)) )))
    tmp=strsplit(file_name, ".", fixed=TRUE)[[1]]
    if (tmp[length(tmp)] %in% c("gz","zip","bz2")){
        tmp=tmp[-length(tmp)]
        tmp[length(tmp)]=paste("zipped",tmp[length(tmp)],sep="_")
    }
    ext=tmp[length(tmp)]
    if(ext %in% c("bam", "BAM")){
        bamindex=file.exists(paste(file_name,".bai", sep=""))
        message("Reading ",ext," alignment ",file_name)
        if(paired){
            scanFlag = scanBamFlag(isPaired=TRUE, isProperPair=TRUE ,
                hasUnmappedMate=FALSE, isUnmappedQuery = FALSE, 
                isFirstMateRead = TRUE, isSecondMateRead = FALSE)
            fields=c("mapq" ,"rname", "pos", "strand","qwidth","isize")
        }else{
            scanFlag = scanBamFlag(isUnmappedQuery = FALSE)        
            fields=c("mapq" ,"rname", "pos", "strand", "qwidth")
        }
        if(bamindex){
            #read only selected chromosomes            
            sel=GRanges(chr.select,IRanges(1, 536870912))            
            scanParam=ScanBamParam(flag=scanFlag, what = fields, which=sel)    
        }else{
            #read all and selected chromosomes/regions later
            scanParam=ScanBamParam(flag=scanFlag, what = fields)
        }
        Reads = scanBam(file=file_name, param=scanParam)
        #list of lists of factors
        Reads = do.call(rbind,lapply(Reads, as.data.frame))
        Reads= Reads[Reads$mapq >= minMapQual,
            grep("mapq", names(Reads), invert=TRUE) ]
        if(!bamindex){
            Reads=Reads[Reads$rname %in% chr.select,]
        }
        totalR=nrow(Reads)
        uniqueR=totalR
        message("Total number of imported ",
            ifelse(paired,"read pairs","reads"),": ", totalR, sep="")
        if(uniquePos){
            message(paste("Extract ",
                ifelse(paired,"read pairs","reads"),
                " with unique position...", sep=""))
            Reads=unique(Reads)
            uniqueR=nrow(Reads)
            message("Number of unique ",
                ifelse(paired,"read pairs","reads"),": ", uniqueR, sep="")
        }
        revStr=Reads$strand=="-"    
        if(paired){
            if(! is.null(fragment_length) ) 
                warning("Ignoring provided fragment length and use ",
                        "paired end information")
            fragment_length=c(mean(abs(Reads$isize)),sd(abs(Reads$isize)) )    
            message("Mean fragment length: ", 
                round(fragment_length[1],3), " nt", sep="")
            message("SD of fragment length: ", 
                round(fragment_length[2],3), " nt", sep="")
            message("Max fragment length: ",
                max(abs(Reads$isize)), " nt", sep="")
            message("Min fragment length: ",
                min(abs(Reads$isize)), " nt", sep="")
        }else{
            Reads$isize=fragment_length[1]
            Reads$isize[revStr]= -fragment_length[1]
            if(length(fragment_length)==1)
                fragment_length=c(1,0.1)*fragment_length
        }
        #for reverse reads, pos is start of second read, and isize is negative
        #obtain start of fragment with:
        Reads$pos[revStr] =Reads$pos[revStr] + Reads$qwidth[revStr] +
            Reads$isize[revStr] 
        if(!is.null(CpGrange)){
            #select reads within CpG range
            if(CpGrange[1]==CpGrange[2])
                message("selecting reads containing ",CpGrange[1]," CpGs")
            else
                message("selecting fragments containing between ",CpGrange[1],
                    " and ",CpGrange[2]," CpGs")
            Reads$nCpG=countOverlaps(GRanges(seqnames=Reads$rname, 
                ranges=IRanges(start=Reads$pos,width=abs(Reads$isize))),CpGpos)
            Reads=Reads[Reads$nCpG>=CpGrange[1] & Reads$nCpG <= CpGrange[2],]
            message("selected ",nrow(Reads),
                ifelse(paired," read pairs"," reads"),
                " within CpG range", sep="")
        }
        #mass point of reads
        Reads = GRanges(seqnames=Reads$rname, 
            ranges=IRanges(start=Reads$pos+abs(round(Reads$isize/2)), width=1),
            strand=Reads$strand)
        message("Calculating short read coverage for genome wide windows...")
        return(list(counts=countOverlaps(Regions, Reads), 
            fragment=c(totalR, uniqueR,fragment_length)))
    }else if (ext %in% c("zipped_wig", "zipped_bw", "wig", "bw")){
        message("Reading coverage file ",file_name)
        if(! is.null(CpGrange)){
            warning("ignoring CpG range parameter for coverage based ",
                "input file",immediate.=TRUE, call.=FALSE)
        }
        sel=GRanges(chr.select,IRanges(1, 536870912))#whole genome
            #this function will warn, if type is not bigwig
            wiggle=rtracklayer::import(file_name, asRangedData=FALSE, which=sel)
        m=as.matrix(findOverlaps(Regions, wiggle, minOverlap=width(Regions)/2))
        #check whether the regions overlap ambigiously
        coverage=numeric(length(Regions))        
        if(any(duplicated(m[,2])))
        {
            warning("Unsufficient resolution in coverage file ",file_name,
                ", returning values several times",immediate.=TRUE, call.=FALSE)
            #return(coverage) #return zeros
        }
        if(any(duplicated(m[,1]))) 
        {
            warning("Defind regions and coverage from ",file_name,
            " do not match uniquely, returning only the first matching value",
                immediate.=TRUE, call.=FALSE)
            #return(coverage) #return zeros
            m=matrix(c(seq_along(coverage),
                findOverlaps(Regions, wiggle,select="first")),
                nrow=length(coverage),2) #return the first match
        }
        if(any(is.na(rowSums(m))) || length(unique(m[,1])) < length(coverage))
        {
            warning("Not all regions are covered in ",file_name,
            ", returning zerro for those regions",immediate.=TRUE, call.=FALSE)
            m=m[!is.na(rowSums(m)),]
        }
        coverage[m[,1]]= mcols(wiggle)$score[m[,2]]
        return(list(counts=coverage,fragment=c(sum(coverage),NA,NA,NA) ))
    }else{
        warning("File format/extension of ",file_name," not supported")
        return(list(counts=numeric(length(Regions)),c(0,0,0,0) ))
    }
}

