
makeGenomeWindows<-function(BSgenome=NULL, chr.select=NULL, window_size=250){
    if(is.null(BSgenome)){stop("Please specify genome")}
    dataset=get(ls(paste("package:", BSgenome, sep="")))
    chr_length=seqlengths(dataset)
    if(! is.null(chr.select) ) {
        if(length(intersect(names(chr_length),chr.select))==0) 
            stop ("specified chromosomes not found in ", BSgenome)
        chr=as.vector(na.omit(
            mixedsort(intersect(names(chr_length),chr.select))))
        if(length(chr)!=length(chr.select)) 
            warning("not all specified chromosomes found in ",BSgenome, 
                "\nusing only ",chr)
    }else{
        chr=as.vector(na.omit(mixedsort(names(chr_length))))
    }
    chr_length=chr_length[match(chr,names(chr_length))]
    filter=chr_length<window_size
    if (any (filter)){
        warning("neglecting chromosomes shorter than window size:\n", 
            paste(chr[filter],collapse=", "))
        chr=chr[!filter]
        chr_length=chr_length[!filter]
    }
    nr_wd=floor(chr_length/window_size)
    wd_start=unlist(lapply(FUN=seq, X=chr_length-window_size, 
        from=1,by=window_size))    
    seqinfo = Seqinfo(chr,chr_length, NA, BSgenome)
    return(Regions=GRanges(seqnames=rep(factor(chr),nr_wd), 
        ranges=IRanges(start=wd_start, width=window_size), seqinfo=seqinfo))
}


estimatePatternDensity <- function(Regions=NULL, pattern="CG", BSgenome=NULL,
        fragment_length=300, fragment_sd=0, fixed=TRUE ){
    ## Get the genomic positions of the sequence pattern
    message("Get genomic positions of \"",pattern,"\" ...",sep="")
    dataset=getBSgenome(BSgenome, masked=FALSE)    
    patternDensity=rep(NA, length(Regions))
    window_size=width(Regions[1])
    #rcpattern=reverseComplement(pattern)
    len=seqlengths(Regions)
    genomeRange=GRanges(names(len), IRanges(1, len), strand="+")
    #genome=getSeq(dataset, name=genomeRange)
    if(fragment_sd>0){
        effect_len=seq(    2,fragment_length+1+2*fragment_sd,by=2)
        effect_len=c(rev(effect_len)-1,effect_len )
        effect=1-pnorm(q=effect_len, mean=fragment_length, sd=fragment_sd)
    }

    for (chr in seqlevels(Regions)){
        message("searching ",chr," for \"",pattern,"\"...")
        sel=seqnames(Regions)==chr
        chr_seq=dataset[[chr]]
        chr_seq<-maskMotif(chr_seq, "N")
        pIdx=start(matchPattern(pattern=DNAString(pattern),
            subject=chr_seq, fixed=fixed))
        message("found ", length(pIdx), " occurances of ",
            pattern , " in ", chr)
        message("estimating expected number of ",pattern, 
            "s per fragment for windows of ", chr,"...")
        if(fragment_sd<=0){
            pInf=IRanges( pIdx-ceiling(fragment_length/2),
                pIdx+floor(fragment_length/2) ) 
            pCov=coverage(pInf)
        }else{
            startPos=pIdx-(ceiling(fragment_length/2+fragment_sd))
            pCov=numeric(seqlengths(Regions)[chr])
            for(i in seq_along(effect)){
                idx=startPos[startPos > -i]+i
                pCov[idx]=pCov[idx]+effect[i]
            }
        }
        v=Views(pCov, start= ranges(Regions[sel]) )
        md=viewSums(v)/window_size
        #set density of regions with more that 50% masked 
        #(or hardmasked by N) to NA
        masked=reduce(nir_list(collapse(masks(chr_seq)))[[1]])
        midx=overlapsAny(Regions[sel],GRanges(chr,masked),
            minoverlap=window_size/2)
        md[midx]=NA
        patternDensity[which(sel)]=md    
        message(" ...done")
    }
    return(patternDensity)    
}
