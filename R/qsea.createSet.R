createQseaSet=function(sampleTable,BSgenome,chr.select,Regions,window_size=250)
{            
    ## Parameter Check
    facI=sapply(sampleTable, is.factor)
    sampleTable[,facI] = sapply(sampleTable[,facI, drop=FALSE], as.character)
    rownames(sampleTable)=sampleTable$sample_name
    checkSampleTab(sampleTable)
    #must have: regions or bsgenome
    if( missing(BSgenome) && (missing(Regions) || class(Regions)!="GRanges" ) )
        stop("Must specify a BSgenome library or a GRanges \"Regions\" object.")
    message("==== Creating qsea set ====")
    # sort chromosomes
    parameters=list(window_size=window_size)    
    if(!missing(chr.select)){
        chr.select=mixedsort(as.character(chr.select))
        message("restricting analysis on ", paste(chr.select, collapse=", "))
    } 
    if (!missing(BSgenome)) {
        parameters[["BSgenome"]] = BSgenome
    }

    ## get Genomic Regions
    if (missing(Regions)){
        message("Dividing selected chromosomes of ",BSgenome , " in ", 
            window_size, "nt windows")
        Regions=makeGenomeWindows(BSgenome, chr.select, window_size)
    }else{
        message("Dividing provided regions in ", window_size, "nt windows")
        Regions=subdivideRegions(Regions, chr.select,window_size, BSgenome)
    }
    chrN=seqlevels(Regions)
    
    cyM=matrix(2,nrow(sampleTable), length(chrN), 
        dimnames=list(sampleTable$sample_name,chrN ))
    sexIdx=match(c("sex", "gender"), names(sampleTable))
    sexIdx=head(sexIdx[!is.na(sexIdx)],1)
    if(length(sexIdx)==0){
        sex=rep("unknown", nrow(sampleTable))
        if(any(c("X", "Y","chrX", "chrY") %in% chrN))
            warning("no column \"sex\" or \"gender\"found in sampleTable, ",
                "assuming heterozygosity for all selected chromosomes")
    }else{
        sex=sampleTable[,sexIdx]
        mIdx=sex %in% c("M","m", "male")
        fIdx=sex %in% c("F","f", "female")
        yIdx=colnames(cyM) %in% c("Y", "chrY")
        xIdx=colnames(cyM) %in% c("X", "chrX")
        if (any(!(mIdx|fIdx)))
            warning("unknown sex specified sampleTable for sample(s) ",
                paste(rownames(cyM)[!(mIdx|fIdx)],collapse=", "),
                "; assuming heterozygosity for all selected chromosomes")
        cyM[mIdx,yIdx|xIdx]=1
        cyM[fIdx,yIdx]=0
    }
    qs=new('qseaSet', sampleTable=sampleTable,
                regions=Regions,
                zygosity=cyM,
                count_matrix=matrix(), 
                cnv=GRanges(), #logFC
                parameters=parameters, 
                libraries=list(), 
                enrichment=list()
    )

}

addNewSamples<-function(qs, sampleTable, force=FALSE, parallel=FALSE){
    #check consistency of sample tables 
    if( ncol(sampleTable) != ncol(getSampleTable(qs)) ||
        ! all.equal(colnames(sampleTable) , colnames(getSampleTable(qs))))
        stop("columns of sampleTable must match sample table in")
    old_samples=getSampleNames(qs)
    facI=sapply(sampleTable, is.factor)
    sampleTable[,facI] = sapply(sampleTable[,facI, drop=FALSE], as.character)
    rownames(sampleTable)=sampleTable$sample_name
    new=character(0)
    for (sa in rownames(sampleTable)){
        if(sa %in% old_samples){
            if(any(sampleTable[sa,]!=getSampleTable(qs,sa), na.rm=TRUE))
                stop("sample ",sa, " is already contained in qs, and differs")
            else message(sa ," already contained in qs")
        }else new=c(new,sa)
    }
    if(length(new)==0)
        stop("no new samples found in sampleTable")
    sampleTable=sampleTable[new,]
    checkSampleTab(sampleTable)
    message("adding ",length(new), " new samples")
    w=FALSE
    if(hasCNV(qs)){
        warning("CNVs for new samples can't be added to qsea set")
        w=TRUE
    }
    if(hasEnrichment(qs) ){
        warning("Enrichment parameters for new samples can't be added to qsea set")
        w=TRUE
    }
    if (w && !force)
        stop("use force=TRUE to force adding new samples ", 
            "by first removeing components that can't be added individually ",
            "for new samples")
    qs=setSampleTable(qs, rbind(getSampleTable(qs), sampleTable))
    qs=setCNV(qs,GRanges()) #remove CNV
    qs=setEnrichment(qs,list()) #remove Enrichment
    chrN=mixedsort(levels(seqnames(getRegions(qs))))
    
    cyM=matrix(2,nrow(sampleTable), length(chrN), 
        dimnames=list(sampleTable$sample_name,chrN ))
    sexIdx=match(c("sex", "gender"), names(sampleTable))
    sexIdx=head(sexIdx[!is.na(sexIdx)],1)
    if(length(sexIdx)==0){
        sex=rep("unknown", nrow(sampleTable))
        if(any(c("X", "Y","chrX", "chrY") %in% chrN))
            warning("no column \"sex\" or \"gender\"found in sampleTable, ",
                "assuming heterozygosity for all selected chromosomes")
    }else{
        sex=sampleTable[,sexIdx]
        mIdx=sex %in% c("M","m", "male")
        fIdx=sex %in% c("F","f", "female")
        yIdx=colnames(cyM) %in% c("Y", "chrY")
        xIdx=colnames(cyM) %in% c("X", "chrX")
        if (any(!(mIdx|fIdx)))
            warning("unknown sex specified sampleTable for sample(s) ",
                paste(rownames(cyM)[!(mIdx|fIdx)],collapse=", "),
                "; assuming heterozygosity for all selected chromosomes")
        cyM[mIdx,yIdx|xIdx]=1
        cyM[fIdx,yIdx]=0
    }
    qs=setZygosity(qs,rbind(getZygosity(qs), cyM))

    if(nrow(getCounts(qs))==0)
        return(qs)
        
    #add counts for new samples
    fragment_length=getParameters(qs, "fragment_length")
    uniquePos=getParameters(qs, "uniquePos")
    minMapQual=getParameters(qs, "minMapQual")
    paired=getParameters(qs, "paired")

    fname_idx=which(names(sampleTable)=="file_name" )[1]
    # get the coverage
    if(parallel) {  
        BPPARAM=bpparam()
        message("Scanning ",bpworkers(BPPARAM) , " files in parallel")
    }else
        BPPARAM=SerialParam()
    
    coverage=unlist(bplapply(X=sampleTable[,fname_idx],FUN=getCoverage, 
            Regions=getRegions(qs),fragment_length=fragment_length,
            minMapQual=minMapQual, paired=paired, uniquePos=uniquePos, 
            BPPARAM=BPPARAM ),FALSE, FALSE)
    
    libraries=rbind(getLibrary(qs,"file_name"), 
        matrix(unlist(coverage[seq(2,length(coverage),by=2)],FALSE,FALSE),
        nrow(sampleTable),6, byrow=TRUE, dimnames=list(sampleTable$sample_name, 
            c("total_fragments", "valid_fragments","library_factor", 
                "fragment_length", "fragment_sd", "offset"))))

    coverage=cbind(getCounts(qs), 
        matrix(unlist(coverage[seq(1,length(coverage),by=2)],FALSE,FALSE), 
        ncol=nrow(sampleTable), byrow=FALSE, 
        dimnames=list(NULL, sampleTable$sample_name)))

    qs=setCounts(qs,count_matrix=coverage)
    qs=setLibrary(qs, "file_name", libraries)
    #addOffset & estimateLibraryFactors
    hasLibF=any(!is.na(getLibrary(qs, "file_name")[,"library_factor"]))
    hasOffset=any(!is.na(getOffset(qs)))
    if(hasLibF){
        qs=addLibraryFactors(qs)
        if(hasOffset)
            qs=addOffset(qs)
    }
    return(qs)
}
#for each window, calculate the sequence preference from low coverage input seq
addSeqPref<-function(qs, seqPref,file_name, fragment_length, paired=FALSE, 
        uniquePos=TRUE, alpha=0.05, pseudocount=5, cut=3){
    if(! missing(seqPref) ){
        if(class(seqPref)!="numeric")
            stop("Please provide sequence preference as log2FC")
        if(length(seqPref) != length(getRegions(qs)))
            stop("length of provided sequence preference",
                " does not match number of windows")
    }else{
        if(missing(file_name))
            stop("please specify column with file names for ",
                "sequence preference estimation")
        if(paired)
            fragment_length=NULL
        if(missing(fragment_length))
            fragment_length=getFragmentLength(qs)[1]
        seqPref=findSeqPref(qs=qs, file_name="input", 
            fragment_length=fragment_length, paired=paired, uniquePos=uniquePos,
            alpha=alpha, pseudocount=pseudocount, cut=cut)    
        qs=setLibrary(qs, "SeqPref",seqPref$libraries)
        seqPref=seqPref$seqPref
    }
    qs=addRegionsFeature(qs, "seqPref",seqPref)
    if(! all(is.na(getOffset(qs,scale="rpkm") ) ))
        warning("Consider recalculating offset based ",
                "on new sequence preference")
    return(qs)

}
getFragmentLength<-function(qs){
            lib=getLibrary(qs,"file_name")
            if(is.null( lib))
                stop("no fragment length found... either specify parameter ",
                    "\"fragment_length\" or run \"addCoverage()\" first.")
            r=rank(lib[,"fragment_length"] + lib[,"fragment_sd"])
            idx=which.min(abs(r-length(r)/2)) #~which.median
            fragment_length=lib[idx,"fragment_length"]
            fragment_sd=lib[idx,"fragment_sd"]
            message("selecting fragment length (", round(fragment_length,2), 
                ") and sd (", round(fragment_sd,2),") from sample ", 
                getSampleNames(qs,idx), " as typical values")
            return(c(fragment_length,fragment_sd))
}

addPatternDensity<-function(qs, pattern,name, fragment_length, fragment_sd,
        patternDensity, fixed=TRUE,masks=c("AGAPS","AMB", "RM", "TRF")[1:2]){
    if(missing(patternDensity) ){
        BSgenome=getGenome(qs)
        if(is.null(BSgenome))
            stop("Pattern density calculation requires BSgenome")
        if(missing(pattern) )    
            stop("please provide sequence pattern for density estimation")
        if(missing(name))
            name=pattern
        if(missing(fragment_length)){
            fl=getFragmentLength(qs)
            fragment_length=fl[1]
            fragment_sd=fl[2]
        }else{
            if(missing(fragment_sd))
                fragment_sd=0
        }                
        patternDensity=estimatePatternDensity(Regions=getRegions(qs), 
            pattern=pattern,BSgenome=BSgenome, fragment_length=fragment_length, 
            fragment_sd=fragment_sd, fixed=fixed, masks=masks) 
    }
    #if(! all(is.na(getOffset(qs,scale="rpkm") ) ))
    #    warning("Consider recalculating offset based on new pattern density")
    if(missing(name)){
        stop("please provide a name for the pattern")
    }
    #add density of pattern
    addRegionsFeature(qs,paste0(name,"_density"), patternDensity )
}

addCoverage<-function(qs, fragment_length, uniquePos=TRUE, minMapQual=1, 
        paired=FALSE, parallel=FALSE){
    sampleTable=getSampleTable(qs)
    Regions=getRegions(qs)
    if(paired)
        fragment_length=NULL
    if(missing(fragment_length))
        stop("for unpaired reads, please specify fragment length")
    fname_idx=which(names(sampleTable)=="file_name" )[1]
    # get the coverage
    if(parallel) {  
        BPPARAM=bpparam()
        message("Scanning ",bpworkers(BPPARAM) , " files in parallel")
    }else
        BPPARAM=SerialParam()
    
    coverage=unlist(bplapply(X=sampleTable[,fname_idx],FUN=getCoverage, 
            Regions=Regions,fragment_length=fragment_length,
            minMapQual=minMapQual, paired=paired, uniquePos=uniquePos, 
            BPPARAM=BPPARAM ),FALSE, FALSE)
    
    libraries=matrix(unlist(coverage[seq(2,length(coverage),by=2)],FALSE,FALSE),
        nrow(sampleTable),6, byrow=TRUE, dimnames=list(sampleTable$sample_name, 
            c("total_fragments", "valid_fragments","library_factor", 
                "fragment_length", "fragment_sd", "offset")))

    coverage=matrix(unlist(coverage[seq(1,length(coverage),by=2)],FALSE,FALSE), 
        ncol=nrow(sampleTable), byrow=FALSE, 
        dimnames=list(NULL, sampleTable$sample_name))
    param=list(uniquePos=TRUE, minMapQual=minMapQual, paired=paired)
    qs=addParameters(qs,param)
    qs=setCounts(qs,count_matrix=coverage)
    setLibrary(qs, "file_name", libraries)
}

addLibraryFactors<-function(qs, factors,...){
    if(!hasCounts(qs))
        stop("No read counts found in qsea set. Run addCoverage first.")
    if(missing(factors)){
        message("deriving TMM library factors for ", 
            length(getSampleNames(qs))," samples" )
        factors=estimateLibraryFactors(qs, ...)        
    }else if (!is.numeric(factors) || (length(factors)!=1 && 
            length(factors)!=length(getSampleNames(qs)) ) ){
        stop("Number of factors does not mathch number of samples.")
    }
    lib=getLibrary(qs, "file_name")
    lib[, "library_factor"]=factors
    setLibrary(qs, "file_name", lib)
}

estimateLibraryFactors<-function(qs,trimA=c(.5,.99), trimM=c(.1,.9),
        doWeighting=TRUE, ref=1, plot=FALSE, nReg=500000){
    if( length(trimM) != 2 || trimM[1]>=trimM[2] || trimM[1]<0 || trimM[2] > 1)
        stop("invalid trimM parameter for TMM: ", paste(trimM))
    if( length(trimA) != 2 || trimA[1]>=trimA[2] || trimA[1]<0 || trimA[2] > 1)
        stop("invalid trimA parameter for TMM: ", paste(trimA))
    n=length(getSampleNames(qs))
    if(n==1) return (1)

    tmm=numeric(n)
    if (missing(ref)) ref=1
    if(is.character(ref))
        ref=which(getSampleNames(qs)==ref)
    if(is.na(ref) | !is.numeric(ref) |ref<1 | ref > n )
        stop("invalid reference sample for TMM (",ref,")")
    others=seq_len(n)[-ref]
    libsz=getLibSize(qs, normalized=FALSE)
    normM=list(scaled=c("zygosity","cnv", "preference"))
    if(ncol(mcols(getRegions(qs)))>=1){
        wd=which(!is.na(mcols(getRegions(qs))[,1]))
        wd=wd[seq(1,length(wd), ceiling(length(wd)/nReg))]
    }else{
        tReg = length(getRegions(qs))
        wd=seq(1,tReg, ceiling(tReg/nReg))
    }
    values=getNormalizedValues(qs,methods=normM,
        windows=wd, 
        samples=getSampleNames(qs) )
    values[values<1]=NA

    if(plot){
        row=ceiling(sqrt(n-1))
        col=ceiling((n-1)/row)
        par(mfrow=c(row, col))
        cr=colorRampPalette(
            c("white","skyblue", "blue", "yellow", "orange","red", "darkred"))
    }
    for(i in others){
        a=log(values[,i])+log(values[,ref])-log(libsz[ref])-log(libsz[i])
        m=log(values[,i])-log(values[,ref])+log(libsz[ref])-log(libsz[i])
        thA=quantile(a, trimA, na.rm=TRUE)
        thM=quantile(m, trimM, na.rm=TRUE)
        sel=which(a>thA[1] & a < thA[2] & m>thM[1] &m<thM[2])
        if(doWeighting){
            w=1/(1/values[sel,ref]+1/values[sel,i]-1/libsz[ref]-1/libsz[i])
            tmm[i]=sum(m[sel]*w)/sum(w)
        }else{
            tmm[i]=mean(m[sel])
        }
        if(plot){
            smoothScatter(a,m, pch=20, colramp=cr, ylab="M", xlab="A", 
                main=paste(getSampleNames(qs, c(i,ref)), collapse=" vs "))
            abline(h=tmm[i], col="red", lwd=2)
            rect(thA[1], thM[1], thA[2], thM[2])
            
        }
    }
    if(any(is.na(tmm))){
        warning("tmm normalization faild")
        return(rep(1,n))
    }else
        return(exp(tmm-mean(tmm)))
}


addEnrichmentParameters<-function(qs, enrichmentPattern, signal,
        windowIdx, min_wd=5,bins=seq(.5,40.5,1)){
    if(missing (windowIdx))
        stop("please specify the windows for enrichment analysis")
    if(missing (signal))
        stop("please provide a signal matrix for the specified windows")
    if(missing(enrichmentPattern))
            enrichmentPattern=getParameters(qs, "enrichmentPattern")
    if(is.null(enrichmentPattern))
            stop("please specify sequence pattern for enrichment analysis")

    enrichment=estimateEnrichmentLM(qs,windowIdx=windowIdx, signal=signal, 
        min_wd=min_wd,bins=bins, pattern_name=enrichmentPattern)
    parameters=fitEnrichmentProfile(enrichment$factors, enrichment$density, 
        enrichment$n, minN=1)
    addParameters(qs,list(enrichmentPattern=enrichmentPattern))
    setEnrichment(qs, c(list(parameters=parameters),enrichment))
}

subdivideRegions<-function(Regions, chr.select,window_size, BSgenome){
    if(missing(chr.select))
        chr.select=mixedsort(seqlevels(Regions))
    Regions=Regions[as.vector(seqnames(Regions)) %in% chr.select]
    #order according to chr.select
    seqlevels(Regions)=chr.select
    #merge overlapping windows
    Regions=reduce(sort(Regions))
    ranges=ranges(Regions)

    #resize (enlarge) to a multiple of window_size
    pos=apply(FUN=function(x) seq(from=x[1], to=x[2], by=window_size), 
        X=as.data.frame(ranges),MARGIN=1)
    if(class(pos)=="matrix"){
        n=nrow(pos)
        pos=as.vector(pos)
        chr=rep(as.character(seqnames(Regions)), each=n)
    }else {
        n=sapply(X=pos,FUN=length)
        pos=unlist(pos)
        chr=rep(as.character(seqnames(Regions)), times=n)
    }
    if(!missing(BSgenome)){
        chr_length=seqlengths(get(ls(paste("package:", BSgenome, sep=""))))
        seqinfo = Seqinfo(names(chr_length),chr_length, NA, BSgenome)
    }else{
        seqinfo=seqinfo(Regions)
    }
    if(!missing(chr.select)){
        seqinfo=seqinfo[chr.select]
    }

    return(GRanges(seqnames=chr, IRanges(start=pos, width=window_size), 
        seqinfo=seqinfo))
}

addOffset<-function(qs,enrichmentPattern , maxPatternDensity=0.01,offset){
    if(missing(offset)){

        if(missing(enrichmentPattern))
            enrichmentPattern=getParameters(qs, "enrichmentPattern")
        if(is.null(enrichmentPattern))
            stop("please specify sequence pattern for enrichment analysis")
        offset=estimateOffset(qs,enrichmentPattern, maxPatternDensity)
    }
    lib=getLibrary(qs, "file_name")    
    lib[,"offset"]=offset
    qs= addParameters(qs,list(enrichmentPattern=enrichmentPattern))
    setLibrary(qs, "file_name", lib)
}




checkSampleTab<-function(sampleTable){
    if(missing(sampleTable) || !is.data.frame(sampleTable))
        stop("Must specify a sample table (data.frame)")
    m=match(c("sample_name", "file_name", "group"),names(sampleTable))
    if(any(is.na(m)))
        stop("sample table must contain \"sample_name\", ",
            "\"file_name\" and \"group\" columns")    
    sname_idx=m[1]
    fname_idx=m[2]
    group_idx=m[3]
    if(length(unique(sampleTable[,sname_idx]))!=nrow(sampleTable))
        stop("Sample names must be unique")
    
    wigFiles=FALSE
    bamFiles=FALSE
    errorMsg=""
    for(sNr in seq_len(nrow(sampleTable))){
        if(file.exists(as.character(sampleTable[sNr,fname_idx]))){
            tmp=strsplit(sampleTable[sNr,fname_idx], ".", fixed=TRUE)[[1]]
            if (tmp[length(tmp)] %in% c("gz","zip","bz2")){
                tmp=tmp[-length(tmp)]
            }
            ext=tmp[length(tmp)]
            if(ext %in% c("wig","bw","bigwig")){
                wigFiles=TRUE
            }else if(ext %in% c("bam","BAM","sam", "SAM")){
                bamFiles=TRUE
            }else{
                errorMsg=paste(errorMsg,"\nFiletype unknown:", 
                sampleTable[sNr,fname_idx])
            }
        }else{errorMsg=paste(errorMsg,"\nFile not found:", 
            sampleTable[sNr,fname_idx])}
    }
    if (errorMsg !=""){
        stop("Input file check failed:", errorMsg)
    }
}

