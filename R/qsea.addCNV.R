
findCNV<-function(sampleTable=NULL,BSgenome,chr.select=NULL,
        file_name="CNV_file_name",fragment_length=NULL, uniquePos=TRUE, 
        paired=FALSE, mu =log2(c(1/2, 2/3, 1, 3/2,2,3)), window_size=1000000, 
        normal_idx=NULL, plot_dir=NULL, MeDIP=FALSE,zygosity=NA,
        parallel=FALSE){
    CpGrange=NULL    
    if(!MeDIP)
        message("== Analyzing Copy Number Alterations ==")
    else{
        message("== Analyzing Copy Number Alterations from MeDIP files ==")
        CpGrange=c(0,0)
    }
    #data files
    CNVfname_idx=which(names(sampleTable)==file_name)
    if(length(CNVfname_idx)<1){
        stop("Column for CNV files (\"",file_name,"\") not found;")
    }
    sampleTable$file_names=sampleTable[,CNVfname_idx]
    checkSampleTab(sampleTable)

    file_names=sampleTable[,CNVfname_idx]
    n=length(file_names)

    #get reference columns
    if(is.null(normal_idx)){
        normal_idx=seq_len(n)
        if(length(grep("([Cc]ontrol|[Nn]ormal)", sampleTable$group))>0 )
            normal_idx=
                normal_idx[grep("([Cc]ontrol|[Nn]ormal)", sampleTable$group)]
    }else if(is.character(normal_idx)){
        normal_idx=match(normal_idx ,sampleTable[,"sample_name"])
    }
    if(any(is.na(normal_idx)))
        stop("samples specified in normal_idx not found")
    message("using median of samples ",
            paste(sampleTable[normal_idx,"sample_name"], collapse=", "),
            " as CNV free reference")
    
    #genomic windows
    CNV_Regions=makeGenomeWindows(BSgenome, chr.select, window_size)
    dataset=getBSgenome(BSgenome, masked=FALSE)    
    CpGpos=GRanges(seqinfo=seqinfo(CNV_Regions))
    if(!is.null(CpGrange)){
        pattern="CG"
        for (chr in seqlevels(CNV_Regions)){
            message("searching ",chr," for \"",pattern,"\"...")
            chr_seq=dataset[[chr]]
            pIdx=start(matchPattern(pattern=DNAString(pattern),
                subject=chr_seq, fixed=TRUE))
            message("found ", length(pIdx), 
                " occurances of ", pattern , " in ", chr)
            CpGpos=c(CpGpos, GRanges(chr, IRanges(pIdx, width=1), 
                seqinfo=seqinfo(CNV_Regions)))
        }
    }
    
    if(parallel) {  
        BPPARAM=bpparam()
        message("Scanning ",bpnworkers(BPPARAM) , " files in parallel")
    }else
        BPPARAM=SerialParam()
    
    #get the read counts in windows
    counts=unlist(bplapply(X=file_names,FUN=getCoverage, Regions=CNV_Regions,
            paired=paired, uniquePos=uniquePos, fragment_length=fragment_length,
            CpGrange=CpGrange, CpGpos=CpGpos,
            BPPARAM=BPPARAM), FALSE, FALSE)
    
    counts=matrix(unlist(counts[seq(1,to=length(counts), by=2)],FALSE, FALSE), 
        ncol=length(file_names), byrow=FALSE)
    #library size normalization using median
    
    #zygosity=getZygosity(qs)
    olCy=match(as.character(seqnames(CNV_Regions)), colnames(zygosity))        


    counts[counts==0]=NA
    nf=apply(X=counts,MARGIN=2, FUN=median, na.rm=TRUE)
    lb=qnbinom(p=.001,mu=nf, size=100)
    ub=qnbinom(p=.999,mu=nf, size=100)
    nSmooth=5
    naM=matrix(NA, nSmooth, ncol(counts))
    ma=t(rbind(naM,counts,naM))
    smoothed=apply(ma, 1, rollapply,width=2*nSmooth+1, median, na.rm=TRUE)
    f=is.na(counts) | t(t(smoothed) < lb|t(smoothed)>ub)
    littleCNVFree=colSums(f)/nrow(f)>.8 
    if(any(littleCNVFree) ){
        f[,littleCNVFree]=is.na(counts[,littleCNVFree])
        warning("Few regions with \"normal\" coverage found for sample ",
            sampleTable[littleCNVFree,"sample_name"],
            ", leading to uncertainty in base coverage estimation.")
    }
    values=counts
    values[f]=NA
    nf=apply(X=values,MARGIN=2, FUN=median, na.rm=TRUE)
    #warn if medians are small
    if(any(nf<100, na.rm=TRUE) )
        warning("Low coverage in CNV files. Consider larger CNV_window_size")
    values=t(t(counts)/zygosity[,olCy]/nf)
    isna=!is.finite(values)|is.na(values)
    
    #CNV analysis
    medianValue=apply(X=values[,normal_idx, drop=FALSE], FUN=median, MARGIN=1)
    filter=medianValue>0
    #run CopyHMM on all samples and write 
    #CNV states (mean logFC of segments) in matrix    
    CNVval=matrix(NA, length(CNV_Regions),n,
        dimnames=list(seq_along(CNV_Regions), sampleTable$sample_name))
    lim=ifelse(is.null(mu),1.5,max(abs(mu))*1.1)
    for (i in seq_len(n)){
        message("== searching for CNVs in ",sampleTable$sample_name[i]," ==")

        CNV=as(CNV_Regions, "RangedData")        
        CNV$valid=filter & values[,i]>0 & medianValue>0
        CNV$copy=log2(values[,i]/medianValue)
        CNV$copy[!is.finite(CNV$copy)|is.na(CNV$copy)]=0
        as=rep(TRUE,length(space(CNV)))
        param=HMMcopy::HMMsegment(CNV, getparam=TRUE, autosomes=as)
        if(!is.null(mu)) param$m=mu
        CNV_seg=HMMcopy::HMMsegment(CNV, param,verbose=FALSE, autosomes=as)    
        if(! is.null(plot_dir)){
            for (chr in chr.select){
                png(paste(plot_dir,"/CNV_",sampleTable$sample_name[i],
                    "_",chr,".png",sep=""),width=600, height=300)
                HMMcopy::plotSegments(CNV, CNV_seg, pch=20, ylab="Copy Number",
                    xlab = "Chromosome Position", 
                    main=paste("CNV ",sampleTable$sample_name[i],chr),chr=chr)
                    #,ylim=c(-lim,lim)) 
                for(j in seq_len(nrow(CNV_seg$mus))) {
                    abline(h = CNV_seg$mus[j ,ncol(CNV_seg$mus)], 
                        col = HMMcopy::stateCols()[j], lwd = 2, lty = 3)
                }
                dev.off()
            }
        }
        CNV_seg$segs$median[CNV_seg$segs$state==3]=mu[3] #"common normal state"
        range=GRanges(seqnames=CNV_seg$segs$chr, 
            ranges=IRanges(start=CNV_seg$segs$start, end=CNV_seg$segs$end), 
            median=CNV_seg$segs$median)
        m=findOverlaps(CNV_Regions, range,select="first")        
        CNVval[,i]=range$median[m]
    }
    CNVval[isna]=NA
    mcols(CNV_Regions)<-CNVval
    names(mcols(CNV_Regions))=sampleTable$sample_name
    return(CNV_Regions)
}


addCNV<-function(qs,file_name, window_size=1000000, paired=FALSE, 
        fragment_length=NULL,cnv=NULL, mu =log2(c(1/2, 2/3, 1, 3/2,2,3)), 
        normal_idx=NULL, plot_dir=NULL, MeDIP=FALSE, parallel=FALSE){    
    if(! missing(cnv) ){
        if(class(cnv)!="GRanges")
            stop("please provide CNVs as GRange object")
        return(setCNV(qs,cnv))
    }#else find CNVs:
    qs=setCNV(qs,findCNV(sampleTable=getSampleTable(qs),
        BSgenome=getGenome(qs),file_name=file_name, chr.select=getChrNames(qs),
        window_size=window_size, paired=paired, fragment_length=fragment_length,
        mu=mu, normal_idx=normal_idx, plot_dir=plot_dir, MeDIP=MeDIP, 
        zygosity=getZygosity(qs), parallel=parallel)) 
    #todo: add option to use regions of interest, not bsgenome    
    if(length(getOffset(qs))>0 && ! all(is.na(getOffset(qs) )) )
        warning("Consider recalculating offset based on new CNV values")
    if("seqPref" %in% names(mcols(getRegions(qs))))
        warning("Consider recalculating sequence ",
                "preference based on new CNV values")
    return(qs)
}

#taken from HMMcopy, added "ylim" parameter
#plotSegments<-function (correctOutput, segmentOutput, 
#        chr = space(correctOutput)[1],...) 
#{
#    if (is.null(segmentOutput$segs)) {
#        warning("Processed segments now found, automatically processing")
#        segmentOutput$segs <- HMMcopy:::processSegments(segments$segs, 
#            space(correctOutput), start(correctOutput), end(correctOutput), 
#            correctOutput$copy)
#    }
#    segs <- segmentOutput$segs
#    correctOutput$state <- segmentOutput$state
#    cols <- HMMcopy::stateCols()
#    if (missing(ylim)){
#        ylim <- quantile(correctOutput$copy, na.rm =TRUE, prob=c(0.01, 0.99))
#    }
#    a <- correctOutput[as.character(chr)]
#    b <- segs[segs$chr == chr, ]
#    plot(start(a), a$copy, col = cols[as.numeric(as.character(a$state))], 
#        ylim = ylim, ...)
#    for (k in seq_len(nrow(b))) {
#        lines(c(b$start[k], b$end[k]), rep(b$median[k], 2), lwd = 3, 
#            col = "green")
#    }
#}


