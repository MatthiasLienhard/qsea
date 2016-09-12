
#internal function to return normalized values
#    return selected samples and regions (idx)


getNormalizedValues<-function(qs, methods, windows=NULL, samples=NULL,
        groupMeans=NULL, verbose=FALSE, minEnrichment=0, chunksize=1e5){
    
    if(is.null(windows))
        windows=seq_along(getRegions(qs))
    if(is.null(samples) && is.null(groupMeans)){
        samples=getSampleNames(qs)
    }else if(is.numeric(samples))
        samples=getSampleNames(qs, samples)
    if(!is.null(groupMeans)){
        if(class(groupMeans) != "list")
            stop("groupMeans must be a list")
        for(i in seq_along(groupMeans) ){
            if(is.numeric(groupMeans[[i]]))
                groupMeans[[i]]=getSampleNames(qs, groupMeans[[i]])
        }
        allSampleIdx=unique(c(samples,unlist(groupMeans)))
    }else{
        allSampleIdx=samples
    }
    nchunk=ceiling(length(windows)/chunksize)
    if(nchunk>1){
        windows=split(windows, 
            cut(seq_along(windows), nchunk, labels = FALSE)) 
        if(verbose)
            message("normalising values in ",nchunk," cunks")
    }
    getChunk<-function(wd){
    if(verbose)
        message("obtaining raw values for ", length(allSampleIdx), 
            " samples in ",length(wd)," windows")
    
    raw_values=getCounts(qs, windows=wd, samples=allSampleIdx)

    #val_table=matrix(NA, length(wd),length(samples)*length(methods),
    #    dimnames=list(NULL, 
    #        paste(samples,rep(names(methods), each=length(samples))
    #        ,sep="_"))
    #    )
    #if(!is.null(groupMeans)){
    #    mean_table=matrix(NA, length(wd),
    #        length(groupMeans)*length(methods),
    #        dimnames=list(NULL, 
    #            paste(names(groupMeans),
    #            rep(names(methods), each=length(groupMeans)),
    #            "means",sep="_"))
    #        )
    #}else{mean_table=NULL}
    #"enrichment", "cnv", "library_size","region_length", "preference"
    #nomalization factor matrix in logscale
    getVal<-function(m){
    #for(m in names(methods)){
        if(verbose)
            message("deriving ",m," values...")
        pseudocount=getNumPar("psC",methods[[m]])    
        if(is.null(pseudocount)) pseudocount=0
        minVal=getNumPar("minCut",methods[[m]])    
        maxVal=getNumPar("maxCut",methods[[m]])    
        qPer=getNumPar("q",methods[[m]])
        values=raw_values+pseudocount
        names(values)=paste(samples,m,sep="_")
        nm=getNormMatrix(qs, methods[[m]],wd,allSampleIdx)
        norm=nm$factors
        offset=nm$offset
        rm(nm)
        #offset=matrix(0, nrow(values),ncol(values))        
        #if("offset" %in% methods[[m]]){
        #    offset=t(t(norm)*getOffset(qs, allSampleIdx))
        #}        

        #transformation to absolute methylation
        if("enrichment" %in% methods[[m]]){
            pattern_name=getEnrichmentPattern(qs)
            if( is.null(pattern_name) || 
                !paste0(pattern_name, "_density") %in% 
                    names(mcols(getRegions(qs))) || 
                ! hasEnrichment(qs))
                    stop("missing parameters for beta values")
                #should not happen at this poin
            cf=mcols(getRegions(qs))[wd,paste0(pattern_name, "_density")]
            sPar=getEnrichmentParameters(qs)
            
            offset=t(offset*t(norm))
            #for(i in seq_along(samples)){
            #    normFM[,i]=normFM[,i]*
            #        .sigmF2(cf,sPar[samples[i],"a"],
            #            sPar[samples[i],"b"],sPar[samples[i],"c"])+pseudocount
            #}
            norm=norm*mapply(.sigmF2, MoreArgs=list(x=cf),
                a=sPar[samples,"a"],
                b=sPar[samples,"b"],
                c=sPar[samples,"c"], USE.NAMES = FALSE)+pseudocount

            if(minEnrichment>0)
                norm[norm<minEnrichment]=NA 
            else
                norm[norm<=0]=NA 

            #for(i in seq_along(allSampleIdx)){
            big=values>5*norm+offset 
                #prevents NaN for big y (dgamma gets 0)
            if(is.null(qPer)){
                values=.estimateMethPois( y=values, c=norm,o=offset)
            }else{
                dimV=dim(values)
                dimN=dimnames(values)
                values=.estimateQuantilePois(y=as.numeric(values),
                        c=as.numeric(norm),o=as.numeric(offset), p=qPer/100)
                dim(values)=dimV
                dimnames(values)=dimN
            }
            values[big]=1
            rm(big)            
        }else{#scale by normalization factor
            if(any(offset!=0))            
                values=t(t(values/norm)-offset)
            else
                values=values/norm
        }
        rm(offset)
        rm(norm)
        if("logit" %in% methods[[m]]){
            values=log2(values/(1-values))
        }else if ("log" %in% methods[[m]]){
            values=log2(values)
        }
        if(!is.null(maxVal)){
            values[values>maxVal]=maxVal
        }
        if(!is.null(minVal)){
            values[values<minVal]=minVal
        }
        #if(!is.null(samples)){
        #    val_table[,paste(samples,m,sep="_")]=
        #        values[,match(samples,allSampleIdx)]
        #}
        #if(!is.null(groupMeans)){
        #    for(i in seq_along(groupMeans))
        #        mean_table[,paste(names(groupMeans)[i],m,"means",sep="_")]=
        #            rowMeans(values[,
        #                which(allSampleIdx %in% groupMeans[[i]]),drop=FALSE], 
        #            na.rm=TRUE)
        #}
        mean_table=NULL
        if(!is.null(groupMeans))
            mean_table=vapply(groupMeans,function(grp){
                rowMeans(values[,which(allSampleIdx %in% grp),drop=FALSE], 
                    na.rm=TRUE)
                }, numeric(length(wd)))
            
        return(list(values[,samples], mean_table))
                
    }#end getVal

    idx=seq(from=1,by=2,length.out=length(methods))
    return(do.call(what=cbind, 
        args=unlist(lapply(names(methods),getVal)
            ,FALSE, FALSE)[c(idx, idx+1)]))
    }#end getChunk
    if(nchunk==1)
        ret=getChunk(windows)
    else 
        ret=rbind(do.call(what=rbind,args=lapply(windows, getChunk)))

    retnames=character(0)
    if(!is.null(samples))
        retnames=paste(samples,
            rep(names(methods), each=length(samples)),sep="_")
    if(!is.null(groupMeans))    
        retnames=c(retnames,paste(names(groupMeans),
            rep(names(methods), each=length(groupMeans)),"means",sep="_"))
    colnames(ret)=retnames
    return(ret)
}


getNumPar<-function(string, methods){
    idx=grep(paste0("^",string,"\\d+(\\.\\d+)?$"), methods)
    if(length(idx)>1)
        warning("more than one value found for ",string)
    if(length(idx)>0)
        return(as.numeric(sub(string,"",methods[idx[1]])))
}


getNormMatrix<-function(qs, methods,windows,samples){

    if(class(methods) == "normMethod") {
        methods=methods[[1]]
        warning("selected first normMethod") #this should not happen
    }
    regInfo=names(mcols(getRegions(qs)))
    if("library_size" %in% methods) # per million reads
        normFM=matrix(getLibSize(qs, samples)/1e6 ,length(windows), 
            length(samples), byrow=TRUE)
    else
        normFM=matrix(1,length(windows), length(samples))
    if("region_length" %in% methods)#per kilobase
        normFM=normFM*getWindowSize(qs)*1e-3    
        #todo: allow different window sizes? change getWindowSize
    if("preference" %in% methods && 
        "seqPref" %in% regInfo)
        normFM=normFM *( 2 ^ getRegions(qs)$seqPref[windows])
    if("cnv" %in% methods && hasCNV(qs) ){
        om=as.matrix(findOverlaps( getRegions(qs)[windows], getCNV(qs) ))
        om=om[!duplicated(om[,1]),, drop=FALSE]
        cnv_os=as.matrix(mcols(getCNV(qs)))[om[,2],samples, drop=FALSE]
        normFM[om[,1],]=normFM[om[,1],]*(2 ^ cnv_os)
        
        #rm(om, cnv_os)
    }
    if("zygosity" %in% methods ){
        m=match(as.vector(seqnames(getRegions(qs)[windows])), 
            colnames(getZygosity(qs)))
        zygos_os=t(getZygosity(qs))[m,samples, drop=FALSE]/2
        normFM=normFM*zygos_os
        #rm(m, zygos_os)
    }
    normFM[normFM<=0]=NA

   if("offset" %in% methods){#this is quite inefficient
        #offset=t(getOffset(qs, samples)*t(normFM))
        offset=getOffset(qs, samples)
    }else{
        #offset=sparseMatrix(x=0,i=1,j=1, w=length(windows),ncol=length(samples))
        offset=rep(0,length(samples))
    }
    return(list(factors=normFM, offset=offset))
}





.sigmF=function(x) x/sqrt(1+x^2)

.sigmF2=function(x,a=0,b=1, c=1,o=0)
    (.sigmF(x/c-a)-.sigmF(-a))*b/(1-.sigmF(-a))+o
#a= -> <- shift, b=f(INF) y scale, c= <--> x scale
.estimateMethPois=function(y,c,o) {
    n=length(y)
    if(length(c)!=n ||length(o)!=n){
        stop("argument length missmatch")
    }
    return(( (y+1)*(pgamma(c+o,y+2) - pgamma(o,y+2))+ 
        o* (pgamma(o,y+1) - pgamma(c+o,y+1) ) )/
        (c*(pgamma(o+c, y+1) - pgamma(o,y+1))) )
}

.delta=function(y,c,o,p,x) 
    (pgamma(o+c*x,y+1)- pgamma(o,y+1))/(pgamma(o+c,y+1)- pgamma(o,y+1)) -p

#use binary search to find the inverse
.estimateMethPoisInv=function(beta,c,o,tol=.01, eps=1e-16, nIter=20) {
    n=length(beta)
    n=length(y)
    if(length(c)!=n ||length(o)!=n){
        stop("argument length missmatch")
    }
    #c=rep(c, length.out=n)
    #o=rep(o, length.out=n)
    todo=which(!(is.na(beta) | is.na(c) | is.na(o) | c==0))
    y=matrix(c(rep(0,n),beta*c+o),n,2)    
    lessThan0=which(.estimateMethPois(y[todo,1],c[todo],o[todo])>beta[todo])
    y[todo[lessThan0],2]=0
    todo=todo[-lessThan0]
    boundToSmall=todo[which(
        .estimateMethPois(y[todo,2],c[todo],o[todo])<beta[todo])]
    while(length(boundToSmall)>0){
        y[boundToSmall,1]=y[boundToSmall,2]
        y[boundToSmall,2]=2*y[boundToSmall,2]
        boundToSmall=boundToSmall[
        .estimateMethPois(y[boundToSmall,2],c[boundToSmall],o[boundToSmall])<
        beta[boundToSmall] 
        & y[boundToSmall,2]<5
        *c[boundToSmall]+o[boundToSmall]]
    }
    moreThanX=which(y[todo,2]>5*c[todo]+o[todo])
    y[todo[moreThanX],]=5*c[todo[moreThanX]]+o[todo[moreThanX]]
    todo=todo[-moreThanX]
    iter=0
    while(length(todo)>0 && iter<nIter){
        iter=iter+1
        y_new=(y[todo,1]+y[todo,2])/2
        beta_new=.estimateMethPois(y_new,c[todo],o[todo])
        tomuch=(beta_new >= beta[todo])
        y[todo[ tomuch],2]=y_new[ tomuch]
        y[todo[!tomuch],1]=y_new[!tomuch]
        todo=todo[y[todo,2]-y[todo,1]>tol]
    }
    return((y[,1]+y[,2])/2)

}


#uses binary search to find the quantile
.estimateQuantilePois<-function(y,c,o, p=.5,tol=.001, eps=1e-16 ){
        niter=ceiling(log2(1/tol))
    n=length(y)
    if(length(c)!=n ||length(o)!=n){
        stop("argument length missmatch")
    }
    #c=rep(c, length.out=n)
    #o=rep(o, length.out=n)
    q=idx=rep(-1,n)
    x=delta_val=matrix(c(0,1),n,2, byrow=TRUE)#upper and lower bounds
    delta_val=matrix(c(.delta(y,c,o,p,x[,1]),.delta(y,c,o,p,x[,2])),n,2)
    todo=seq_len(n)
    inNA=which(is.na(y) | is.na(c) | is.na(o)|c==0 )
    overexp=which(pgamma(o+c, y+1)<eps | delta_val[,2] <= 0 ) 
        #^way more observed reads (y) than expected (c+o) -->methylation=1
    underexp=which(delta_val[,1]>=0) 
        #^less observed reads than offset -->methylation=0
    if(length(c(inNA , overexp, underexp))>0)
        todo=todo[-c(inNA, overexp, underexp)]
    for(i in seq_len(niter)){
        xnew=(x[todo,1]+x[todo,2])/2
        delta_new=.delta(y[todo],c[todo],o[todo],p,xnew)
        tomuch=(delta_new<=0)
        x[todo[ tomuch],1]=xnew[ tomuch]
        x[todo[!tomuch],2]=xnew[!tomuch]
        delta_val[todo[ tomuch],1]=delta_new[ tomuch]
        delta_val[todo[!tomuch],2]=delta_new[!tomuch]
    }
    q[inNA]=NA
    q[overexp]=1
    q[underexp]=0
    r=todo
    w=delta_val[r,2]/(delta_val[r,2]-delta_val[r,1]) #weight
    q[r]=w*x[r,1]+(1-w)*x[r,2] #weighted average of upper and lower estimate
    return(q)
}


estimateEnrichmentLM<-function(qs,windowIdx, signal, min_wd=5,
    bins=seq(.5,40.5,2), trim=.05, pattern_name){
    #check: offset and pattern density must be set    
    if(missing(pattern_name))
        pattern_name=getEnrichmentPattern(qs)
    if(is.null(pattern_name))
        stop("please specify sequence pattern name")
    if(! paste0(pattern_name,"_density") %in% names(mcols(getRegions(qs))))
        stop("no ",pattern_name,
            " density found. Please run \"addPatternDensity\" first")
    if(missing(windowIdx))
        windowIdx=seq_along(getRegions(qs))
    m=length(windowIdx)
    
    samples=getSampleNames(qs)
    n=length(samples)
    if(missing (signal)){
            stop("plese provide a calibration signal or trimQuantile")
    }
    if (class(signal)=="numeric") {
        if(length(signal)==1)
            type="allSame"
        else if(length(signal)!=m)
            stop("number of windows (",m, ") does not match number of values (",
                length(signal),") in \"signal\"")
        else type="onePerWd"
    }else if (class(signal)=="matrix") {
        if(n != ncol(signal) )
            stop("number of samples (",n,") does not match number of columns (",
                ncol(signal),") in \"signal\"")
        if(m != nrow(signal) )
            stop("number of windows (",m,") does not match number of rows (",
                nrow(signal),") in \"signal\"")
            type="individual"
    }else{stop("signal must be provided as numeric vector or matrix, but is ",
        class(signal))}
    patternD=mcols(getRegions(qs))[windowIdx, paste0(pattern_name,"_density")]
        
    nbin=length(bins)-1
    binWd=numeric(length(getRegions(qs)))
    cf=matrix(NA,nbin,n, dimnames=list(NULL,samples))
    #for each density class, estimate cf from provided signal
    bin_median=bin_n=numeric(nbin)
    norm_method=normMethod("beta") 
    norm_method$beta=norm_method$beta[-grep("enrichment", norm_method$beta)]
    vals=getNormalizedValues(qs,methods=norm_method, windows=windowIdx, 
        samples=samples)
    getSignal<-function(signal,wd,sa,type){
        switch(type,
            allSame=signal,
            onePerWd=signal[wd],
            individual=signal[wd,sa])}

    for(i in seq_len(nbin)){
        dcw=which(patternD>=bins[i] & patternD < bins[i+1])
        bin_median[i]=median(patternD[dcw])
        bin_n[i]=length(dcw)

        if(length(dcw)>=min_wd*(1-2*trim)){
            for(j in seq_len(n))
                cf[i,j]=mean(vals[dcw,j]/
                    getSignal(signal,dcw,j,type), na.rm=TRUE, trim=trim)

        }
    }    
    return(list(factors=cf, density=bin_median, n=bin_n, 
            pattern_name=pattern_name))

}


fitEnrichmentProfile<-function(factors, density, n, minN=1,...){
    #fit the sigmoidal function
    RSSfun=function(par,pd,factor, w=1) 
        sum((.sigmF2(x=pd,par[1], par[2],par[3])-factor)^2*w)
    #weight w ~ 1/SEM
    #SEM= sd/sqrt(n)
    #assuming sd is independent of pattern density pd --> w=sqrt(n)
    if(class(factors)=="numeric")
        factors=matrix(factors)
    sigmPar=matrix(NA,ncol(factors),3, 
        dimnames=list(colnames(factors), c("a", "b", "c")))
    upperQ=apply(X=factors, MARGIN=2, FUN=quantile,p=.75, na.rm=TRUE)
    for(j in seq_len(ncol(factors))){
        CF_used=which(!is.na(factors[,j]) &is.finite(factors[,j]) & n>=minN)
        #for(i=5:length(CF_used))
        sigmPar[j,]=optim(par=c(1, upperQ[j],10), fn=RSSfun, 
            factor=factors[CF_used,j],pd=density[CF_used],w=sqrt(n[CF_used]), 
            lower=c(-1, upperQ[j]/2,1), upper=c(3, upperQ[j]*2,20), 
            method="L-BFGS-B")$par
    }
    return(sigmPar)
}


estimateOffset<-function(qs,enrichmentPattern, maxPatternDensity=0.01){
    
    if(! paste0(enrichmentPattern,"_density") %in% 
        names(mcols(getRegions(qs)))){
        stop("no ",enrichmentPattern,
            " density found. Please run \'addPatternDensity\' with name=\"",
            enrichmentPattern,"\" first.")
    }
    message("selecting windows with low ",enrichmentPattern,
        " density for background read estimation")
    
    patternD=mcols(getRegions(qs))[, paste0(enrichmentPattern,"_density")]

    bg=which(patternD <= maxPatternDensity)
    nna=!is.na(patternD)
    nValid=sum(nna)
    #sum(apply(FUN=any, X=!is.na(mcols(getRegions(qs))) , MARGIN=1 ))
    fraction=length(bg)/ nValid
    if(length(bg)<100) 
        stop("not enough windows with enrichment pattern density of at most ", 
        maxPatternDensity, 
        " per fragment found to estimate background read distribution")
    message(round(fraction*100,3),
        "% of the windows have enrichment pattern density of at most ", 
        maxPatternDensity, 
        " per fragment and are used for background reads estimation")
    bgCounts=getNormalizedValues(qs,methods=normMethod("nrpkm"), 
        windows=bg, samples=getSampleNames(qs))
    #upperQ=apply(bgCounts, 2,quantile,p=.5,.75) #remove stangly high values?
    #limit=upperQ[1]+5*(upperQ[2]-upperQ[1])
    offset=colMeans(bgCounts, na.rm=TRUE)
    return(offset)
}

findSeqPref<-function(qs,file_name, fragment_length, paired, 
        uniquePos, alpha, pseudocount,minMapQual, cut,
        parallel=FALSE){
    samples=getSampleTable(qs)
    Regions=getRegions(qs)
    
    fname_idx=which(names(samples)==file_name )
    if(length(fname_idx)!=1) stop("no column ", file_name, 
        " found in sample table")
    #read files (not all need to be present)
    sIdx=which(!is.na(samples[,fname_idx]) )
    if(length(sIdx)==0) stop("no files found in column ", file_name)
    
    message("estimating sequence preference from ",length(sIdx),
        " files in \"",file_name, "\" column...")
     if(parallel)
        BPPARAM=bpparam()
    else
        BPPARAM=SerialParam()

    
    #get the read counts in windows
    coverage=unlist(bplapply(X=samples[sIdx,fname_idx],FUN=getCoverage, 
            Regions=Regions,fragment_length=fragment_length,
            minMapQual=minMapQual, paired=paired, uniquePos=uniquePos, 
            BPPARAM=BPPARAM ), FALSE, FALSE)
    
    libraries=matrix(unlist(coverage[seq(2,length(coverage),by=2)],FALSE,FALSE),
        nrow(samples), 6,byrow=TRUE, dimnames=list(samples$sample_name, 
            c("total_fragments", "valid_fragments","library_factor", 
                "fragment_length", "fragment_sd", "offset")))

    coverage=matrix(unlist(coverage[seq(1,length(coverage),by=2)],FALSE,FALSE), 
        ncol=nrow(samples), byrow=FALSE, 
        dimnames=list(NULL, samples$sample_name))

    message("estimating sequence preference...")
    ws=width(Regions)[1]
    if(hasCNV(qs)){
        #correct for CNV
        CNV=getCNV(qs)
        for (i in seq_along(sIdx)){
            altered=which(mcols(CNV)[,sIdx[i]]!=0)
            if(length(altered)>0){
                ol=findOverlaps(Regions, CNV[altered], 
                    minoverlap=ceiling(ws/2)+1,select="first")
                altered_wd=which(!is.na(ol))
                if(length(altered_wd)>0){
                    coverage[altered_wd,i]=coverage[altered_wd,i]/
                        2^(mcols(CNV)[altered[ol[altered_wd]],sIdx[i]])
                }
            }            
        }
    }
    els=colSums(coverage)
    libraries$effective=els
    mls=median(libraries[,"valid_fragments"])
    
    seqPref=rowSums(coverage)
    #seqPref=colSums(t(coverage)/els*mls) #scale to median library size
    ci=qpois(c(.01,.99),mean(seqPref))
    #cut the tails to calculate mean of distribution
    dm=mean(seqPref[seqPref>ci[1] & seqPref<ci[2] ]) 
    logSeqPref=log2((seqPref+pseudocount)/(dm+pseudocount) )    
    ci=qpois(c(.025,.975),dm)
    logSeqPref[seqPref>ci[1] & seqPref<ci[2] ]=0
    logSeqPref[abs(logSeqPref)>cut]=NA    
    return(list(seqPref=logSeqPref, libraries=libraries))
    #compute sum y
    #estimate E(y)
    #compute p(y~Pois(E(y))
    #if p-value<alpha seqPref= log2(y)-log2(E(y))
}

    
