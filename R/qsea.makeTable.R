makeTable<-function(qs,glm,norm_methods="counts",samples,groupMeans, keep, ROIs,
    annotation, minPvalSummarize, CNV=FALSE, verbose=TRUE, minEnrichment=3, 
    chunksize=1e5){
     if(missing(qs) || class(qs) != "qseaSet" )
        stop("please specify a qseaSet\n")    
    if(missing(samples))
        samples=NULL
    samples=checkSamples(qs, samples)
    if(missing(groupMeans))
        groupMeans=NULL
    else{
        if(class(groupMeans) != "list")
            stop("\"groupMeans\" must be a list of sample ids")
        groupMeans=lapply(groupMeans, checkSamples, qs=qs)
    }

    if(CNV && hasCNV(qs)==FALSE ) 
        stop("qs does not contain CNV information. Run addCNV first")
    #check norm methods
    if(class(norm_methods)=="character")
        norm_methods=normMethod(norm_methods)
    if(class(norm_methods)!="normMethods") 
        stop("Invalid norm_methods. Please see normMethod() ",
            "function for predefined and user defined normalization methods.")
    pattern_name=getEnrichmentPattern(qs)

    if( "enrichment" %in% unlist(norm_methods) && ( is.null(pattern_name) || 
        (! paste0(pattern_name,"_density") %in% names(mcols(getRegions(qs))))||
        !hasEnrichment(qs)))
        stop("Transformation to absolute methylation requires calibration of ",
            "enrichment profiles. Please run addEnrichmentParameters() first.")
    if(! missing(glm) && class(glm) != "qseaGLM")
        stop("GLM object must be of type \"qseaGLM\"")
    if(missing(keep))
        keep=seq_along(getRegions(qs))
    if(is.logical(keep)){
        if(length(keep)==1) keep=seq_along(getRegions(qs))
        else keep=which(keep)
    }
    if(length(keep)==0)
        stop("length(keep)==0: no regions selected")
    if(!missing(annotation) && (class(annotation)!="list" || 
        any(sapply(annotation,class)!="GRanges")||
        is.null(names(annotation))))
            stop("annotation must be a named list of GRanges objects")
    

    if(!missing(ROIs) && class(ROIs) =="GRanges" ){
        reducedM=NULL
        if(! missing(minPvalSummarize))
            reducedM=getReducedModel(glm, minPvalSummarize)
        if(!is.null(reducedM)){
            message("selecting the window with minimal p-vaule in ",
                minPvalSummarize," per ROI")
            glm_wd=getWindows(glm)
            keep=intersect(keep, glm_wd[!is.na(reducedM$LRT_pval)])
            if(length(keep)==0)
                stop("no regions left")
            m=IRanges::as.data.frame(findOverlaps(getRegions(qs)[keep], ROIs))
            if(ncol(mcols(ROIs))>1)
                roi_id=apply(X=as.data.frame(mcols(ROIs)), FUN=paste, 
                    MARGIN=1,collapse="%_%")
            else if(ncol(mcols(ROIs))==1)
                roi_id=as.character(mcols(ROIs)[,1])
            else roi_id=seq_len(mcols(ROIs))
            roi_list=split(x=m$queryHits,f=roi_id[m$subjectHits]) 
                #for each gene, the indices of windows
            nonEmpty = sapply(X = roi_list, FUN = length) > 0 
            p_idx=match(keep,as.numeric(glm_wd))    
            sel= unlist(sapply(X = roi_list[nonEmpty], FUN = function(x) {
                x[which.min(
                    reducedM$LRT_pval[p_idx[x]])]
            }))
            keep=keep[sel]
            ord=order(keep)#sort roi_tab wrt window number
            keep=keep[ord]
            roi_idx=match(names(roi_list[nonEmpty]),roi_id)
            if(ncol(mcols(ROIs))>0){
                roi_tab=as.data.frame(mcols(ROIs))[roi_idx, ,drop=FALSE]
            }else{
                roi_tab=as.data.frame(ROIs)[roi_idx,c(2,3)]
            }
            roi_tab=roi_tab[ord,, drop=FALSE]
        }else{
            message("selecting specified ROIs")
            m=IRanges::as.data.frame(findOverlaps(getRegions(qs)[keep], ROIs))
            roi_tab=as.data.frame(ROIs)[m[,2],, drop=FALSE]
            keep=keep[m[,1]]        
            if(length(keep)==0)
                stop("no regions left")
            roi_tab=roi_tab[,c(-1,-4), drop=FALSE] #remove chr and width
        }
        n=length(keep)
        names(roi_tab)=paste("ROI",names(roi_tab),sep="_")
    }else{
        n=length(keep)
        roi_tab=as.data.frame(matrix(NA,n,0))
    }
    
    if(n==0) stop("no windows within ROIs...")
    #message(n," windows selected")
    #test_tab=table containing the test results (p-value, logFC)
    if(!missing(glm)){
        test_tab=as.data.frame(matrix(NA,n,3*length(getContrastNames(glm))))
        names(test_tab)=paste(rep(getContrastNames(glm), each=3), 
            c("log2FC", "pvalue", "adjPval"), sep="_")
        idx=match(keep, getWindows(glm))    
        for(i in getContrastNames(glm) ){
            message("adding test results from ", i)
            #add effect
            contr=getReducedModel(glm, i)
            test_tab[,paste(i, "log2FC", sep="_")]=contr$effect[idx]            
            ##add pvalue
            test_tab[,paste(i, "pvalue", sep="_")]=contr$LRT_pval[idx]
            ##add adjusted pvalue    
            test_tab[,paste(i,"adjPval", sep="_")]=contr$LRT_adjP[idx]
        }
    }else
        test_tab=as.data.frame(matrix(NA,n,0))
    
    #count_tab contains the sample values, normalized by norm_method
    if(length(groupMeans)>0 || length(samples)>0)
        count_tab=getNormalizedValues(qs, methods=norm_methods ,windows=keep,
            samples=samples,groupMeans=groupMeans,verbose=verbose, 
            minEnrichment=minEnrichment, chunksize=chunksize)
    else
        count_tab=as.data.frame(matrix(NA,n,0))
    cnv_tab=as.data.frame(matrix(NA,n,0))
    if(CNV && ! missing(samples)){
        cnv_tab=matrix(NA,n,length(samples), dimnames=
            list(NULL,sub("CNV","CNVlogFC",names(mcols(getCNV(qs)))[samples])))
        message("adding CNV logFC values")
        om=as.matrix(findOverlaps( getRegions(qs)[keep], getCNV(qs) ))
        om=om[!duplicated(om[,1]),]
        cnv_tab[om[,1],]=as.matrix(mcols(getCNV(qs)))[om[,2],samples]
        colnames(cnv_tab)=paste0(colnames(cnv_tab), "_CNV")
    }
    if(!missing(annotation)){
        if(verbose)
            message("adding annotation")
        anno_tab=as.data.frame(matrix("",n,length(annotation), 
            dimnames=list(NULL,names(annotation))), stringsAsFactors=FALSE)
        for(reg in names(annotation)){
            ol=as.data.frame(
                findOverlaps(getRegions(qs)[keep], annotation[[reg]]))
            m=ncol(mcols( annotation[[reg]][1]))
            if(m==0)
                ol$anno_names=rep("yes",nrow(ol) )
            else if(m==1)
                ol$anno_names=as.character(
                    mcols(annotation[[reg]][ol[, 2]])[,,drop=TRUE])
            else
                ol$anno_names=do.call(mapply,
                    c(as.list(mcols(annotation[[reg]][ol[,2]])),FUN=paste))
                #ol$anno_names=apply(
                #   X=as.data.frame(mcols(annotation[[reg]][ol[,2]])),
                #   MARGIN=1,FUN=paste,collapse="_")


            ol=ol[!duplicated(ol[,c(1,3)]),]
            anno_names=sapply(split(x=ol$anno_names, f=ol$queryHits),
                paste,collapse=", ")
            anno_tab[as.integer(names(anno_names)),reg]=anno_names
        }
    }else
        anno_tab=as.data.frame(matrix(NA,n,0))

    if(verbose)        
        message("creating table...")
    #remove width and strand
    reg_tab=as.data.frame(getRegions(qs))[keep,-c(4:5)]
    colnames(reg_tab)[1:3]=c("chr", "window_start", "window_end")
    result_tab=cbind(reg_tab,roi_tab, anno_tab,test_tab, count_tab,cnv_tab)
    rownames(result_tab)=seq_len(nrow(result_tab))
    return(result_tab)
}

#isSignificant
#returns a index vector with the regions that are different regarding the 
#specified criteria. 
#the resulting vector can then be passed to the keep parameter in makeTable, 
#to create a table of significant regions.
#idea: maybe include ROIs chr usw and rename function to selectRegions?

isSignificant<-function(glm,contrast=NULL, fdr_th=NULL, pval_th=NULL,
        absLogFC_th=NULL, direction="both"){
    if(class(glm) != "qseaGLM") stop("please specify a qseaGLM object")
    if(is.null(fdr_th) && is.null(pval_th) &&is.null(absLogFC_th) ){
        message("No thresholds defined. Selecting regions with fdr < 0.1")
        fdr_th=0.1
    }
    if(length(getContrastNames(glm))==0){
        stop("no contrasts found. Run addContrast() first.")
    }
    if (missing (contrast))    
        contrast=getContrastNames(glm)[1]
    if(is.null(getReducedModel(glm, contrast)))
        stop("cannot find specified contrast ",contrast)
    if(is.numeric(contrast))    
        contrast=getContrastNames(glm)[contrast]
    sig=rep(TRUE, length(getWindows(glm)))
    message("selecting regions from contrast ", contrast)

    direction=match.arg(direction, 
        c("both","hypo", "hyper","more", "less", "gain", "loss"))
    reducedModel=getReducedModel(glm, contrast)
    if(! is.null(fdr_th))
        sig[is.na( reducedModel$LRT_pval)|reducedModel$LRT_adjP>fdr_th]=FALSE
    if(! is.null(absLogFC_th))
        sig[is.na(reducedModel$effect) | 
            abs(reducedModel$effect)<absLogFC_th]=FALSE
    if(! is.null(pval_th))
        sig[is.na( reducedModel$LRT_pval) | 
            reducedModel$LRT_pval>pval_th]=FALSE    
    if(direction %in% c("hypo", "less", "loss"))
        sig[is.na(reducedModel$effect) | reducedModel$effect>0]=FALSE
    if(direction %in% c("hyper", "more", "gain"))
        sig[is.na(reducedModel$effect) | reducedModel$effect<0]=FALSE
    sig[is.na(sig)]=FALSE
    message("selected ", sum(sig), "/",length(sig)," windows")
    return(getWindows(glm)[sig])
}


regionStats<-function(qs,
    subsets=list(covered=which(rowSums(getCounts(qs))>=20)), ROIs=list(),
    minoverlap=1L, maxgap=-1L )
{
    if(missing(qs) || class(qs)!="qseaSet")
        stop("qs must be a qseaSet object")
    if(class(ROIs)=="GRanges")
        ROIs=list(ROI=ROIs)
    if(class(ROIs)!="list")
        stop("ROIs must be a (list of) Granges object(s)")
    if(class(subsets)=="numeric")
        subsets=list(subset=subsets)
    if(class(subsets)!="list")
        stop("subsets must be a list of index sets")
    if(is.null(names(ROIs))| is.null(names(subsets)) )
        stop("lists must be named")
    rn=c("all Regions", names(ROIs))
    cn=c("total", names(subsets))
    stats=matrix(NA, nrow=length(ROIs)+1, ncol=length(subsets)+1 , 
        dimnames=list(rn, cn) )
    reg=getRegions(qs)
    stats[1,]=c(length(reg), sapply(subsets, length))
    message("getting numbers of overlaps with total qs")
    stats[-1,1]=sapply(X=ROIs,
        FUN=function(x,query,...){sum(overlapsAny(query,x,...))},
        query=reg , minoverlap =minoverlap,maxgap=maxgap)
    for (i in names(subsets)){
        message("getting numbers of overlaps with ",i)
        reg=getRegions(qs)[subsets[[i]]]
        stats[-1,i]=sapply(X=ROIs,
            FUN=function(x,query,...){sum(overlapsAny(query,x,...))},
            query=reg, minoverlap =minoverlap,maxgap=maxgap)
    }
    return(stats)
}

