
#######
# PCA #
#######

getPCA<-function(qs, chr=getChrNames(qs),ROIs, minRowSum=20, keep,norm_method=
        normMethod(logRPM=c("log", "library_size","cnv","preference","psC10")),
        topVar=1000, samples=seq_len(nrow(getSampleTable(qs)))){
    if(missing(keep)) 
        keep=which(rowSums(getCounts(qs)) >= minRowSum )
    else
        keep=intersect(keep, which(rowSums(getCounts(qs)) >= minRowSum ))
    if(class(norm_method)=="character"){
        norm_method=normMethod(norm_method)
    }
    keep=keep[as.character(seqnames(getRegions(qs)[keep])) %in% chr]
    if(is.null(samples))
        samples=    getSampleNames(qs)
    else if (is.numeric(samples))
        samples=getSampleNames(qs, samples)
    if(! missing (ROIs) ){
        if(ncol(as.data.frame(mcols(ROIs)))==0){
            mcols(ROIs)=paste("ROI", seq_along(ROIs))
        }
        m=findOverlaps(query=getRegions(qs)[keep],subject=ROIs, select="first")
        keep=keep[!is.na(m)]
        if(ncol(as.data.frame(mcols(ROIs)))==1)
            names=as.data.frame(mcols(ROIs[m[!is.na(m)]]))[,1]
        else
            names=apply(X=as.data.frame(mcols(ROIs[m[!is.na(m)]])),
                FUN=paste,MARGIN=1, collapse="_")
    }else{
        names=paste0(seqnames(getRegions(qs)[keep]), ":", 
            start(getRegions(qs)[keep]), "-", end(getRegions(qs)[keep]))
    }
    if(length(keep)<=10){
        stop("not enough regions left:",length(keep)," Regions, minimum is 10")
    }
    vals=getNormalizedValues(qs,methods=norm_method, 
        windows=keep, samples=samples)
    missing=apply(X=is.na(vals), MARGIN=1, FUN=any)
    vals=vals[!missing,]-rowMeans(vals[!missing,])
    if(any(!is.finite(vals))){
        stop("Infinite values due to log or logit transformation. ",
            "Please specify minVal and maxVal.")
    }
    if(!is.null(topVar) && nrow(vals)>topVar){
        rv=apply(vals,1,var)
        th=sort(rv,decreasing=TRUE)[topVar]
        vals=vals[rv>=th,]
        names=names[rv>=th]
    }
    svdVals=svd(vals)
    new('qseaPCA', svd=svdVals, sample_names=samples, 
        factor_names=as.character(names))
}
        
plotPCA.qsea=function(object,plotComponents=c(1,2), fgColor="black", 
        bgColor="white", legend=NULL, plotLabels=TRUE, radius=5, 
        labelOffset=.5,labelPos=1, labelAdj=NULL, labelColor="black", 
        cex=1, ...) {
    pca=object
    pc1=plotComponents[1]
    pc2=plotComponents[2]
    svd=getSVD(pca)
    n=length(svd$d)
    if(length(radius)==1) radius=rep(radius, n)
    x=svd$v[,pc1]*svd$d[pc1]
    y=svd$v[,pc2]*svd$d[pc2]
    symbols(x=x, 
        y=y, fg=fgColor,bg=bgColor, pch=20, 
        xlab=paste0("PC",pc1), ylab=paste0("PC",pc2), circles=radius,
        inches=min(radius)/100 , ... )
    if(plotLabels)
        text(svd$v[,pc1]*svd$d[pc1], svd$v[,pc2]*svd$d[pc2], 
            label=getSampleNames(pca), pos=labelPos,adj=labelAdj, 
            offset=labelOffset,col=labelColor, cex=cex)
    if(!is.null(legend))
        legend(x=legend, legend=getSampleNames(pca), fill=bgColor, cex=cex)
    invisible(list(x=x, y=y))

}
setMethod('plotPCA', signature(object='qseaPCA'), plotPCA.qsea)

setGeneric('plotPCAfactors', function(object,...) 
    standardGeneric('plotPCAfactors'))

plotPCAfactors.qsea=function(object,plotComponents=c(1,2), fgColor="black", 
        bgColor="white", plotTopLabels=100,labelsOfInterest,radius=1, 
        labelOffset=.5,labelPos=1,labelColor="black", cex=1 , ... ){
    pca=object
    pc1=plotComponents[1]
    pc2=plotComponents[2]    
    svd=getSVD(pca)
    n=nrow(svd$u)
    plotLabels=numeric(0)
    if(!missing(labelsOfInterest))
        plotLabels=match(labelsOfInterest, getFactorNames(pca))
    x=svd$u[,pc1]*svd$d[pc1]
    y=svd$u[,pc2]*svd$d[pc2]
    if(plotTopLabels>0){
        if(plotTopLabels>n)
            plotLabels=seq_len(n)
        else{
            dist=sqrt(x^2+y^2)
            th=sort(dist, decreasing=TRUE)[plotTopLabels]
            plotLabels=unique(sort(c(plotLabels, which(dist>=th) ) ) )
        }
    }
    if(length(radius)==1) radius=rep(radius, n)
    symbols(x=x, y=y, fg=fgColor,bg=bgColor, pch=20, xlab=paste0("PC",pc1), 
        ylab=paste0("PC",pc2), circles=radius, inches=min(radius)/100, ...)
    if(length(plotLabels)>0)
        text(x=x[plotLabels], y=y[plotLabels], 
            label=getFactorNames(pca)[plotLabels],pos=labelPos,
            offset=labelOffset, col=labelColor, cex=cex)
    invisible(list(x=x, y=y))

}
setMethod('plotPCAfactors', signature(object='qseaPCA'), plotPCAfactors.qsea)

#########
#CNV    #
#########



plotCNV<-function(qs, dist=c("euclid", "cor")[1], clust_method="complete", 
    chr= getChrNames(qs), samples=getSampleNames(qs),cex=1, 
    labels=c(TRUE,TRUE, TRUE, TRUE), naColor="darkgrey" , indicateLogFC=TRUE){

    if(missing(qs) | class(qs) != "qseaSet" )
        stop("No qseaSet specified!")
    if(length(getCNV(qs))==0)
        stop("qs does not contain CNV information. Run \"addCNV\" first")
    n=length(samples)
    if(is.numeric(samples))
        samples=getSampleNames(qs)[samples]

    colF=colorRamp(c("lightgreen","green", "darkgreen", "darkblue", "darkred", 
        "red", "orange"))
    tickpos=0
    chrorder=numeric(0)
    chrpos=numeric(0)

    for (c in chr){
        idx=which(seqnames(getCNV(qs))==c)
        l=length(tickpos)
        cl=length(idx)
        chrpos=c(chrpos, tickpos[l]+cl/2)
        tickpos=c(tickpos, tickpos[l]+cl)
        chrorder=c(chrorder,idx)
    }
    chrpos=chrpos/max(tickpos)
    tickpos=tickpos/max(tickpos)
    CNVval=as.matrix(mcols(getCNV(qs)))[chrorder,samples]
    if(dist=="cor"){
        dist_ma=as.dist(1-cor(CNVval, use="pair"))
        dist_ma[is.na(dist_ma)]=0
    }else if(dist=="euclid"){
        dist_ma=dist(t(CNVval))
    }else stop("distance must be \"euclid\" or \"cor\"")
    #tree
    clust=hclust(as.dist(dist_ma), method=clust_method)
    space=max(nchar(as.character(samples )))*cex*.6-1
    lim=quantile(abs(CNVval),.99, na.rm=TRUE)
    lim=max(abs(CNVval), na.rm=TRUE)
    CNVval[is.na(CNVval)]=lim*42/41

    #par( mai=c(3.1,2,2.1,1)/6 )

    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
        widths=c(1,4), heights=c(1,8))
    #plot(1,axes=FALSE,xlim=c(0,.1),xlab="", ylab="")
    sav=par(mar=c(2,1,2,1))
    #legend
    image(matrix(0:40),col=rgb(colF(0:40/40),maxColorValue=255), axes=FALSE)
    axis(side=1,labels=c("deletion","neutral", "amplification" ),
        at=c(0.1,0.5,.9), cex.axis=cex)
    if(indicateLogFC){
        markAt=1
        if(lim<1) markAt=.5
    
        text((markAt*c(-1,0,1)/lim+1)/2,c(.5,.5,.5),markAt*c(-1,0,1) , cex=cex)
    }
    #values    
    par( mar=c(3+2*cex*labels[1],2,2+2*cex*labels[3],1) )
    
    plot(as.dendrogram(clust),horiz=TRUE,axes=FALSE,leaflab="none",yaxs="i")
    #image
    par( mar=c(3+2*cex*labels[1],
        1+space*labels[2],
        2+2*cex*labels[3],
        2+space*labels[4]))# b, l, t, r
    image(x=seq(0,1,length.out=nrow(CNVval)+1), z=CNVval[,clust$order],
        zlim=c(-lim,lim*42/41), axes=FALSE,, xlim=c(0,1), xlab="",
        col=c(rgb(colF(0:40/40),maxColorValue=255),naColor))
    abline(v=tickpos)
    if(labels[1])
        axis(side=1,labels=chr,at=chrpos,cex.axis=1,las=2, tick=FALSE, 
            cex.axis=cex)
    axis(side=1,labels=rep("", length(tickpos)),at=tickpos,cex.axis=cex,
        las=2, tick=TRUE)
    if(labels[2])
        axis(side=2,labels=samples[clust$order],at=(seq_len(n)-1)/(n-1),
            cex.axis=cex,las=2)
    if(labels[3])
        axis(side=3,labels=chr,at=chrpos,cex.axis=1,las=2, 
            tick=FALSE, cex.axis=cex)
    axis(side=3,labels=rep("", length(tickpos)),at=tickpos,cex.axis=cex,las=2, 
            tick=TRUE)
    if(labels[4])
        axis(side=4,labels=samples[clust$order],at=(seq_len(n)-1)/(n-1),
            cex.axis=cex,las=2)
    par(sav)
    par(mfrow=c(1,1))
    invisible(dist_ma)
}
###########
#Coverage #
###########


plotCoverage<-function(qs,test_results, chr, start, end, samples,samples2, 
        norm_method="nrpkm", yoffset, xlab="Position", ylab="MeDIP seq", 
        col="black", main, reorder="non", indicate_reorder=TRUE, distfun=dist, 
        clustmethod="complete", scale=TRUE, steps=TRUE, space=0.05, 
        baselines=TRUE, scale_val, scale_unit=NULL, logFC_pc=.1, cex=1, 
        smooth_width, smooth_function=mean, regions, regions_lwd=1, 
        regions_col="black", regions_offset, regions_names, regions_dash=0){
    if(! missing(regions)){
        if(class(regions)!="list"){
            regions=list(ROIs=regions)
        }
        if(class(norm_method)=="character"){
            norm_method=normMethod(norm_method)
        }
        if(! all(sapply(regions, class)=="GRanges")){
            stop("\"regions\" should be a \"GenomicRanges\" object")
        }
        if(! missing(regions_lwd) && class(regions_lwd) != "list"){
            regions_lwd=rep(list(regions_lwd), length(regions))
            names(regions_lwd)=names(regions)
        }
        if(! missing(regions_col) && class(regions_col) != "list"){
            regions_col=rep(list(regions_col), length(regions))
            names(regions_col)=names(regions)
        }
        if(! missing(regions_dash) && class(regions_dash) != "list"){
            regions_dash=rep(list(regions_dash), length(regions))
            names(regions_dash)=names(regions)
        }
        if( missing(regions_offset) || class(regions_offset) != "list"){
            regions_offset=rep(list(regions_offset), length(regions))
            names(regions_offset)=names(regions)
        }
        if( missing(regions_names) || class(regions_names) != "list"){
            regions_names=rep(list(regions_names), length(regions))
            names(regions_names)=names(regions)
        }
    }
    if(length(col)<length(samples)) 
        col=rep(col,ceiling(length(samples)/length(col) ))
    if(missing(main)) main=paste0(chr,":",start, "-", end)
    if(!missing(samples2) && length(samples) != length(samples2) ) 
        stop("length(samples) != length(samples2)")
    if(missing(test_results) && reorder=="minP")
        stop("please specify \"test_results\" for reorder=\"minP\"")
    if(is.numeric(reorder) && (reorder< start || reorder > end))
        stop("please specify reordering position within ",
            "the specified plotting region")
    message("selecting specified Region")

    roi_tab=makeTable(qs=qs, samples=samples, 
        ROIs=GRanges(chr, ranges=IRanges(start, end)),norm_methods=norm_method)
    ret=list(regions=roi_tab[,1:4])    
    roi_val=roi_tab[,grep(paste0("_",norm_method), names(roi_tab)), drop=FALSE]
    if(!missing(samples2)){
        roi_tab2=makeTable(qs=qs, samples=samples2, 
            ROIs=GRanges(chr, ranges=IRanges(start, end)),
            norm_methods=norm_method)
        roi_val2=roi_tab2[,
            grep(paste0("_",names(norm_method)), names(roi_tab2)), drop=FALSE]
        roi_val=log2((roi_val+logFC_pc)/(roi_val2+logFC_pc))
        colnames(roi_val)=
            paste0("log2(", colnames(roi_val),"/",colnames(roi_val2) ,")")
    }
    roi_val[is.na(roi_val)]=0
    ret$values=roi_val
    if(! missing(smooth_width)){
        message("smoothing values")
        for(i in seq_len(ncol(roi_val))){
            roi_val[,i]=zoo::rollapply(data=roi_val[,i], 
                width=smooth_width, fill=0,FUN=smooth_function )
        }        
    }
    ord=seq_along(samples)
    ordPos=NA
    if(reorder == "clust"){
        clust=hclust(distfun(t(roi_val)),method=clustmethod)
        ord=clust$order
        ret[["clustering"]]=clust
    }else if(reorder=="max"){
        maxrow=which.max(rowSums(roi_val))
        ord=order(roi_val[maxrow,])
        ret$ordPos=(roi_tab$window_start[maxrow]+
            roi_tab$window_end[maxrow])/2
    }else if(reorder=="minP"){
        test_tab=makeTable(qs, test_results, 
            ROIs=GRanges(chr, ranges=IRanges(start, end)))
        ret$regions=test_tab[,-(5:7)]
        minProw=which.min(test_tab[,grep("PValue", colnames(test_tab))])
        ord=order(roi_val[minProw,])
        ret$ordPos=(roi_tab$window_start[minProw]+
            roi_tab$window_end[minProw])/2
    }else if (is.numeric(reorder)){
        ordRow=sum(roi_tab$window_start<=reorder)
        ord=order(roi_val[ordRow,])
        ret$ordPos=(roi_tab$window_start[ordRow]+
            roi_tab$window_end[ordRow])/2
    }else if(reorder != "non"){
        warning("unknown reorder method")
    }        
    roi_val=roi_val[,ord, drop=FALSE]
    col=col[ord]

    n=ncol(roi_val)
    m=nrow(roi_val)
    maxval=max(roi_val)
    if(missing(yoffset)) yoffset=maxval/4
    roi_val=t(t(roi_val)+(seq_len(n)-1)*yoffset)
    xval=roi_tab[,2] #the window start
    if(steps){
        roi_val=roi_val[as.vector(t(matrix(seq_len(m),m,2,byrow=FALSE))),,
            drop=FALSE]
            #double each line
        xval=as.vector(t(as.matrix(roi_tab[,2:3, drop=FALSE])))
            #start and end
    }
    xlim=range(xval, na.rm=TRUE, finite=TRUE)
    space_val=(xlim[2]-xlim[1])*space
    xlim[1]=xlim[1]-space_val
    if(scale){
        xlim[2]=xlim[2]+space_val
    }
    ylim=range(roi_val, na.rm=TRUE, finite=TRUE)
    ret$xlim=xlim
    ret$ylim=ylim
    reg_bins=list()
    if(!missing(regions) && length(regions)>0){
        for(reg_n in names(regions)){
            regions[[reg_n]]=
                regions[[reg_n]][
                    overlapsAny(regions[[reg_n]],
                    GRanges(chr, IRanges(start, end))) ]
            if(!missing(regions_names[[reg_n]]) && 
                length(regions_names[[reg_n]])==1 && 
                regions_names[[reg_n]] %in% colnames(mcols(regions[[reg_n]])))
                    regions_names[[reg_n]]=
                        mcols(regions[[reg_n]])[,regions_names[[reg_n]]]
            if(!missing(regions_names[[reg_n]]) && 
                length(regions_names[[reg_n]])<length(regions[[reg_n]]) )
                    regions_names[[reg_n]]=
                        rep(regions_names[[reg_n]],
                            ceiling(length(regions[[reg_n]])/
                                length(regions_names[[reg_n]]) ))
            if(length(regions_col[[reg_n]])<length(regions[[reg_n]])) 
                regions_col[[reg_n]]=
                    rep(regions_col[[reg_n]],
                        ceiling(length(regions[[reg_n]])/
                            length(regions_col[[reg_n]]) ))
            reg_bins[[reg_n]]=
                disjointBins(regions[[reg_n]], ignore.strand =TRUE)
            if(missing(regions_offset[[reg_n]])) 
                regions_offset[[reg_n]]=yoffset*1.5
            if(length(regions[[reg_n]])>0)            
                ylim[1]=ylim[1]-(max(reg_bins[[reg_n]])+1)*
                    regions_offset[[reg_n]]
            else
                ylim[1]=ylim[1]-regions_offset[[reg_n]]
        }
    }
    plot(-1,-1,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,axes=FALSE,main=main)
    if (scale){
        if(missing(scale_val))
            scale_val=round(maxval,ceiling(-log10(maxval)))
        scale_x=rep(xlim[2],6)-c(.8,.6,.7,.7,.8,.6 )*space_val
        scale_y=c(0,0,0,1,1,1)*scale_val
        points(x=scale_x, y=scale_y,type="l")
        text(scale_x[2], mean(scale_y), paste(scale_val, scale_unit), pos=4)
    }
    box()
    axis(1)
    for(i in seq_len(n))
    {
        if(baselines)
            points(range(xval, na.rm=TRUE, finite=TRUE),rep(yoffset*(i-1), 2),
                type="l", col="lightgrey", lty=4)
        points(xval, roi_val[,i, drop=FALSE], type="l", col=col[i])
    }
    text(rep(xlim[1],n) + (xlim[2]-xlim[1])*(space-.01),(seq_len(n)-1)*yoffset, 
        colnames(roi_val), col=col,adj = c(1, NA), cex=cex)
    if(indicate_reorder && ! is.null(ret$ordPos)) 
        arrows(ret$ordPos, 0, ret$ordPos, ylim[2]*.9, col="blue")
    if(!missing (regions)&& length(regions)>0){
        offsum=0
        for(reg_n in names(regions)){
            if(length(regions[[reg_n]])>0)    {
                for(i in seq_along(regions[[reg_n]])){
                    ypos=offsum-regions_offset[[reg_n]]*reg_bins[[reg_n]][i]
                    points(c(start(regions[[reg_n]])[i],
                        end(regions[[reg_n]])[i]), rep(ypos,2), 
                        type="l", col=regions_col[[reg_n]][i], 
                        lwd=regions_lwd[[reg_n]])
                    if(! is.null(regions_dash[[reg_n]]) && 
                        regions_dash[[reg_n]]>0){
                        points(c(start(regions[[reg_n]])[i],
                            start(regions[[reg_n]])[i]), 
                            ypos+regions_offset[[reg_n]]*
                                regions_dash[[reg_n]]*c(-1,1), 
                            type="l", col=regions_col[[reg_n]][i], 
                            lwd=regions_lwd[[reg_n]])
                        points(c( end(regions[[reg_n]])[i],
                            end(regions[[reg_n]])[i]), 
                            ypos+regions_offset[[reg_n]]*
                                regions_dash[[reg_n]]*c(-1,1), 
                            type="l", col=regions_col[[reg_n]][i], 
                            lwd=regions_lwd[[reg_n]])
                    }
    
                }
                if(!missing(regions_names[[reg_n]]) )
                    text((start(regions[[reg_n]])+end(regions[[reg_n]]) )/2, 
                        offsum-regions_offset[[reg_n]]*reg_bins[[reg_n]], 
                        regions_names[[reg_n]] , adj=c(.5,1.2),
                        col=regions_col[[reg_n]])
                off_roi=(max(reg_bins[[reg_n]])+1)*regions_offset[[reg_n]]
            }else{
                off_roi=regions_offset[[reg_n]]
            }
            message(reg_n)
            text(xlim[1]+ (xlim[2]-xlim[1])*(space-.01),offsum-off_roi/2, 
                reg_n,col=regions_col[[reg_n]], adj=c(1,.5))
            offsum=offsum-off_roi
        }
    }
    invisible(ret)
}

plotEnrichmentProfile<-function(qs,sample, sPoints=seq(0,30,1), 
    fitPar=list(lty=2,col="green"),cfPar=list(lty=1),densityPar,meanPar,...){
    if(!hasEnrichment(qs))
        stop("no enrichment analysis found")
    if(missing(sample) || length(sample)!=1 )
        stop("please select one sample")
    plotPar=list(...)
    defVal=list(xlab=paste(getEnrichmentPattern(qs),"density") , 
        ylab="expected rpkm", main=sample)
    plotPar=c(plotPar, defVal[! names(defVal) %in% names(plotPar) ])
    if(is.null(plotPar$xlim))
        plotPar$xlim=range(sPoints)
    maxY=1
    recycle=c( "lty", "lwd", "col")
    cfPar=c(cfPar, 
    plotPar[(! names(plotPar) %in% names(cfPar)) & 
        names(plotPar) %in% recycle ])
    factors=getEnrichmentFactors(qs, sample)+getOffset(qs,sample)
    if(is.null(plotPar$ylim))
        maxY=c(maxY,max(factors,na.rm=TRUE))
    if(! missing(densityPar)| ! missing (meanPar)){
        vals=getNormalizedValues(qs, methods=normMethod("nrpkm"), 
            windows=NULL, samples=sample)[,1]
        cfactors=
            mcols(getRegions(qs))[,paste0(getEnrichmentPattern(qs),"_density")]
    }    
    if(!missing(meanPar)){
        breakPoints=(sPoints[-1]+sPoints[-length(sPoints)])/2
        group=findInterval(x=cfactors, vec=breakPoints)
        lval=split(vals, group)
        means=sapply(lval, mean, na.rm=TRUE)
        maxY=max(c(maxY,means))
    }
    par=getEnrichmentParameters(qs, sample)
    if(!is.null(par)){        
        fitted=.sigmF2(sPoints,par[1], par[2], par[3] )+getOffset(qs,sample)
        maxY=max(c(maxY,fitted))
    }
    if(is.null(plotPar$ylim))
        plotPar$ylim=c(0,maxY)
    if(! missing(densityPar)){
        densityPar[["y"]]=vals
        densityPar[["x"]]=cfactors
        densityPar=
            c(densityPar,plotPar[!names(plotPar) %in% names(densityPar)])
        do.call(smoothScatter, densityPar)
    }else(
        do.call(plot, c(x=NA, plotPar))
    )
    retList=list(x=sPoints)
    if(!missing(meanPar)){
        f=is.finite(means)
        do.call("lines", c(list(x=sPoints[f], y=means[f]),meanPar))
        retList$window_means=means
    }
    
    if(!is.null(par)){        
        f=is.finite(sPoints)
        do.call("lines",c(list(x=sPoints[f], y=fitted[f]), fitPar))
        retList$fitted=fitted
    }
    dens=getEnrichmentDensity(qs)
    f=is.finite(factors)& is.finite(dens)
    do.call("lines", c(list(x=dens[f], y=factors[f]),cfPar))
    invisible(retList)

}

plotEPmatrix<-function(qs, sa=getSampleNames(qs),
    nrow=ceiling(sqrt(length(sa))), ncol=ceiling(length(sa)/nrow), ...){
    param=list(...)
    ma=matrix(c(rep(1,ncol),(seq_len((nrow*ncol)))+3,rep(3,ncol)),
        nrow+2,ncol,byrow=TRUE)
    ma=cbind(rep(0,nrow+2),ma, rep(0,nrow+2))
    ma[seq_len(nrow)+1,1]=2
    wd=c(.05,rep(.9/ncol,ncol),.05)
    ht=c(.05,rep(.9/nrow,nrow),.05)
    par( mar=c(1,1,1,1))
    lo=layout(ma, widths=wd, heights=ht)
    #layout.show(lo)
    frame()
    text(0.5,0.5,labels=paste("CpG Enrichment Analysis"),cex=3)
    frame()
    text(0.5,0.5,labels=paste("expected rpkm at full methylation"),
        cex=2.5, srt=90)
    frame()
    text(0.5,0.5,labels=paste("CpG density [#/fragment]"),cex=2.5, srt=0)
    i=0
    if(is.null(param$ylim))
        param$ylim=c(0,
            max(getEnrichmentFactors(qs,sa, minN=10), na.rm=TRUE) )
    retList=list()
    for(s in sa){
        i=i+1
        message(s)
        retList[[s]]=do.call(plotEnrichmentProfile,
            c(list(qs, s, axes=FALSE, cex.main=2, 
            cfPar=list(col="black", lty=1, lwd=1.5),
            fitPar=list(col="darkgreen")),param) )    
        box()
        if(i%%ncol==1) axis(2, cex.axis=1.5)
        if(i>length(sa)-ncol) axis(1, cex.axis=1.5)
        
    }
    par(mfrow=c(1,1))
    invisible(retList)

}


addLine=function(x,y,bins=NULL,fun=median ,...){
    if(is.null(bins)){
        bins=seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=20)
    }
    n=length(bins)-1
    med=numeric(n)
    for(i in seq_len(n)){
        med[i]=fun(y[x>=bins[i]&x<bins[i+1]], na.rm=TRUE)
    }
    lines((bins[-n-1]+bins[-1])/2,med,...)
    invisible(med)
}


