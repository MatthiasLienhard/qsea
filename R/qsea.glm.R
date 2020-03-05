
fitNBglm=function(qs,design,link="log",keep, disp_method="region_wise",
    norm_method="rpkm",init_disp=0.5 ,verbose=TRUE,minRowSum=10,pseudocount=1, 
    disp_iter=3,nChunks=NULL, parallel=FALSE ){

    if(! disp_method %in% 
        c("initial","common", "region_wise", "cutAtQuantiles"))
        stop("disp_method should be in \"initial\",\"common\",",
            "\"cutAtQuantiles\" or \"region_wise\"")
    sampleIdx=row.names(design)
    sampleIdx=checkSamples(qs,sampleIdx)
    fitDisp=disp_method!="initial"
    if(is(norm_method,"character"))
        norm_method=normMethod(norm_method)
    if(is(norm_method,"normMethods") & length(norm_method)==1)
        norm_method=norm_method[[1]]
    else stop("please provide one valid norm_method")
    if( ! is(design,"matrix")) #todo: some more checks: nrow, full rank...!!!
        stop("please provide valid design matrix")
    if(missing(qs) || ! is(qs, "qseaSet"))
        stop("please provide a qseaSet")
    if(link !="log")
        stop("currently, \"link\" must be \"log\"")
    message("selecting regions with at least ",minRowSum," reads...")
    
    if(missing(keep)) keep=seq_along(getRegions(qs))
    if(length(init_disp) !=1 && length    (init_disp) !=length(keep) )
        stop("for the initial dispersion, ",
            "provide either one value for all regions, ",
            "or one value for each region")
    if(parallel) 
        BPPARAM=bpparam()
    else
        BPPARAM=SerialParam()
    sel=rowSums(getCounts(qs, windows=keep, samples=sampleIdx)) >= minRowSum
    keep=keep[sel]
    if(length(init_disp)>1)
        init_disp=init_disp[sel]
    if(length(keep)==0)
        stop("no regions left")

    #get normalization factors
    message("deriving normalization factors...")
    nm=getNormMatrix(qs, methods=norm_method,windows=keep,sampleIdx)
    if("enrichment" %in% norm_method){
       pattern_name=getEnrichmentPattern(qs)
       cf=mcols(getRegions(qs))[keep,paste0(pattern_name, "_density")]
       sPar=getEnrichmentParameters(qs)
       nm$factors=nm$factors*mapply(.sigmF2, MoreArgs=list(x=cf),
            a=sPar[sampleIdx,"a"],
            b=sPar[sampleIdx,"b"],
            c=sPar[sampleIdx,"c"], 
            o=nm$offset,
            USE.NAMES = FALSE)+pseudocount
    }
    normF=nm$factors
    rm(nm)
    # filter invalid regions
    if(any(!is.na(normF) | normF <=0)){
        normF[normF<=0]=-Inf
        #isfinite=apply(is.finite(normF), 1,all)
        isfinite=is.finite(rowSums(normF))
        normF=normF[isfinite,]
        keep=keep[isfinite]
        if(length(init_disp)>1)
            init_disp=init_disp[isfinite]
    }
    if(length(keep)==0)
        stop("no regions left")

    #get raw counts
    if(verbose) message("obtaining raw values for ",length(sampleIdx),
        " samples in ",length(keep)," windows")        
    y=getCounts(qs, windows=keep, samples=sampleIdx)+pseudocount    
    if(pseudocount==0)
        y[y==0]=1

    coef=NULL
    if(fitDisp)    {
        message("initial fit to estimate dispersion")
        fullModel=fitNBglmMatrix(design, y=y, nf=normF, linkf=link, 
            disp=init_disp,fitDisp=TRUE, maxit=disp_iter,nChunks=nChunks,
            BPPARAM=BPPARAM)
        if(disp_method=="region_wise"){
            disp=fullModel$dispersion
        }else if(disp_method=="common"){
            disp=mean(fullModel$dispersion)
        }else if(disp_method=="cutAtQuantiles"){
            disp=fullModel$dispersion
            dispQ=quantile(disp, c(.25,.75))
            disp=pmax(disp, dispQ[1])
            disp=pmin(disp, dispQ[2])
        }else{
            warning("unknown dispersion method")
            disp=fullModel$dispersion
        }
    }else{
        fullModel=list(coefficients=NULL)
        disp=init_disp
    }
    message("fitting GLMs")    
    fullModel=fitNBglmMatrix(design, y=y, nf=normF, coef=fullModel$coefficients,
            linkf=link, disp=disp,fitDisp=FALSE, nChunks=nChunks,
            BPPARAM=BPPARAM)
    fullModel$dispersion=disp
    return(new('qseaGLM', 
        fullModelDesign=design, 
        contrast=list(), 
        parameters=list(link=link, norm_method=norm_method,
            pseudocount=pseudocount),
        fullModel=fullModel,
        windows=keep
    ))
}


addContrast=function(qs,glm,contrast,coef,name,verbose=TRUE,
    nChunks=NULL, parallel=FALSE){
    sampleIdx=getSampleNames(glm)
    if(missing(name))
        stop("please specify a name for the contrast")
    if( missing(contrast) == missing(coef) )
        stop("either contrast or coeficient must be set")
    if(!missing(contrast)){ 
        if(is(contrast,"vector"))
            contrast=as.matrix(contrast)
        if(dim(contrast)[2]>1)
            stop("currently, only one contrast is supported")
        design=limma::contrastAsCoef(getDesignMatrix(glm), contrast)$design 
        coef=1
        effect=getFullModel(glm)$coefficients %*% contrast
    }else{
        design=getDesignMatrix(glm)
        effect=getFullModel(glm)$coefficients[, coef]
    }
    if(is(coef,"character"))
        coef=match(coef, getCoefNames(glm))
    if(any(is.na(coef)))
        stop("contrast not found in design")
    design=design[,-coef, drop=FALSE]
    param=getParameters(glm)
    if(parallel) BPPARAM=bpparam() 
    else   BPPARAM=SerialParam()
    
    nm=getNormMatrix(qs, methods=param$norm_method,
        windows=getWindows(glm),sampleIdx)
    if("enrichment" %in% param$norm_method){
       pattern_name=getEnrichmentPattern(qs)
       cf=mcols(getRegions(qs))[getWindows(glm),paste0(pattern_name, "_density")]
       sPar=getEnrichmentParameters(qs)
       nm$factors=nm$factors*mapply(.sigmF2, MoreArgs=list(x=cf),
            a=sPar[sampleIdx,"a"],
            b=sPar[sampleIdx,"b"],
            c=sPar[sampleIdx,"c"], 
            o=nm$offset,
            USE.NAMES = FALSE)+param$pseudocount
    }
    normF=nm$factors
    rm(nm)

    y=getCounts(qs, windows=getWindows(glm), samples=sampleIdx)+
        param$pseudocount
    if(param$pseudocount==0)
        y[y==0]=1 
    #y==0 for all samples is problematic --> quick n dirty fix ^
    #pmax(y,1)does not work for large matrices

    fit=fitNBglmMatrix(design, y=y, nf=normF, linkf=param$link,
        disp=getFullModel(glm)$dispersion, fitDisp=FALSE,
        nChunks=nChunks,BPPARAM=BPPARAM)

    df=fit$rdf-getFullModel(glm)$rdf
    LR=fit$deviance-getFullModel(glm)$deviance
    fit$LRT_pval=pchisq(LR,df, lower.tail = FALSE)
    fit$LRT_adjP=p.adjust(fit$LRT_pval, method="BH")
    fit$effect=effect
    setReducedModel(glm, name=name, fit)

}


fitNBglmMatrix=function(x,y,disp, nf,linkf="log", maxit=60, coef=NULL,
    eps=10e-6,  fitDisp=FALSE,nChunks=NULL,BPPARAM){
    #this is an attempt to increase performance of glm.fit for this special case
    #by vector operations
    
    #link function
    #nf=norm$factors
    #offset=offset+pseudocount
    
    #y[y<offset]=offset[y<offset]
    #offset[offset>0]=0
    fam=list(linktype=linkf, eps=eps)
    if(linkf=="log"){
        fam$link=function(mu,nf) log(mu/nf)
        fam$linkInv=function(eta,nf) pmax(exp(eta)*nf, .Machine$double.eps)
        fam$dLinkInv=function(eta,nf) pmax(exp(eta)*nf, .Machine$double.eps)

    #link=linear is not functional at the moment, maybe change to softmax?
    }else if (linkf=="linear"){
        fam$link=function(mu,nf) mu/nf
        #linkInv=function(x,nf, o) pmax(x*nf+t(o+t(nf)), eps)
        fam$linkInv=function(eta,nf) pmax(eta*nf,eps)
        fam$dLinkInv=function(eta,nf) nf
    }else(stop ("unknown link function"))
    fam$variance=function(mu, disp) mu+mu^2*disp
    fam$dev.resids=function(y, mu, disp, eps=.Machine$double.eps^0.25) {
        theta=1/pmax(disp,eps)
        2*(y*log(pmax(1, y)/mu)-(y+theta)*log((y+theta)/(mu+theta)))
    }
    n=nrow(y)       
    nvars <- ncol(x) 
    if(n<1e5) 
        BPPARAM=SerialParam()
    if(is.null(nChunks))
        nChunks=BiocParallel::bpnworkers(BPPARAM)
    #todo=seq_len(n)
    if(nChunks>1){
        message("Fitting GLMs in ",nChunks , " chunks")
        binF= cut(seq_len(n), nChunks, labels = FALSE)
    }else {
        binF=1
    }
    todoL=split(seq_len(n), binF)
    if(!is.null(coef))
        coef=lapply(split(coef, binF), matrix, ncol=nvars)
    else
        coef=rep(list(NULL),nChunks)

    if(length(disp)>1)
        disp=split(disp, binF)       
    rdf=nrow(x)-qr(x)$rank

    fitted=unlist(bpmapply(todo=todoL,coef=coef, disp=disp, #those change
            FUN=.fitNBglmMatrixPart, 
            MoreArgs=list(#those do not change
                y=y,nf=nf,x=x,fam=fam,rdf=rdf, maxit=maxit, fitDisp=fitDisp), 
            BPPARAM=BPPARAM , SIMPLIFY=FALSE),FALSE, FALSE)
    coef=t(matrix(unlist(fitted[seq(1,length(fitted), 4)], FALSE, FALSE),
        nrow=nvars))
    if(fitDisp )
        disp=unlist(fitted[seq(2,length(fitted), 4)], FALSE, FALSE)
    else
        disp=NULL
    conv=unlist(fitted[seq(3,length(fitted), 4)], FALSE, FALSE)
    dev=unlist(fitted[seq(4,length(fitted), 4)], FALSE, FALSE)
            
    return(list(coefficients=coef,dispersion=disp, 
        converged=conv, rdf=rdf, deviance=dev))

}


.fitNBglmMatrixPart<-function(x,y,nf,todo, coef, disp,rdf, fam, maxit,fitDisp){ 
    n=length(todo)
    idxO=todo[1]-1
    todo=todo-idxO
    nvars <- ncol(x)
    disp_old=disp=rep(disp, length.out=n)    
    dev_old = dev = rep(Inf, n)
    if(is.null(coef)){
        coef_old=NULL
        mu=y[todo+idxO,,drop=FALSE]
        mu[mu<1/6]=1/6    
        eta=fam$link(mu, nf[todo+idxO,,drop=FALSE])    
        coef=matrix(0,n,nvars)
    }else{
        coef_old=coef
        eta = t(x %*% t(coef))
        mu = fam$linkInv(eta, nf[todo+idxO,,drop=FALSE])
        dev=dev_old=rowSums(fam$dev.resids(y[todo+idxO,,drop=FALSE], 
            mu, disp))
    }
    conv=rep(NA,n)
    newConv=numeric(0)

    #IRLS
    #message("fitting GLMs")
    for (iter in 1L:maxit) {
        message("iteration ", iter,": " )
        mu.eta.val=fam$dLinkInv(eta, nf[todo+idxO,,drop=FALSE])
        #z = eta + (y[todo,] - mu)/mu.eta.val
        w = sqrt((mu.eta.val^2)/fam$variance(mu, disp[todo]))
        zw =w* (eta + (y[todo+idxO,,drop=FALSE] - mu)/mu.eta.val)
        FUN<-function(i){#Cdqrls is taken from stats package
                ret=numeric(nvars)
                fit = .Call("Cdqrls",x*w[i,], zw[i,], tol=fam$eps, 
                    check = FALSE, PACKAGE = "qsea")
                ret[fit$pivot]=fit$coefficients
                ret
        }
        coef[todo, ]=t(vapply(seq_along(todo) , FUN, numeric(nvars)))
        rm(w)
        rm(zw)
        rm(mu.eta.val)
        eta = t(x %*% t(coef[todo,, drop=FALSE]))
        mu = fam$linkInv(eta ,nf[todo+idxO,,drop=FALSE])
        dev[todo] = rowSums(fam$dev.resids(y[todo+idxO,, drop=FALSE], 
            mu, disp[todo])) 
        #sometimes glms increase in deviance... adapt stepsize
        increase=which(dev[todo]>dev_old[todo])
        step=0
        while(length(increase)>0 & step<8){
            coef[todo[increase],]=(coef[todo[increase],,drop=FALSE]+
                coef_old[todo[increase],,drop=FALSE])/2
            eta[increase,] = t(x %*% t(coef[todo[increase],, drop=FALSE]))
            mu[increase,] = fam$linkInv(eta[increase,,drop=FALSE] ,
                nf[todo[increase]+idxO,,drop=FALSE])
            dev[todo[increase]] = 
                rowSums(fam$dev.resids(y[todo[increase]+idxO,, drop=FALSE], 
                    mu[increase,,drop=FALSE], disp[todo[increase]])) 
            increase=increase[dev[todo[increase]]>dev_old[todo[increase]]]
            step=step+1    #if deviance still increased mark glm as "divergent"
        }
        if(!fitDisp){    #remove converged regions    
            newConv=which(abs(dev[todo]-
                dev_old[todo])/(0.1+abs(dev[todo]))<fam$eps) 
            conv[todo[newConv]]=TRUE
            conv[todo[increase]]=FALSE
            newConv=union(newConv, increase)
            message("...",sum(conv, na.rm=TRUE),"/",length(conv),
                " regions converged")
            if (length(todo)-length(newConv)<=0) {
                break
            }else if(length(newConv)>0){
                    todo=todo[-newConv]
                    eta=eta[-newConv,,drop=FALSE]
                    mu=mu[-newConv,,drop=FALSE]
            }
        }
        dev_old = dev
        coef_old=coef

    }
    if(fitDisp){
        disp=1/theta_md_matrix(y[todo+idxO,, drop=FALSE],mu,rdf)        
    }else{disp=NA}
    if(fam$linktype=="log")
        coef=coef/log(2)
    return(list(coefficients=t(coef),dispersion=disp, 
        converged=conv,  deviance=dev))
}

theta_mm_matrix=function (y, mu, dfr,iter = 10,eps = .Machine$double.eps^0.25) 
{
    n = ncol(y)
    t0 = n/rowSums((y/mu - 1)^2)
    it = 0
    del = 1    
    todo=seq_len(nrow(y))
    while ((it = it + 1) < iter && length(todo)>0) {
        t0[todo]= abs(t0[todo])
        del=(rowSums(((y[todo,,drop=FALSE]-
                mu[todo,,drop=FALSE])^2/(mu[todo,,drop=FALSE] + 
                mu[todo,,drop=FALSE]^2/t0))) - dfr)/
            rowSums((y[todo,,drop=FALSE] - 
                mu[todo,,drop=FALSE])^2/(mu[todo,,drop=FALSE] + t0[todo])^2)
        t0[todo] = t0[todo] - del
        todo=todo[del >= eps]
    }
    t0[t0<0]=-Inf
    t0
}

theta_md_matrix=function (y, mu, dfr,iter = 10, eps = .Machine$double.eps^0.25) 
{
    n = ncol(y)
    t0 = n/rowSums((y/mu - 1)^2)

    #a = 2 * rowSums(y * log(pmax(1, y)/mu)) - dfr
    a = 2 * rowSums(y * log(y/mu)) - dfr
    it = 0
    del = 1
    todo=seq_len(nrow(y))
    while ((it = it + 1) < iter && length(todo)>0) {
        t0[todo] = abs(t0[todo])
        tmp = log((y[todo,,drop=FALSE] + t0[todo])/
            (mu[todo,,drop=FALSE] + t0[todo]))
        top = a[todo] - 2 * rowSums((y[todo,,drop=FALSE] + t0[todo]) * tmp)
        bot = 2*rowSums(((y[todo,,drop=FALSE] - mu[todo,,drop=FALSE])/
            (mu[todo,,drop=FALSE] + t0[todo]) - tmp))
        del = top/bot
        del[is.na(del)]=0
        t0[todo] = t0[todo] - del
        todo=todo[del >= eps]
    }
    t0[t0<0] = -Inf
    t0
}


