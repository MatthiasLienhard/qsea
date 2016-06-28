
fitNBglm=function(qs,design,link="log",keep, disp_method="region_wise",
    norm_method="rpkm",init_disp=0.5 ,verbose=TRUE,minRowSum=10,pseudocount=1){

    if(! disp_method %in% 
        c("initial","common", "region_wise", "cutAtQuantiles"))
        stop("disp_method should be in \"initial\",\"common\",",
            "\"cutAtQuantiles\" or \"region_wise\"")
    sampleIdx=row.names(design)
    fitDisp=disp_method!="initial"
    if(class(norm_method)=="character")
        norm_method=normMethod(norm_method)
    if(class(norm_method)=="normMethods" & length(norm_method)==1)
        norm_method=norm_method[[1]]
    else stop("please provide one valid norm_method")
    if(class(design) !="matrix") #todo: some more checks: nrow, full rank...!!!
        stop("please provide valid design matrix")
    if(missing(qs) || class(qs) != "qseaSet")
        stop("please provide a qseaSet")
    if(link !="log")
        stop("currently \"link\" must be \"log\"")        
    message("selecting regions with at least ",minRowSum," reads...")
    
    if(missing(keep)) keep=seq_along(getRegions(qs))
    if(length(init_disp) !=1 && length    (init_disp) !=length(keep) )
        stop("for the initial dispersion, ",
            "provide either one value for all regions, ",
            "or one value for each region")
    
    sel=rowSums(getCounts(qs, windows=keep, samples=sampleIdx)) >= minRowSum
    keep=keep[sel]
    if(length(init_disp)>1)
        init_disp=init_disp[sel]
    if(length(keep)==0)
        stop("no regions left")

    #get normalization factors
    message("deriving normalization factors...")
    nm=getNormMatrix(qs, methods=norm_method,windows=keep,sampleIdx, 
        pseudocount,cfCut=0)
    normF=nm$factors
    offset=nm$offset
    #if("offset" %in% norm_method){
    #    offset=t(getOffset(qs, sampleIdx)*t(normF))
    #}

    normF[normF<=0]=-Inf
    isfinite=apply(is.finite(normF), 1,all)
    normF=normF[isfinite,]

    offset=offset[isfinite,]
    
    keep=keep[isfinite]
    if(length(init_disp)>1)
        init_disp=init_disp[isfinite]

    if(length(keep)==0)
        stop("no regions left")

    #get raw counts
    if(verbose) message("obtaining raw values for ",length(sampleIdx),
        " samples in ",length(keep)," windows")        
    y=getCounts(qs, windows=keep, samples=sampleIdx)
    coef=NULL
    if(fitDisp)    {
        message("initial fit to estimate dispersion")    
        fullModel=fitNBglmMatrix(design, y=y, nf=normF,offset, link=link, 
            disp=init_disp,fitDisp=TRUE, pseudocount=pseudocount, maxit=2)
        if(disp_method=="region_wise"){
            init_disp=fullModel$dispersion
        }else if(disp_method=="common"){
            init_disp=mean(fullModel$dispersion)
        }else if(disp_method=="cutAtQuantiles"){
            init_disp=fullModel$dispersion
            dispQ=quantile(init_disp, c(.25,.75))
            init_disp=pmax(init_disp, dispQ[1])
            init_disp=pmin(init_disp, dispQ[2])
        }else{    
            warning("unknown dispersion method")
            init_disp=fullModel$dispersion
        }
        coef=fullModel$coefficients
    }
    message("fitting GLMs")    
    fullModel=fitNBglmMatrix(design, y=y, nf=normF,offset, coef=coef, link=link,
            disp=init_disp,fitDisp=FALSE, pseudocount=pseudocount)
    fullModel$dispersion=init_disp
    return(new('qseaGLM', 
        fullModelDesign=design, 
        contrast=list(), 
        parameters=list(link=link, norm_method=norm_method,
            pseudocount=pseudocount),
        fullModel=fullModel,
        windows=keep
    ))
}


addContrast=function(qs,glm,contrast,coef,name,verbose=TRUE ){
    sampleIdx=getSampleNames(glm)
    if(missing(name))
        stop("please specify a name for the contrast")
    if( missing(contrast) == missing(coef) )
        stop("either contrast or coeficient must be set")
    if(!missing(contrast)){ 
        if(class(contrast)=="vector")
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
    if(class(coef)=="character")
        coef=match(coef, getCoefNames(glm))
    if(any(is.na(coef)))
        stop("contrast not found in design")
    design=design[,-coef, drop=FALSE]
    param=getParameters(glm)
    nm=getNormMatrix(qs, methods=param$norm_method,
        windows=getWindows(glm),sampleIdx, param$pseudocount,cfCut=0)

    normF=nm$factors
    offset=nm$offset
    y=getCounts(qs, windows=getWindows(glm), samples=sampleIdx)
    fit=fitNBglmMatrix(design, y=y, nf=normF,offset=offset, 
        link=param$link, disp=getFullModel(glm)$dispersion, fitDisp=FALSE, 
        pseudocount=param$pseudocount)
    df=fit$rdf-getFullModel(glm)$rdf
    LR=fit$deviance-getFullModel(glm)$deviance
    fit$LRT_pval=pchisq(LR,df, lower.tail = FALSE)
    fit$LRT_adjP=p.adjust(fit$LRT_pval, method="BH")
    fit$effect=effect
    setReducedModel(glm, name=name, fit)

}


fitNBglmMatrix=function(x,y,disp, nf,offset,link="log", maxit=60, coef=NULL,
    eps=10e-6, pseudocount=1, fitDisp=FALSE){
    #this is an attempt to increase performance of glm.fit for this special case
    #by vector operations

    #link function
    #nf=norm$factors
    offset=offset+pseudocount
    y=y+pseudocount
    #y[y<offset]=offset[y<offset]
    #offset[offset>0]=0
    if(link=="log"){
        link=function(mu,nf, o) log2(mu/(nf+o))
        #linkInv=function(x,nf, o) pmax(2^(x)*nf+o,eps)
        linkInv=function(eta,nf, o) 2^(eta)*(nf+o)
        dLinkInv=function(eta,nf, o) 2^(eta)*(nf+o)
        
    }else if (link=="linear"){
        link=function(mu,nf, o) (mu-o)/nf
        #linkInv=function(x,nf, o) pmax(x*nf+o, eps)
        linkInv=function(eta,nf, o) pmax(eta*nf+o,0.1*o)
        dLinkInv=function(eta,nf, o) nf
    }
    else(stop ("unknown link function"))
    variance=function(mu, disp) mu+mu^2*disp
    rdf=nrow(x)-qr(x)$rank

    dev.resids=function(y, mu, disp, eps=.Machine$double.eps^0.25) {
        theta=1/pmax(disp,eps)
        2*(y*log(pmax(1, y)/mu)-(y+theta)*log((y+theta)/(mu+theta)))
    }
    n=nrow(y)
    nvars <- ncol(x)

    disp_old=disp=rep(disp, length.out=n)
    
    dev_old = dev = rep(Inf, length.out=n)
    
    if(is.null(coef)){
        coef_old=NULL
        mu=y
        mu[mu<1/6]=1/6    
        eta=link(mu, nf, offset)    
        coef=matrix(0,n,dim(x)[2])
    }else{
        coef_old=coef
        eta = t(x %*% t(coef))
        mu = linkInv(eta, nf, offset)
        dev=dev_old=rowSums(dev.resids(y, mu, disp))
    }
    conv=rep("no",n)
    newConv=numeric(0)
    fail=numeric(0)
    todo=seq_len(n)
    #IRLS
    #message("fitting GLMs")
    for (iter in 1L:maxit) {

        message("iteraiton ", iter,": " )
        mu.eta.val=dLinkInv(eta, nf, offset)
        z = eta + (y[todo,] - mu)/mu.eta.val
        w = sqrt((mu.eta.val^2)/variance(mu, disp[todo]))
        for(i in seq_along(todo) ){
            fit = .Call("Cdqrls", x*w[i,], z[i,] * w[i,], eps, 
                check = FALSE, PACKAGE = "qsea")
            coef[todo[i],fit$pivot] = fit$coefficients
            #Cdqrls is taken from stats package
            #todo: make the loop in the c function to increase performance
        }
        eta = t(x %*% t(coef[todo,, drop=FALSE]))
        mu = linkInv(eta ,nf, offset)
        dev[todo] = rowSums(dev.resids(y[todo,, drop=FALSE], mu, disp[todo])) 
        #sometimes glms increase in deviance... adapt stepsize
        increase=which(dev[todo]>dev_old[todo])
        step=0
        while(length(increase)>0 & step<8){
            coef[todo[increase],]=(coef[todo[increase],,drop=FALSE]+
                coef_old[todo[increase],,drop=FALSE])/2
            eta[increase,] = t(x %*% t(coef[todo[increase],, drop=FALSE]))
            mu[increase,] = linkInv(eta[increase,,drop=FALSE] ,
                nf[increase,,drop=FALSE], offset[increase,,drop=FALSE])
            dev[todo[increase]] = 
                rowSums(dev.resids(y[todo[increase],, drop=FALSE], 
                    mu[increase,,drop=FALSE], disp[todo[increase]])) 
            increase=increase[dev[todo[increase]]>dev_old[todo[increase]]]
            step=step+1    #if deviance still increased mark glm as "divergent"
        }
        if(!fitDisp){    #remove converged regions    
            newConv=which(abs(dev[todo]-
                dev_old[todo])/(0.1+abs(dev[todo]))<eps) 
            conv[todo[newConv]]="yes"
            conv[todo[increase]]="divergent"
            newConv=union(newConv, increase)
            message("...",sum(conv=="yes"),"/",length(conv),
                " regions converged")
            if (length(todo)-length(newConv)<=0) {
                break
            }else if(length(newConv)>0){
                    todo=todo[-newConv]
                    eta=eta[-newConv,,drop=FALSE]
                    mu=mu[-newConv,,drop=FALSE]
                    nf=nf[-newConv,,drop=FALSE]
                    if(!is.null(dim(offset)))
                        offset=offset[-newConv, ,drop=FALSE]
            }
        }
        dev_old = dev
        coef_old=coef

    }
    if(fitDisp){
        disp=1/theta_md_matrix(y,mu,rdf)        
    }else{disp=NA}
    return(list(coefficients=coef,dispersion=disp, 
        converged=conv, rdf=rdf, deviance=dev))
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
                mu[todo,])^2/(mu[todo,,drop=FALSE] + t0[todo])^2)
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

    a = 2 * rowSums(y * log(pmax(1, y)/mu)) - dfr
    it = 0
    del = 1
    todo=seq_len(nrow(y))
    while ((it = it + 1) < iter && length(todo)>0) {
        t0[todo] = abs(t0[todo])
        tmp = log((y[todo,] + t0[todo])/(mu[todo,,drop=FALSE] + t0[todo]))
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


