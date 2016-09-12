

#### qsea object class definition ###
setClass(Class = 'qseaSet', 
    representation = representation(sampleTable='data.frame',
                    count_matrix='matrix',
                    zygosity='matrix',
                    regions='GenomicRanges', 
                    cnv='GenomicRanges',    
                    parameters='list', 
                    libraries='list', 
                    enrichment='list'
                    ),
    prototype = prototype(),
    validity = function(object){
        if(FALSE){stop("")}
        return(TRUE)
    }
)

### show ###
setMethod('show', signature='qseaSet', definition=function(object) {
    message("qseaSet")
    message("=======================================")
    message(nrow(getSampleTable(object) )," Samples: ")
    message(paste(getSampleNames(object),collapse=", ") )
    msg=paste0(length(getRegions(object) ), " Regions in ",
        length(getChrNames(object))," chromosomes")
    if(!is.null(getParameters(object)[["BSgenome"]]))
        msg=paste0(msg, " of ", getParameters(object)[["BSgenome"]])
    message(msg)
    #add some more output
})

#### Methods for extracting slots ###
# sample table
setGeneric('getSampleTable', 
        function(object, ...) standardGeneric('getSampleTable'))
setMethod('getSampleTable', 'qseaSet', 
        function(object, samples=seq_len(nrow(object@sampleTable))) 
            object@sampleTable[samples,] )

# short to get sample names
setGeneric('getSampleNames', 
        function(object,...) standardGeneric('getSampleNames'))
setMethod('getSampleNames', 'qseaSet', 
        function(object, samples=seq_len(nrow(object@sampleTable))) 
                object@sampleTable$sample_name[samples])

# get list with groups
setGeneric('getSampleGroups', function(object,...) 
        standardGeneric('getSampleGroups'))
setMethod('getSampleGroups', 'qseaSet', 
        function(object, samples=seq_len(nrow(object@sampleTable)), 
                group="group") 
                split(object@sampleTable$sample_name[samples], 
                f=object@sampleTable[samples,group]))

# chr names
setGeneric('getChrNames', function(object) 
        standardGeneric('getChrNames'))
setMethod('getChrNames','qseaSet', function(object) 
        as.vector(na.omit(mixedsort(levels(seqnames(object@regions))))) )

# bsgenome
setGeneric('getGenome', function(object) 
        standardGeneric('getGenome'))
setMethod('getGenome','qseaSet', function(object) 
        object@parameters$BSgenome )

setGeneric('getZygosity', function(object) 
        standardGeneric('getZygosity'))
setMethod('getZygosity','qseaSet', function(object) 
        object@zygosity )

setGeneric('setZygosity', function(object,...) 
        standardGeneric('setZygosity'))
setMethod('setZygosity','qseaSet', function(object, zygosityMatrix) {
        if(missing(zygosityMatrix))
            stop("provide a zygosity Matrix")
        #add checks of row and column names
        object@zygosity=zygosityMatrix 
        object
})

# parameters
setGeneric('getParameters', function(object) 
    standardGeneric('getParameters'))
setMethod('getParameters', 'qseaSet', function(object) 
    object@parameters)

setGeneric('hasEnrichment', function(object) 
    standardGeneric('hasEnrichment'))
setMethod('hasEnrichment', 'qseaSet', function(object) 
    !is.null(object@enrichment$parameters))

setGeneric('setEnrichment', function(object,...) 
    standardGeneric('setEnrichment'))
setMethod('setEnrichment', 'qseaSet', function(object, enrichment) {
    object@enrichment=enrichment
    object
})

setGeneric('getEnrichmentPattern', function(object) 
    standardGeneric('getEnrichmentPattern'))
setMethod('getEnrichmentPattern', 'qseaSet', function(object) 
    object@enrichment$pattern_name)

setGeneric('getEnrichmentDensity', function(object) 
    standardGeneric('getEnrichmentDensity'))
setMethod('getEnrichmentDensity', 'qseaSet', function(object) 
    object@enrichment$density)

setGeneric('getEnrichmentFactors', function(object,...) 
    standardGeneric('getEnrichmentFactors'))
setMethod('getEnrichmentFactors', 'qseaSet', 
    function(object, samples=getSampleNames(object), minN=0) 
        object@enrichment$factors[object@enrichment$n >= minN,
            samples,drop=FALSE] )

setGeneric('getEnrichmentParameters', function(object,...) 
    standardGeneric('getEnrichmentParameters'))
setMethod('getEnrichmentParameters', 'qseaSet', 
    function(object, samples=getSampleNames(object)) 
        object@enrichment$parameters[samples,,drop=FALSE] )


# short to get window size
setGeneric('getWindowSize', function(object) 
        standardGeneric('getWindowSize'))
setMethod('getWindowSize', 'qseaSet', function(object) 
        object@parameters$window_size) #todo:allow different window sizes?

# Regions
setGeneric('getRegions', function(object) standardGeneric('getRegions'))
setMethod('getRegions', 'qseaSet', function(object) object@regions)

setGeneric('addRegionsFeature', 
    function(object,...) standardGeneric('addRegionsFeature'))
setMethod('addRegionsFeature', 'qseaSet', function(object, name, feature){
    feature=data.frame(feature)
    names(feature)=name
    mcols(object@regions)=cbind(mcols(object@regions), feature)
    object
})


# Library Size
setGeneric('getLibSize', 
    function(object, ...) standardGeneric('getLibSize'))
setMethod('getLibSize', 'qseaSet', 
    function(object,samples=getSampleNames(object),normalized=TRUE) {
        if (normalized)
            factors=object@libraries[["file_name"]][samples,"library_factor"]
        else
            factors=1
        if(any(is.na(factors)))
            stop("No library factors found. Run addLibraryFactors first")
        object@libraries[["file_name"]][samples, "valid_fragments"]*factors
})

setGeneric('setLibrary', function(object,...) standardGeneric('setLibrary'))
setMethod('setLibrary', 'qseaSet', function(object, colName, libTable){
    object@libraries[[colName]]=libTable
    object
})

setGeneric('getLibrary', function(object,...) standardGeneric('getLibrary'))
setMethod('getLibrary', 'qseaSet', function(object, colName)
    object@libraries[[colName]]
)

# Offset
setGeneric('getOffset', function(object, ...) standardGeneric('getOffset'))
setMethod('getOffset', 'qseaSet', 
    function(object,samples=seq_len(nrow(object@sampleTable)), scale="rpkm") {
        if(scale=="rpkm"){
            object@libraries[["file_name"]][samples, "offset"]
        }else if(scale=="rpw"){#reads per "normal" window
            getOffset(object, samples, "rpkm")*
                getWindowSize(object)*10^-9*
                getLibSize(object,samples)
        }else if(scale=="fraction"){#estimated fraction of total reads
            getOffset(object, samples, "rpkm")*
                getWindowSize(object)*10^-9*
                length(getRegions(object))
        }else stop("unknown scale parameter. Use one of rpkm, rpw or fraction")
    })


normMethod<-function(methods=NULL,...){
    ownMethods=list(...)
    #check user defined methods    
    if(length(ownMethods)>0 && (is.null(names(ownMethods)) || 
        any(names(ownMethods) %in% "")) || any(duplicated(names(ownMethods))))
        stop("unnamed normalization method definition")    
    #"q" requires "enrichment"


    snm=list(#predefined normalization schemas
        reads="",
        counts="",
        beta=c("enrichment", "cnv", "library_size","region_length", 
            "preference", "offset", "zygosity"),
        logitbeta=c("logit","enrichment", "cnv", "library_size",
            "region_length", "preference", "offset", "zygosity"),
        betaMedian=c("enrichment", "cnv", "library_size","region_length", 
            "preference","q50", "offset", "zygosity"),
        betaLB=c("enrichment", "cnv", "library_size","region_length", 
            "preference","q2.5", "offset", "zygosity"),
        betaUB=c("enrichment", "cnv", "library_size","region_length", 
            "preference","q97.5", "offset", "zygosity"), 
        rpm="library_size",
        nrpm=c("cnv", "library_size", "preference", "zygosity"),
        rpkm=c("library_size", "region_length", "zygosity"),
        nrpkm=c("cnv", "library_size","region_length", "preference", 
            "zygosity"),
        lognrpkm=c("log", "cnv", "library_size","region_length", "preference", 
            "zygosity")
    )
    #valid normalization methods
    vnm=c("enrichment", "cnv", "zygosity", "library_size","region_length", 
            "preference",  "offset", "log", "logit","")#+ psC& qXX.YY (regexpr)
    numbered_vnm=c("psC", "q", "minCut", "maxCut", "scaleF")
        #issues: q>100 possible, negative values not possible
    if(missing(methods))
        methodList=list()
    else if(all(methods %in% names(snm) ))
        methodList=snm[methods]
    else stop("invalid normalization method:",
        paste(methods[! methods %in% names(snm)],collapse=", "))
    methodList=c(methodList,ownMethods)
    if(length(methodList)==0) stop("no methods defined")
    allMethods=unique(unlist(methodList))
    good= (allMethods %in% vnm)
    for(nV in numbered_vnm)
        good[regexpr(paste0(nV,"\\d+(\\.\\d+)?$"),allMethods) ==1 ]=TRUE
    if(!all(good))
        stop("invalid normalization method(s): ",paste(allMethods[! good],
            collapse=", "))
    structure(methodList, class = "normMethods")
}


# CNV
setGeneric('hasCNV', function(object) standardGeneric('hasCNV'))
setMethod('hasCNV', 'qseaSet', function(object) length(object@cnv)>0)

setGeneric('getCNV', function(object, ...) standardGeneric('getCNV'))
setMethod('getCNV', 'qseaSet', function(object,samples=getSampleNames(object) )
    if(hasCNV(object)) return (object@cnv[,samples]) else return(NULL) )

setGeneric('setCNV', function(object, ...) standardGeneric('setCNV'))
setMethod('setCNV', 'qseaSet', function(object, cnv){
    sn=names(mcols(cnv))
    if(!all(getSampleNames(object)==sn)) 
        stop("sample names of CNV and qseaSet do not match")
    if(hasCNV(object)) warning("overwriting CNV")
    object@cnv=cnv
    object
})

# Counts
setGeneric('getCounts', function(object,...) 
    standardGeneric('getCounts'))
setMethod('getCounts', 'qseaSet', function(object,samples=NULL, windows=NULL ){
        if(!is.null(samples) && !is.null(windows)) 
            return(object@count_matrix[windows,samples, drop=FALSE])
        if(!is.null(samples)) 
            return(object@count_matrix[,samples, drop=FALSE])
        if(!is.null(windows)) 
            return(object@count_matrix[windows,, drop=FALSE])
        object@count_matrix        
    }
)

setGeneric('setCounts', function(object,...) 
    standardGeneric('setCounts'))
setMethod('setCounts', 'qseaSet', function(object,count_matrix=matrix() ){
    object@count_matrix=count_matrix
    object
})

#### test object class definition ###
setClass(Class = 'qseaGLM', 
    representation = representation(
        fullModelDesign='matrix',
        contrast='list', 
        parameters='list',
        fullModel='list',
        windows='numeric'
        ),
    prototype = prototype(),
    validity = function(object){
        return(TRUE)
})

setGeneric('getReducedModel', 
    function(object, ...) standardGeneric('getReducedModel'))
setMethod('getReducedModel', 'qseaGLM', function(object, contrast=1) 
    object@contrast[[contrast]])

setGeneric('setReducedModel', 
    function(object, ...) standardGeneric('setReducedModel'))
setMethod('setReducedModel', 'qseaGLM', function(object, name, fit) {
    object@contrast[[name]]=fit
    object
})

setGeneric('getFullModel', 
    function(object, ...) standardGeneric('getFullModel'))
setMethod('getFullModel', 'qseaGLM', function(object) 
    object@fullModel)

setGeneric('getCoefNames', function(object) standardGeneric('getCoefNames'))
setMethod('getCoefNames', 'qseaGLM', function(object) 
    colnames(object@fullModelDesign))

setGeneric('getWindows', function(object) standardGeneric('getWindows'))
setMethod('getWindows', 'qseaGLM', function(object) 
    object@windows)

setGeneric('getContrastNames', 
    function(object) standardGeneric('getContrastNames'))
setMethod('getContrastNames', 'qseaGLM', function(object) 
    names(object@contrast))

setMethod('getSampleNames', 'qseaGLM', function(object) 
    rownames(object@fullModelDesign))

setGeneric('getDesignMatrix', 
    function(object) standardGeneric('getDesignMatrix'))
setMethod('getDesignMatrix', 'qseaGLM', function(object) 
    object@fullModelDesign)

setMethod('getParameters', 'qseaGLM', function(object) 
    object@parameters)

### show ###
setMethod('show', signature='qseaGLM', definition=function(object) {
    message("qsea GLM object")
    message("=======================================")
    message("full model: ",object@parameters$link,"(y/norm) ~ ", 
        paste(gsub("\\(|\\)","",colnames(object@fullModelDesign)), 
        collapse=" + "))
    if(length(getContrastNames(object))>0)
        message("contrasts: ", paste(getContrastNames(object), collapse=", "))
})

setClass(Class = 'qseaPCA', 
    representation = representation(
        svd='list',
        sample_names='character', 
        factor_names='character'
        ),
    prototype = prototype(),
    validity = function(object)
        return(TRUE)
)


### show ###
setMethod('show', signature='qseaPCA', definition=function(object) {
    message("qsea PCA object")
    message("=======================================")
    message("PCA from ",length(object@sample_names)," samples and ",
        length(object@factors_names), " genomic windows")
})

setGeneric('getSVD', function(object) standardGeneric('getSVD'))
setMethod('getSVD', 'qseaPCA', function(object) object@svd)

setMethod('getSampleNames', 'qseaPCA', function(object) object@sample_names)

setGeneric('getFactorNames',function(object) standardGeneric('getFactorNames'))
setMethod('getFactorNames', 'qseaPCA', function(object) object@factor_names)


