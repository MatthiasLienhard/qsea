getExampleQseaSet<-function(CpG=TRUE,CNV=TRUE,repl=2, 
    doSampling=TRUE,enrichmentAnalysis=TRUE, expSamplingDepth=50000, 
    bgfraction=.1){
    n=2*repl

    ws=500
    sampleTable=data.frame(
        sample_name=paste0("Sim",seq_len(repl), rep(c("T", "N"),each=repl)),
        file_name=NA,
        group=rep(c("Tumor", "Normal"), each=repl), 
        stringsAsFactors=FALSE)
    rownames(sampleTable)=sampleTable$sample_name

    reg=GRanges("chr1", IRanges(seq(1,5000000,ws), width=ws),
        seqinfo=Seqinfo("chr1", 5000000, genome="Examplegenome"))
    zyM=matrix(2,nrow(sampleTable), 1, 
        dimnames=list(sampleTable$sample_name,"chr1" ))
    
    qs=new('qseaSet', sampleTable=sampleTable,
                regions=reg,
                zygosity=zyM,
                count_matrix=matrix(), 
                cnv=GRanges(), 
                parameters=list(window_size=ws, BSgenome="example_genome"), 
                libraries=list(), 
                enrichment=list())
    if(CpG)
        qs=addPatternDensity(qs,name="CpG", patternDensity=cpg)
    if(CNV){
        CNVws=500000
        CNVreg=GRanges("chr1", IRanges(seq(1,5000000,CNVws),width=CNVws))
        cnvT=c(0,0,1,1,-2,-1,-1,0,0,0)
        cnvN=rep(0,10)
        mcols(CNVreg)=matrix(c(rep(cnvT,repl), rep(cnvN, repl)),
            ncol=n, dimnames=list(NULL, getSampleNames(qs)))
        qs=addCNV(qs, cnv=CNVreg)
    }
    if(doSampling){
        expSamplingDepth=rep(expSamplingDepth, length.out=n)
        bgfraction=rep(bgfraction, length.out=n)
        par=c(1,8,6)
        gr=getSampleGroups(qs)
        qs=setLibrary(qs, "file_name", 
            libTable=data.frame(
                row.names=getSampleNames(qs),
                total_fragments=expSamplingDepth,
                valid_fragments=expSamplingDepth, 
                library_factor=1,
                fragment_length=300,
                fragment_sd=30,
                offset=0))
        nm=normMethod("nrpkm")[[1]]
        nfT=getNormMatrix(qs, methods=nm, windows=seq_along(reg), 
            samples=gr$Tumor)$factors
        nfN=getNormMatrix(qs, methods=nm, windows=seq_along(reg), 
            samples=gr$Normal)$factors
        
        
        cpg[is.na(cpg)]=0
        trueMethN[is.na(trueMethN)]=0
        trueMethT[is.na(trueMethT)]=0
        bgReads=cbind(nfT, nfN)
        bgReads=t(t(bgReads)/colSums(bgReads)*
                expSamplingDepth*bgfraction)
        
        signalT=.sigmF2(cpg,par[1], par[2], par[3] )*trueMethT*nfT
        signalN=.sigmF2(cpg,par[1], par[2], par[3] )*trueMethN*nfN
        
        meanVals=cbind(signalT, signalN)
        meanVals=t(t(meanVals)/colSums(meanVals)*
                expSamplingDepth*(1-bgfraction))+bgReads

        counts=matrix(
            rnbinom(n=n*length(reg), size=10,mu=meanVals),
            ncol=n,
            dimnames=list(NULL, getSampleNames(qs)))
        qs=setCounts(qs, count_matrix=counts)
        libTab=getLibrary(qs, "file_name")
        libTab$total_fragments=colSums(counts)
        libTab$valid_fragments=colSums(counts)
        qs=setLibrary(qs, "file_name", libTab)
        qs=addOffset(qs, "CpG",.5)
        if(enrichmentAnalysis){
            f=which(cpg>.5)
            signal=matrix(c(rep(trueMethT,repl), rep(trueMethN, repl))
                ,ncol=n)
            signal[signal<.8]=NA
            qs=addEnrichmentParameters(qs,pattern_name="CpG", 
                signal=signal[f,], windowIdx=f )
        }
    }
    qs
}

#cpg=getRegions(epitreatQS)$CpG_density[5001:15000]
#extab=makeTable(epitreatQS, norm_methods="beta", 
#    groupMeans=getSampleGroups(epitreatQS), keep=5001:15000)
#trueMethN=extab$Normal_beta_means
#trueMethT=extab$Tumour_beta_means
#devtools::use_data(cpg,trueMethN, trueMethT, internal = TRUE, overwrite=T)

