

getCoverage<-function(file_name=NULL,Regions=NULL, fragment_length=NULL,
    minMapQual=1, uniquePos=TRUE, paired=FALSE, CpGrange=NULL, CpGpos=NULL){

    if(!paired && is.null(fragment_length) )
        stop("for single end reads, please provide an estimate ",
            "for the fragment length")
    chr.select=as.vector(na.omit(mixedsort(levels(seqnames(Regions)) )))
    tmp=strsplit(file_name, ".", fixed=TRUE)[[1]]
    if (tmp[length(tmp)] %in% c("gz","zip","bz2")){
        tmp=tmp[-length(tmp)]
        tmp[length(tmp)]=paste("zipped",tmp[length(tmp)],sep="_")
    }
    ext=tmp[length(tmp)]
    if(ext %in% c("bam", "BAM")){
        batch=GRanges(seqlevels(Regions),IRanges(1, seqlengths(Regions)))
        batchNames=paste0(seqnames(batch),":",start(batch),"-",end(batch))
        batchOffset=c(findOverlaps(batch, Regions, select="first"), 
            length(Regions)+1)
        bamindex=file.exists(paste(file_name,".bai", sep=""))
        message("Reading ",ext," alignment ",file_name)
        if(paired){
            scanFlag = scanBamFlag(isPaired=TRUE, isProperPair=TRUE ,
                hasUnmappedMate=FALSE, isUnmappedQuery = FALSE, 
                isFirstMateRead = TRUE, isSecondMateRead = FALSE)
            fields=c("rname", "pos", "strand","qwidth","isize")
        }else{
            scanFlag = scanBamFlag(isUnmappedQuery = FALSE)        
            fields=c("rname", "pos", "strand", "qwidth")
        }
        if(bamindex){
            #read only selected chromosomes
            scanParam=ScanBamParam(flag=scanFlag, what = fields, which=batch, 
                mapqFilter=minMapQual)    
        }else{
            message("warning: reduced efficiency, no bam index found")
            #read all and selected chromosomes/regions later
            scanParam=ScanBamParam(flag=scanFlag, what = fields, 
                mapqFilter=minMapQual)
        }

        

        ReadsL = scanBam(file=file_name, param=scanParam)            
        if(!bamindex){#is very inefficient and should be avoided for big files
            fields=names(ReadsL[[1]])            
            nfields=length(fields)            
            ReadsL=unlist(lapply(ReadsL[[1]],split,  
                f=ReadsL[[1]]$rname, drop=TRUE),FALSE, FALSE)
            nchr=(length(ReadsL)/nfields)
            chr_read=sapply(ReadsL[seq_len(nchr)],
                function(x) as.character(x[1]))
            idx=lapply(which(chr_read %in% chr.select), seq, by=nchr, 
                length.out=nfields)
            names(idx)=chr_read[which(chr_read %in% chr.select)]
            idx=idx[chr.select]
            ReadsL=lapply(idx,function(x) ReadsL[x])
            ReadsL=lapply(ReadsL, 'names<-',fields)
         }
        
        totalR=vapply(ReadsL, function(x) length(x[[1]]), 0, USE.NAMES=FALSE)
        message("Number of imported sequencing fragments: ", sum(totalR))
        message("Calculating short read coverage for genome wide windows.")
        if(!is.null(CpGrange)){
            if(CpGrange[1]==CpGrange[2])
                message("selecting fragments containing ",CpGrange[1]," CpGs")
            else
                message("selecting fragments containing between ",CpGrange[1],
                    " and ",CpGrange[2]," CpGs")
        }
        
        processBam<-function(idx){        
            #message("processing ", batchNames[idx]," in process ",Sys.getpid())

            if(uniquePos){
                keep=!duplicated(as.data.frame(ReadsL[[idx]]))
            }else{keep=TRUE}

            revStr=ReadsL[[idx]]$strand=="-"   
            pos=ReadsL[[idx]]$pos 
            if(paired){
                isize=abs(ReadsL[[idx]]$isize)
                pos[revStr] =pos[revStr] + ReadsL[[idx]]$qwidth[revStr] -
                    isize[revStr]    
            }else{
                isize=fragment_length[1]
                pos[revStr] =pos[revStr] + ReadsL[[idx]]$qwidth[revStr] -
                    fragment_length[1]
            }
            rm(revStr)
            #now: pos is begin of fragment
            if(!is.null(CpGrange)){
                #select reads within CpG range
                nCpG=countOverlaps(GRanges(seqnames=ReadsL[[idx]]$rname, 
                    ranges=IRanges(start=pos,width=isize)),
                    CpGpos)
                keep=keep & nCpG >= CpGrange[1] & nCpG <= CpGrange[2]
                rm(nCpG)
            }
            if( totalR[idx] > 0 && sum(keep) > 0 ){
                if(paired){
                    isize=isize[keep]
                    fragment=c(mean(isize),var(isize), 
                        min(isize), max(isize) )                 
                }else{
                    isize=fragment_length[1]
                    fragment=NULL
                }
                #mass point of reads
                massP = GRanges(seqnames=ReadsL[[idx]]$rname[keep], 
                    ranges=IRanges(start=pos[keep]+round(isize/2), width=1))
            }else{
                # no valid reads for this chromosome
                return(list(counts=numeric(batchOffset[idx+1]-batchOffset[idx]), fragment=NULL))
            }
            return(list(counts=countOverlaps(
                Regions[batchOffset[idx]:(batchOffset[idx+1]-1)], massP), 
                fragment=fragment))
        }
        #BPPARAM_in= bpstart(BPPARAM_in)
        #workers=bpworkers(BPPARAM_in)
        #if(workers>1)
        #    message("using ",workers, " cores")
        #count_list=unlist(
        #    bplapply(seq_along(ReadsL), processBam,BPPARAM=BPPARAM_in ),
        #    FALSE, FALSE)
        #BPPARAM_in= bpstop(BPPARAM_in)
        count_list=unlist(
            lapply(seq_along(ReadsL), processBam),FALSE, FALSE)
        rm(ReadsL)
        countIdx=seq(1,by=2, length.out=length(batch))
        uniqueR=vapply( count_list[countIdx], sum,0)
        message("Number of selected sequencing fragments: ", sum(uniqueR))
        if(paired){
            n=sum(uniqueR)
            stats=count_list[countIdx+1]
            nN=uniqueR>0
            stats=as.data.frame(stats[nN], fix.empty.names=FALSE)
            stats[is.na(stats)]=0    #set variance for contigs with 1 read to 0        
            mu=sum(stats[1,]*uniqueR[nN]/n)
            sd=sqrt(1/(n-1)*(sum(uniqueR[nN]*(mu-stats[1,])^2) + 
                sum((uniqueR[nN]-1)*stats[2,])))
            message("Mean fragment length: ", 
                round(mu,3), " nt", sep="")
            message("SD of fragment length: ", 
                round(sd,3), " nt", sep="")
            message("Max fragment length: ",
                max(as.numeric(stats[4,])), " nt", sep="")
            message("Min fragment length: ",
                min(as.numeric(stats[3,])), " nt", sep="")

            return(list(
                counts=unlist(count_list[countIdx],FALSE, FALSE), 
                library=c(sum(totalR),n,NA,mu,sd,NA)))
        }else{       
            if(length(fragment_length)==1)
                fragment_length=c(1,0.1)*fragment_length
            return(list(
                counts=unlist(count_list[countIdx],FALSE, FALSE), 
                library=c(sum(totalR), sum(uniqueR),NA, 
                    fragment_length[1:2],NA)))
        }

        
    }else if (ext %in% c("zipped_wig", "zipped_bw", "wig", "bw")){
        message("Reading coverage file ",file_name)
        if(! is.null(CpGrange)){
            warning("ignoring CpG range parameter for coverage based ",
                "input file",immediate.=TRUE, call.=FALSE)
        }
        sel=GRanges(chr.select,IRanges(1, 536870912))#whole genome
            #this function will warn, if type is not bigwig
            wiggle=rtracklayer::import(file_name, asRangedData=FALSE, which=sel)
        m=as.matrix(findOverlaps(Regions, wiggle, minOverlap=width(Regions)/2))
        #check whether the regions overlap ambigiously
        coverage=numeric(length(Regions))        
        if(any(duplicated(m[,2])))
        {
            warning("Unsufficient resolution in coverage file ",file_name,
                ", returning values several times",immediate.=TRUE, call.=FALSE)
            #return(coverage) #return zeros
        }
        if(any(duplicated(m[,1]))) 
        {
            warning("Defind regions and coverage from ",file_name,
            " do not match uniquely, returning only the first matching value",
                immediate.=TRUE, call.=FALSE)
            #return(coverage) #return zeros
            m=matrix(c(seq_along(coverage),
                findOverlaps(Regions, wiggle,select="first")),
                nrow=length(coverage),2) #return the first match
        }
        if(any(is.na(rowSums(m))) || length(unique(m[,1])) < length(coverage))
        {
            warning("Not all regions are covered in ",file_name,
            ", returning zerro for those regions",immediate.=TRUE, call.=FALSE)
            m=m[!is.na(rowSums(m)),]
        }
        coverage[m[,1]]= mcols(wiggle)$score[m[,2]]
        return(list(counts=coverage,fragment=c(sum(coverage),NA,NA,NA) ))
    }else{
        warning("File format/extension of ",file_name," not supported")
        return(list(counts=numeric(length(Regions)),c(0,0,0,0) ))
    }
}


