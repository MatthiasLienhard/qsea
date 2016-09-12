context("glm")
library(MASS)
#test QSEAs matrix-version of glm.fit with stats::glm.fit

qs <- getExampleQseaSet()
sel=1:100
nm=qsea:::getNormMatrix(qs, methods=normMethod("beta")[[1]], windows=sel, 
samples=getSampleNames(qs))
design=model.matrix(~group, getSampleTable(qs))
y=getCounts(qs,getSampleNames(qs),sel)
fitM=qsea:::fitNBglmMatrix(design,y=y,disp=.1, nf=nm$factors,
    link="log", maxit=60, coef=NULL, eps=1e-11, fitDisp=FALSE)
fam=negative.binomial(link="log", theta=1/0.1)
fitR=list()
fitR$coefficients=matrix(NA, length(sel), 2)
fitR$deviance=numeric(length(sel))

for(i in 1:length(sel)){
    fitR_i=glm.fit(x=design, y=pmax(y[i,],1),
        family=fam, offset=log(nm$factors[i,]))
        fitR$deviance[i]=fitR_i$deviance
    fitR$coefficients[i,]=fitR_i$coefficients/log(2)
}

test_that("fitNBglmMatrix and glm.fit agree in deviance", {
    dev=abs(fitR$deviance-fitM$deviance)/(0.1+abs(fitR$deviance))>10e-6
    expect_false(any(dev))        
})

test_that("fitNBglmMatrix and glm.fit agree in coefficients", {
    cor1=cor(fitR$coefficients[,1], fitM$coefficients[,1])
    cor2=cor(fitR$coefficients[,2], fitM$coefficients[,2])
    expect_true(cor1>.999 & cor2>.999)        
})

