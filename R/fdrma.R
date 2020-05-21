fdrma <-
function(slc0cel,method=c("fdn","myrma"), mod=c("MP","MM","LS"),no.cores = floor(detectCores())) {
  require(affy); require(MASS)
  if(method[1]=="myrma") { my.PM=pm(normalize.AffyBatch.quantiles(slc0cel)) } else
    if(method[1]=="fdn"){
      slc0celbg=bg.correct.rma(slc0cel)
      pmslc0 = pm(slc0celbg)
      #    pmslc0 <- apply(pmslc0, 2, bg.adjust)
      my.PM= normalization(pmslc0)
      z= floor (min((my.PM))); if(z<0) my.PM=my.PM-z
      #    for(i in 1:ncol(my.PM)) my.PM[,i]= my.PM[rank(pmslc0[,i]),i]
    }
  pG = probeNames(slc0cel)
  pn = unique(pG)
  
  G = length(pn)
  n = ncol(my.PM)
  xx = log2(my.PM)
  
  cl = makeCluster(no.cores)
  registerDoSNOW(cl)
  z = split(data.frame(xx),pG)
  fun1=function(x){ xxx = medpolish (as.matrix(x),trace.iter=F);xxx$overall+xxx$col }
  
  fun2=function(x) {
    i0 = nrow(x); nc=ncol(x)
    tt11 = data.frame( y=c(as.matrix(x)), x=factor(rep(colnames(x),rep(i0,nc) )),
                       p=factor(rep(1:i0,nc)))
    for(i1 in 1:1000) {
      tt4=try( rlm(y~p+x,data=tt11,method="MM",init="ls")$coef)
      if(is.numeric(tt4[1])) break
    }
    
    fun3=function(x) {
      i0 = nrow(x); nc=ncol(x)
      tt11 = data.frame( y=c(as.matrix(x)), x=factor(rep(colnames(x),rep(i0,nc) )),
                         p=factor(rep(1:i0,nc)))
      tt3 = lm(y~p+x,data=tt11)$coef
      tt3[1]+c(0,tt3[-(1:i0)])+mean(c(0,tt3[2:i0]))
    }
    
    tt4[1]+c(0,tt4[-(1:i0)])+median(c(0,tt4[2:i0]))
  }
  fun= if(mod=="MP") fun1 else if(mod=="MM") fun2 else fun3
  zz3=parSapply(cl,z,fun)
  stopCluster(cl)
  t(zz3)
}
