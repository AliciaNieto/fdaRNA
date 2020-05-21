outlier.detection <-
function(ar, nor=T, c.T=F, alpha=.01, ap=1, mn=1:ncol(ar), t.m="random",MS=F,SM=F, clases=1, Cla=list(1:ncol(ar)), ce=1, ms='Multidimensional Scaling', sm='Spectral Map',xlb=' ',ylb=' ',...){
  if(nor){ar=normalization(ar)}
  dimension=dim(ar); G<-dimension[1];

  if (c.T){mult.val=c.T}
  else{require(mvtnorm)
       constant=rep(0,ap);Co=cov(ar);
       if ((dimension[2] %% 2)==1){c3=0}else{c3=integer(0)}
      for (j in 1:ap){
        z3=c(prof.funct(rmvnorm(G,  sigma = Co))[[3]],c3)
        constant[j]=quantile(z3/median(z3),probs=1-alpha)
      }
       mult.val=median(constant)
    }
  
  do<-ar[,mn]
  O=ncol(do);  prof.funct(do)->z
  
  re=rank(z[[1]], ties.method=t.m); 
  f=floor(O/2);U=rep(0,f);D=rep(0,f); Distance=matrix(0,3,f); z3=z[[3]]; Distance[3,]=z3
  U[1]=which(re==1);D[1]=which(re==2); 
  
  for(i in 2:f){
    u=which(re==(2*i-1));
    d=which(re==(2*i));
    if(dist(t(do[,c(U[i-1],u)]))<dist(t(do[,c(D[i-1],u)]))){U[i]=u;D[i]=d}else{U[i]=d;D[i]=u}
    }
   Distance[1,]=mn[U];Distance[2,]=mn[D];
   mMed=mult.val*median(z3)
  result=which(z3>mMed)
  
  ds=matrix(0,G,O);
  if (length(z[[2]])>1){for(i in 1:O){ds[,i]=do[,i]-apply(do[,z[[2]]],1,mean)}}else {for(i in 1:O){ds[,i]=do[,i]-do[,z[[2]]]}}
  lr=length(result); out=integer(0)
  if(lr>0){
    for (i in 1:lr){
      if(norm(as.matrix(ds[,U[result[i]]]),'F')>norm(as.matrix(ds[,D[result[i]]]),'F'))
        {out=c(out,Distance[1,result[i]])}
      else {out=c(out,Distance[2,result[i]])}
    }
  }
  
  if(MS){
    rd = dist(t(do),diag=TRUE); rs = cmdscale(rd)
    plot(rs,type="n",xlab=xlb,ylab=ylb,main=ms,...)
    names = paste(mn,sep="")
    for(Nclases in 1:clases){
      text(rs[Cla[[Nclases]],], names[Cla[[Nclases]]], cex=ce,col=Nclases)}
  }
  
  if(SM){
    require(mpm)
    dos<-mpm(data.frame(do), logtrans = FALSE)
    dos$row.names <- 1:nrow(do)
    par(cex=0.00001)
    pdf(file = NULL) 
    plot(dos,scale="eigen",row.size=1)->zz
    dev.off()
    par(cex=.8)
    names = paste(mn[-1],sep="")
    plot(zz$Columns,cex=0,xlab=xlb,ylab=ylb,main=sm,...)
    for(Nclases in 1:clases){
      if(Nclases==1){cla.v=(Cla[[Nclases]]-1)[-1]; text(zz$Columns[cla.v,], names[cla.v], cex=ce,col=Nclases)}
      else {text(zz$Columns[Cla[[Nclases]]-1,], names[Cla[[Nclases]]-1], cex=ce,col=Nclases)}
    }
  }
  
  list(Distance.values=Distance, benchmark.outliers =mMed, id.out =out, constant.Tukey=mult.val)
}
