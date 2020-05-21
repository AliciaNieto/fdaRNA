plot.order <-
function(ar,x=1:nrow(ar), tp='l', xl=range(x), yl=c(min(ar[x,]), max(ar[x,])), xlb='',ylb='',m='',...){
  sar = apply(ar, 2, sort); O=ncol(sar); 
  r=rainbow(O, s = 1, v = 1, start = 0, end = .68, alpha = 1)
  po=prof.funct(sar)[[1]]*O; 
  plot(x, sar[x,1],type=tp,xlim=xl,ylim=yl,xlab=xlb,ylab=ylb,main=m,...)
  for (j in 1:O){w=which(po==j); if(length(w)!=0) {matlines(x,sar[x,w], type='l',col=r[O-j+1])}}
}
