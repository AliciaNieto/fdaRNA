normalization <-
function (x){
   sar = apply(x, 2, sort); pom=prof.funct(sar)[[2]];
  if (length(pom)>1) xxm=sort(rowMeans(x[,pom])) else xxm=sort(x[,pom])
  xr <- c(apply(x, 2, rank))
  array(approx(1:nrow(x), xxm, xr)$y, dim(x), dimnames(x))}
