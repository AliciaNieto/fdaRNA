prof.funct <-
function(x){
  mg=as.matrix(dist(t(x))) 
  m=mg; k=ncol(x); for (i in 1:k) {m[i,i:k]=0}
  r=matrix(0,k,k); c=matrix(0,k,k); for (i in 1:k){r[i,]=i; c[,i]=i}
  p=rep(0,k); e=0
  if((k %% 2)==0) {fu=k-2; di=rep(0,k/2)} else {fu=k-1; di=rep(0,fu/2)} 
  cont=0
  for (i in which(((1:fu) %% 2)!=0)) 
  {cont=cont+1; di[cont]=max(m)
    w=which.max(m); ro=r[w][1]; co=c[w][1]; u=min(mg[-c(e,ro),ro]); du=min(mg[-c(e,co),co])
    if (u>du) 
    {p[ro]=i/k; p[co]=(i+1)/k} 
    else 
    {
      if (u==du) 
      {s=sample(c(1,0)); p[ro]=(i+s[1])/k; p[co]=(i+s[2])/k} 
      else 
      {p[ro]=(i+1)/k; p[co]=i/k}
    }
    m[ro,]=0; m[,co]=0; m[,ro]=0; m[co,]=0; e=c(e,ro,co)
  }
  p[which(p==0)]=1; if((k %% 2)==0){di[cont+1]=max(m)}
  list(depth.val = p, deepest.ele = which(p==max(p)), distant.val = di)
}
