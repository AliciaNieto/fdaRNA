normalization.boxplot <-
function(ar, namesN=T, before=F, after=T, par.r=F, par.c=F, arg1="Before", arg2="After"){
  if(namesN){names(ar)<-1:dim(ar)[2]}
  if(before){
    if(after){
      if(par.r){par(mfrow=c(1,2))}
      if(par.c){par(mfrow=c(2,1))}
    }
    boxplot(log(ar+1),main=arg1)}
  if(after){boxplot(log(normalization(ar)+1),main=arg2)}
  par(mfrow=c(1,1))
}
