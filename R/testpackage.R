#' The (tnl) function.
#' @export
#' @param k Vector of quantiles.
#' @param n Sample size.
#' @param l The class parameter of Tnl.
#' @description Perform exact tnl.
#' @return list of an exact probability and distribution function.
#'
#'@examples
#'  ##You need to select all of \code{k}, \code{n} and \code{l}
#' tnl(k=2,n=7,l=1)
#'\dontrun{
#'$prob
#'[1] 0.1538462
#'
#'$tt
#'[1] 0.2307692
#'}
#'
#' @import partitions
#' @importFrom utils txtProgressBar setTxtProgressBar
tnl=function(k,n,l){
  if (any(k < l))
    stop(paste("k must be >= l", "\n",""))
  if (any(n < (2*l+1)))
    stop(paste("n must be >= (2*l+1)", "\n",""))
  x=NULL
  for(j in 0:n) {x=rbind(x,t(partitions::compositions(j,n)))}
  ss=0; prob=NULL
  pb=txtProgressBar(min=0,max=n,style=3,width=50)
  for(v in 1:n) {
    setTxtProgressBar(pb,v)
    ss=ss+1
    nn=nrow(x)
    count=0
    for (kk in 1:nn) {
      zz2=NULL
      m=x[kk,]
      for(i in 1:l) {if(sum(m[1:(i+l)])>=i) zz2[i]=1 else zz2[i]=0}
      for(i in (l+1):(n-l)) {
        if(sum(m[1:(i-l)])<i&sum(m[1:(i+l)])>=i) zz2[i]=1 else zz2[i]=0}
      for(i in (n-l+1):n) {if(sum(m[1:(i-l)])<i) zz2[i]=1 else zz2[i]=0}
      if(sum(zz2)==v) count=count+1
    }
    prob[ss]=count
  }
  close(pb)
  prob=prob/choose(2*n,n)
  tt=sum(prob[1:k])
  result=list(prob=prob[k],tt=tt)
  return(result)
}
