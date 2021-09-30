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
tnl=function(k,n,l){
  if (any(k < l))
    stop(paste("k must be >= l", "\n",""))
  if (any(n < (2*l+1)))
    stop(paste("n must be >= (2*l+1)", "\n",""))
  x=NULL
  for(j in 0:n) {x=rbind(x,t(partitions::compositions(j,n)))}
  ss=0; prob=NULL

  for(v in 1:n) {

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

  prob=prob/choose(2*n,n)
  tt=sum(prob[1:k])
  result=list(prob=prob[k],tt=tt)
  return(result)
}

#' The Tnl simulation function.
#' @export
#' @param k Vector of quantiles.
#' @param n Sample size.
#' @param l The class parameter of Tnl.
#' @param trial The number of trials.
#' @description Perform simulated tnl.
#' @return list of simulated probability and distribution function.
#' @examples
#' \dontrun{You need to select all of \code{k}, \code{n} and \code{l}}
#' \dontrun{\code{trial}=100000 by default}
#'tnl.sim(k=2,n=7,l=1)
#'\dontrun{
#'$prob
#'[1] 0.15476
#'
#'$tt
#'[1] 0.23158
#'}
#' @import plyr
#' @import stats
tnl.sim=function(k,n,l,trial=100000){
  if (any(k < l))
    stop(paste("k must be >= l", "\n",""))
  if (any(n < (2*l+1)))
    stop(paste("n must be >= (2*l+1)", "\n", ""))
  stest= function(a,b,l) {
    a=sort(a)
    b=sort(b)
    t=0
    for(i in 1:l)         {if (b[i]<a[i+l]) t=t+1}
    for(i in (l+1):(n-l)) {if (b[i]>=a[i-l]&b[i]<a[i+l]) t=t+1}
    for(i in (n-l+1):n)   {if (b[i]>=a[i-l]) t=t+1}
    return(t)
  }
  x=y=NULL
  statistic=NULL

  for(i in 1:trial){

    x=rnorm(n);  y=rnorm(n)
    statistic[i]=stest(x,y,l)
  }

  statistict=(plyr::count(statistic)/trial)$freq
  statistict=c(rep(0,(l-1)),statistict)
  tt=sum(statistict[1:k])
  result=list(prob=statistict[k],tt=tt)
  return(result)
}


#' The distribution function
#' @export
#' @rdname tnl
#' @param k Vector of quantiles.
#' @param n Sample size.
#' @param l The class parameter of Tnl.
#' @param trial The number of trials.
#' @description Gives the distribution function
#' @return The value of the distribution against the specified \code{k}.
#' @examples
#' \dontrun{You need to select all of \code{k}, \code{n} and \code{l}}
#' \dontrun{\code{trial}=100000 by default}
#' ptnl(k=2,n=6,l=2,trial = 100000)
#' \dontrun{
#' $k
#' [1] 2
#'
#' $n
#' [1] 6
#'
#' $l
#' [1] 2
#'
#' $method
#' [1] "exact"
#'
#' $ptnl
#' [1] 0.03030303}
ptnl=function(k,n,l,trial = 100000){
  if (any(k < l))
    stop(paste("k must be >= l", "\n",""))
  if (any(n < (2*l+1)))
    stop(paste("n must be >= (2*l+1)", "\n",""))
  if (n<=10) {ptnl=tnl(k,n,l)$tt;method="exact"}
  else {ptnl=tnl.sim(k,n,l,trial)$tt;method="simulation"}
  result=list(k=k,n=n,l=l,method=method,ptnl=ptnl)
  return(result)
}

