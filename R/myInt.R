#' Internal function used to compute a special case of integral faster
#'
#' @param fun basic function
#' @param K number of spline basis functions
#' @param tol error tolerance of the integral
#' @param start start point of the interval
#' @param end end point of the interval
#'
#' @author Hongming Pu \email{phmhappier@@163.com}
#'

myInt<-function(fun,K,tol=1e-6,start=0,end=1){
  sumFunAll<-function(x){
    f<-fun(x)
    if(is.vector(f) && K!=1){
      return(f%*%t(f))
    }
    return(t(f)%*%f)
  }
  N<-1
  h<-end-start
  x<-c(start,end)
  T1<-h/2*sumFunAll(x)
  S<-T1
  repeat{
    h<-h/2
    x<-start+(2*1:N-1)*h
    T2<-T1/2+h*sumFunAll(x)
    I<-(4*T2-T1)/3
    if(max(abs(I-S))<tol)break
    N<-2*N
    T1<-T2
    S<-I
  }
  return(I)
}
