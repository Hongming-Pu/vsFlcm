#' This function standarizes the predictors
#'
#' @param data data to be standarized
#'
#' @author Hongming Pu \email{phmhappier@@163.com}
#'
stanFit<-function(data){
  ms<-function(x){return(sum(x*x)/length(x))}
  dat.matrix<-as.matrix(data)
  sds<-sqrt(apply(dat.matrix,2,ms))
  n<-nrow(dat.matrix)
  m1<-matrix(1,nrow=n,ncol=1)
  dat.matrix<-dat.matrix/(m1%*%sds)
  return(list(mat=dat.matrix,sds=sds))
}
