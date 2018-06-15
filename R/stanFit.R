#' This function standarizes the predictors
#'
#' @param data data to be standarized
#'
#' @author Hongming Pu \email{phmhappier@@163.com}
#'
stanFit<-function(data, intercept){
  ms<-function(x){return(sum(x*x)/length(x))}
  dat.matrix<-as.matrix(data)
  means<-apply(dat.matrix,2,mean)
  #only when intercept=TRUE the data will be centralized
  if(intercept){means[length(means)]<-0}
  else{means<-rep(0,length(means))}
  n<-nrow(dat.matrix)
  m1<-matrix(1,nrow=n,ncol=1)
  dat.matrix<-dat.matrix-m1%*%means
  if(intercept){
  sds<-sqrt(apply(dat.matrix,2,ms))}
  else{sds<-matrix(1,nrow=1,ncol=ncol(dat.matrix))}


  dat.matrix<-dat.matrix/(m1%*%sds)
  return(list(mat=dat.matrix,sds=sds,means=means))
}
