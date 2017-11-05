#' This function standarizes the predictors
#'
#' @param data data to be standarized
#'
#' @author Hongming Pu \email{phmhappier@@163.com}
#'
stanFit<-function(data){
  dat.matrix<-as.matrix(data)
  means<-apply(dat.matrix,2,mean)
  sds<-sqrt(apply(dat.matrix,2,var))
  n<-nrow(dat.matrix)
  m1<-matrix(1,nrow=n,ncol=1)
  dat.matrix<-(dat.matrix-m1%*%means)/(m1%*%sds)
  return(list(mat=dat.matrix,means=means,sds=sds))
}
