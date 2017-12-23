#' predict for a new dataset based on a fitted concurrent functional linear model
#'
#' @param obj object returned from \code{\link{vsflcm}}
#' @param data dataset which fitted values should be computed for
#' @param sub.reg Logical, whether the subjects should be recognized based on the knowledge from the training data
#'
#' @return a vector of fitted values
#'
#' @author Hongming Pu \email{phmhappier@@163.com}
#'
#' @importFrom stats predict
#' @import dplyr
#'
#' @export
#'

predict.flcmgl<-function(obj,data,sub.reg=TRUE){
  n<-nrow(data)
  ts<-(as.matrix(data[obj$id.time])-obj$t.min)/(obj$t.max-obj$t.min)

  bas<-obj$basis.fun(ts)
  bas.beta<-bas%*%obj$beta
  sd.mat<-rep(1,n) %*% t(obj$sds)
  values.pre= data %>% subset(select=obj$predictors) %>%
    '/'(sd.mat) %>%
    '*'(bas.beta) %>% apply(1, sum)

  fun<-function(x){
    temp<-which(x)
    if(length(temp==1)){
      return(temp)
    }else{return(obj$n.sub+1)}
  }
  if(sub.reg && obj$fpc.on){
      thetaCs<-matrix(0,nrow=obj$K,ncol = obj$n.sub+1)
      thetaCs[,1:(obj$n.sub)]<-obj$thetaCs
      values.resi.pre<-rep(0,length(values.pre))
      m1<-matrix(rep(data[,obj$id.sub],obj$n.sub),nrow=n,ncol = obj$n.sub)
      m2<-matrix(rep(obj$subs,n),ncol = obj$n.sub,nrow=n,byrow = TRUE)
      jud<-(m1==m2)
      inds<-apply(jud, 1, fun)
      values.resi.pre<-apply(t(bas)*thetaCs[,inds],2,sum)
    values.pre<-values.pre+values.resi.pre
  }
  return(values.pre)
  }
