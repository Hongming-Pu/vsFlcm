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
#' @import dplyr
#'
#' @export
#'

predict.flcmgl<-function(obj,data,sub.reg=TRUE){
  n<-nrow(data)
  ts<-(as.matrix(data[obj$id.time])-obj$t.min)/(obj$t.max-obj$t.min)

  bas<-obj$basis.fun(ts)
  bas.beta<-bas%*%obj$beta
  mean.mat<-rep(1,n) %*% t(obj$means)
  sd.mat<-rep(1,n) %*% t(obj$sds)
  values.pre= data %>% subset(select=obj$predictors) %>%
    '-'(mean.mat) %>% '/'(sd.mat) %>%
    '*'(bas.beta) %>% apply(1, sum)
  if(sub.reg && obj$fpc.on){
    values.resi.pre<-rep(0,length(values.pre))
    for (i in c(1:n)) {
      temp<-which(obj$subs== data[i,obj$id.sub])
      if(length(temp==1)){##find this subject
        values.resi.pre[i]<-t(bas[i,])%*%res$theta%*%res$cs[((temp-1)*res$K0+1):(temp*res$K0)]
      }
    }
    values.pre<-values.pre+values.resi.pre
  }
  return(values.pre)
  }
