#' tune key parameters(K,lambda,lam.nuc,lam.smo) using cross validation in vsflcm, if you don't want to tune the parameter by yourself, this is for you.
#'
#' @param lam.smos tuning parameter vector for \code{lam.smo} to choose from
#' @param lambdas tuning parameter vector for \code{lambda} to choose from
#' @param lam.nucs tuning parameter vector for \code{lam.nuc} to choose from
#' @param seed.split seed to randomly split the data set into training set and test set
#' @param num.fold \code{num.fold}-fold cross validation
#' @param thre when automaticly chosing K, K will be chosen as the smallest one which ensures that the model can explain \code{thre} of the total variance.
#' @param formula formula for the regression. should have form \code{ y ~ V1 + V2 + ... + Vk}. Note: don't contain \code{id.time}
#' @param data  data frame
#' @param id.time  the variable that represents time
#' @param id.sub variable giving subject ID
#' @param t.min minimum value to be evaluated on the time domain. if `NULL`, taken to be minium observed value.
#' @param t.max maximum value to be evaluated on the time domain. if `NULL`, taken to be minium observed value.
#' @param K number of spline basis functions for coefficients and fpc. If default this function will automaticly choose K. Otherwise you need to provide a value.
#' @param spline.fun spline basis functions. If not default, it should be a function taking a vector(corresponding to time) belonging to [0,1] as input and returning a \code{length(vector)*K} spline matrix where each row corresponds to a time point in the vector and each column corresponds to a spline basis function
#' @param K0 number of FPCs
#' @param K.max max \code{K} can be chosen automaticly, only useful when \code{K} is NULL
#' @param delta Only useful when \code{method.obj}='nonconvex'. When the square of the norm of the coefficient is smaller than \code{delta}, L1 penalty will smoothly turned into L2 penalty
#' @param method.optim the optimization method to be used. For details and other options, please refer to \code{optim}
#' @param times Only useful when \code{method.obj}='nonconvex'. the algorithm randomly chooses \code{times} initial points to optimize and select the best
#' @param fpc.on Logical. Whether fpc will be used
#' @param method.obj the methodology chosen. Should be 'nuclear' or 'nonconvex'. Typically, we recommend 'nuclear'
#' @param intercept logical, whether intercept is used
#' @param spline.fun2d second derivative of \code{spline.fun}, if \code{spline.fun} is default then the second derivative will be computed automaticly and you don't need to provide values for \code{spline.fun2d}. Otherwise you need to provide a function with input and output dimensions in accordance with \code{spline.fun} if you want to model the smoothness.
#' @param maxit \code{maxit} of \code{optim}
#'
#' @author Hongming Pu \email{phmhappier@@163.com}
#'
#' @export
#'
tune.par.vsflcm<-function(formula,data=NULL,id.time=NULL, id.sub=NULL,t.min=NULL,t.max=NULL, seed.split=4, num.fold=3,thre=0.9,
         fpc.on=TRUE,intercept=TRUE, K=NULL,lambdas=10^seq(-2,2,1),lam.smos=c(0,10^seq(-2,2,1)),lam.nucs=c(0,10^seq(-2,2,1)),
         K0=2,spline.fun='B-spline',method.obj='nuclear', K.max=10,
         delta=1e-1,method.optim="BFGS",times=1,spline.fun2d=NULL,maxit=200)

  {

  id.out = as.character(formula)[2]
  set.seed(seed.split)
  subjs = unique(data[id.sub][,1])
  groups = sample(subjs, length(subjs)) %>%
    split(., ceiling(seq_along(.)/ (length(.)/ num.fold)  ))

  #choose K as the smallest one that satisfies the model can explain thre of variance
  if(is.null(K)){
  K.est<-4}
  else{K.est<-K}
  res<-vsflcm(formula,data=data,
                K=K.est,lambda=0,lam.smo=0,lam.nuc=0,
                spline.fun=spline.fun,K0=K0,method.obj=method.obj,id.sub=id.sub,t.min=t.min,t.max=t.max,id.time=id.time, intercept=intercept,
                delta=delta,method.optim=method.optim,times=times,fpc.on=fpc.on,spline.fun2d=spline.fun2d,maxit=maxit)
  ys<-as.matrix(data[,id.out])
  var.whole<-var(ys)
  err.model<-function(obj){
    yHat<-predict(obj, data = data, sub.reg=TRUE)
    return(sum((data[,id.out]-yHat)^2)/(length(ys)-1))
  }
  if(is.null(K)){

  while(err.model(res)>(1-thre)*var.whole && K.est<K.max){

    K.est<-K.est+1
    res<-vsflcm(formula,data=data,
                 K=K.est,lam.smo=0,lam.nuc=0,lambda=0,
                 spline.fun=spline.fun,K0=K0,method.obj=method.obj,id.sub=id.sub,t.min=t.min,t.max=t.max,id.time=id.time, intercept=intercept,
                 delta=delta,method.optim=method.optim,times=times,fpc.on=fpc.on,spline.fun2d=spline.fun2d,maxit=maxit)
  }
  cat("choose K as ",K.est,'\n')
  }

  err.cv<-function(lambda,lam.nuc,lam.smo){
    err.arr<-rep(NA,num.fold)
  for (fold in 1:num.fold) {
    data.train = filter(data, !(data[id.sub][,1] %in% groups[[fold]]))
    data.test = filter(data, data[id.sub][,1] %in% groups[[fold]])
    res.temp<-vsflcm(formula,data=data.train,
                      K=K.est,lambda=lambda,lam.nuc=lam.nuc,lam.smo=lam.smo,
                      spline.fun=spline.fun,K0=K0,method.obj=method.obj,id.sub=id.sub,t.min=t.min,t.max=t.max,id.time=id.time, intercept=intercept,
                      delta=delta,method.optim=method.optim,times=times,fpc.on=fpc.on,spline.fun2d=spline.fun2d,maxit=maxit)
      yHat = predict(res.temp, data = data.test, sub.reg=FALSE)
      err.arr[fold] = mean((data.test[id.out]-yHat)^2, na.rm = TRUE)
  }
    return(mean(err.arr))
  }

  lambda.est<-lambdas[1]
  lam.nuc.est<-lam.nucs[1]
  lam.smo.est<-lam.smos[1]
  err<-err.cv(lambda.est,lam.nuc.est,lam.smo.est)
  for(i in 1:length(lambdas)){
    lambda<-lambdas[i]
    err2<-err.cv(lambda,lam.nuc.est,lam.smo.est)
    if(err2<err){
      err<-err2
      lambda.est<-lambda
    }
  }
  cat('choose lambda as ',lambda.est,'\n')

  if(fpc.on && method.obj=='nuclear'){
  for(i in 1:length(lam.nucs)){
    lam.nuc<-lam.nucs[i]
    err2<-err.cv(lambda.est,lam.nuc,lam.smo.est)
    if(err2<err){
      err<-err2
      lam.nuc.est<-lam.nuc
    }
  }
  cat('choose lam.nuc as ',lam.nuc.est,'\n')
  }

  for(i in 1:length(lam.smos)){
    lam.smo<-lam.smos[i]
    err2<-err.cv(lambda.est,lam.nuc.est,lam.smo)
    if(err2<err){
      err<-err2
      lam.smo.est<-lam.smo
    }
  }
  cat('choose lam.smo as ',lam.smo.est,'\n')

  res<-vsflcm(formula,data=data,
                    K=K.est,lambda=lambda.est,lam.nuc=lam.nuc.est,lam.smo=lam.smo.est,
                    spline.fun=spline.fun,K0=K0,method.obj=method.obj,id.sub=id.sub,t.min=t.min,t.max=t.max,id.time=id.time, intercept=intercept,
                    delta=delta,method.optim=method.optim,times=times,fpc.on=fpc.on,spline.fun2d=spline.fun2d,maxit=maxit)
  return(res)
}
