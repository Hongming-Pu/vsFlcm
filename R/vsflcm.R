#' vb_concurrent
#'
#' Implements group lasso with fpc for the linear functional concurrent model
#'
#' @param formula formula for the regression. should have form \code{ y ~ V1 + V2 + ... + Vk}. Note: don't contain \code{id.time}
#' @param data  data frame
#' @param id.time  the variable that represents time
#' @param id.sub variable giving subject ID
#' @param t.min minimum value to be evaluated on the time domain. if `NULL`, taken to be minium observed value.
#' @param t.max maximum value to be evaluated on the time domain. if `NULL`, taken to be minium observed value.
#' @param K number of spline basis functions for coefficients and fpc
#' @param spline.fun spline basis functions. If not default, it should be a function taking a vector(corresponding to time) belonging to [0,1] as input and returning a \code{length(vector)*K} spline matrix where each row corresponds to a time point in the vector and each column corresponds to a spline basis function
#' @param lambda Constant multiplying the L1 penalty term
#' @param K0 number of FPCs
#' @param delta When the square of the norm of the coefficient is smaller than \code{delta}, L1 penalty will smoothly turned into L2 penalty
#' @param method.optim the optimization method to be used. For details and other options, please refer to \code{optim}
#' @param times the algorithm randomly chooses \code{times} initial points to optimize and select the best
#' @param fpc.on Logical. Whether fpc will be used
#'
#' @author Hongming Pu \email{phmhappier@@163.com}
#'
#' @importFrom splines bs
#'
#' @export
#'

vsflcm<-function(formula,data=NULL,id.time=NULL,
                 id.sub=NULL,t.min=NULL,t.max=NULL,
                 K=5,spline.fun='B-spline',lambda=0.5,K0=2,
                 delta=1e-5,method.optim='BFGS',times=1,fpc.on=TRUE){
  if (is.null(id.sub)) {
    stop("Please specify the subject ID")}
  if(is.null(id.time)){
    stop("Please specify the time ID")}
  if(is.null(data)){
    stop('data missed, illegal input')}
  variable.L<-as.character(formula)[2]
  variable.R<-as.character(formula)[3]
  res<-list()
  variable.L<-gsub(" ","",variable.L)
  variable.R<-gsub(" ","",variable.R)
  res$predictors<-strsplit(variable.R,'+',fixed = TRUE)[[1]]
  datNamesFull<-colnames(data)
  if(FALSE %in% (res$predictors%in% datNamesFull)){
    stop('illegal input, please check the formula')}
  if(id.time %in% res$predictors){
    stop('id.time should not be in the formula')}
  if (! id.time %in%datNamesFull){
    stop('id.time not in the column names of data')}
  if(! id.sub %in% datNamesFull){
    stop('id.sub not in the column names of data')}
  if(spline.fun=='B-spline' && K<4){
    stop('If you choose the default B-spline, then you should choose K>=4')
  }
  if(is.null(t.max)){t.max<-max(data[id.time])}
  if(is.null(t.min)){t.min<-min(data[id.time])}
  res$t.min<-t.min
  res$t.max<-t.max
  res$id.time<-id.time

  data.sort<-data[order(data[id.sub]),c(variable.L,id.sub,id.time,res$predictors)]
  #standardize
  stan.pre<-stanFit(data.sort[res$predictors])
  res$means<-stan.pre$means
  res$sds<-stan.pre$sds
  num.pre<-length(res$means)#number of predictors
  pre.mat<-stan.pre$mat
  ts<-as.matrix((data.sort[id.time]-t.min)/(t.max-t.min))
  ys<-as.matrix(data.sort[variable.L])
  subs<-as.matrix(data.sort[id.sub])

  #assign the indexes according to the subjects
  res$n.sub<-1
  indexs<-c(1)
  res$n.total<-length(subs)
  if(res$n.total<=2){
    warning('too few rows')
  }
  for(i in c(2:res$n.total)){
    if(subs[i]!=subs[i-1]){
      res$n.sub<-res$n.sub+1
      indexs[res$n.sub]<-i
    }
  }
  indexs[res$n.sub+1]<-res$n.total+1

  #provide basis and design matrix
  knots.bspline<-quantile(ts,probs=seq(0,1,length=K-2))[-c(1,K-2)]
  basis.fun<-function(x){
  return(bs(c(0,1,x),knots=knots.bspline,intercept=TRUE,degree=3)[-c(1,2),])
  }
  if(spline.fun!='B-spline'){basis.fun<-spline.fun}
  res$basis.fun<-basis.fun
  bas.pen.mat<-myInt(basis.fun,K)
  bas.mat<-basis.fun(ts)
  design.mat<-kronecker(pre.mat,t(rep(1,K)))*kronecker(t(rep(1,num.pre)),bas.mat)

  obj.fun.raw<-function(beta,theta,cs,fpc=TRUE){
    #beta (K*num.pre) vector
    #theta K*K0 matrix
    #cs K0*res$n.sub vector
    ys.hat<-design.mat%*%beta
    #penalty for beta
    beta.pen<-0
    for(i in c(1:num.pre)){
      beta.pen.i<-t(beta[((i-1)*K+1):(i*K)])%*%bas.pen.mat%*%beta[((i-1)*K+1):(i*K)]
      if(beta.pen.i<delta){
        beta.pen<-beta.pen+sqrt(delta)/2+beta.pen.i/(2*sqrt(delta))
      }else{beta.pen<-beta.pen+sqrt(beta.pen.i)}
    }
    beta.pen<-beta.pen*lambda


    #fpc residual
    ys.resi.hat<-rep(0,length(ys))
    for(i in c(1:res$n.sub)){
      ys.resi.hat[indexs[i]:(indexs[i+1]-1)]<-bas.mat[indexs[i]:(indexs[i+1]-1),]%*%(theta%*%cs[((i-1)*K0+1):(i*K0)])
    }

    ys.resi<-ys-ys.hat
    if(fpc){ys.resi<-ys.resi-ys.resi.hat}
    return(sum(ys.resi*ys.resi)+beta.pen)
  }

  gra.fun.raw<-function(beta,theta,cs,fpc=TRUE){
    #gradient function
    ys.hat<-design.mat%*%beta
    ys.resi<-ys-ys.hat
    if(! fpc){
    der.beta<- t(ys.resi)%*%design.mat}
    der.cs<-rep(0,K0*res$n.sub)
    der.theta<-matrix(0,nrow=K,ncol=K0)

    for(i in c(1:res$n.sub)){
      a<-indexs[i]
      b<-indexs[i+1]-1
      a.cs<-(i-1)*K0+1
      b.cs<-i*K0
      temp.cs<-bas.mat[a:b,]%*%theta
      ys.resi[a:b]<-ys.resi[a:b]-temp.cs%*%cs[a.cs:b.cs]
      der.cs[a.cs:b.cs]<-t(ys.resi[a:b])%*%temp.cs
      der.theta<-der.theta+t(t(ys.resi[a:b])%*%bas.mat[a:b,])%*%t(cs[a.cs:b.cs])
    }

    #der.beta<- t(ys.resi)%*%design.mat
    if(! is.vector(der.theta)){der.theta<-as.vector(der.theta)}
    if(fpc){der.beta<- t(ys.resi)%*%design.mat}
    der.theta<-as.vector(der.theta)
    der.beta<-der.beta*(-2)
    der.theta<-der.theta*(-2)
    der.cs<-der.cs*(-2)

    for(i in c(1:num.pre)){
      a<-(i-1)*K+1
      b<-i*K
      temp<-bas.pen.mat%*%beta[((i-1)*K+1):(i*K)]
      beta.pen.i<-sum(beta[((i-1)*K+1):(i*K)]*temp)
      if(beta.pen.i<delta){
        der.beta[a:b]<-der.beta[a:b]+lambda*sqrt(beta.pen.i/delta)*temp
      }else{
      der.beta[a:b]<-der.beta[a:b]+lambda*temp/sqrt(beta.pen.i)}
    }
    if(! fpc){
    der.theta<-rep(0,n2-n1)
    der.cs<-rep(0,n3-n2)}
    return(c(der.beta,der.theta,der.cs))
  }

  n1<-K*num.pre
  n2<-n1+K*K0
  n3<-n2+K0*res$n.sub
  pars.vec2split<-function(pars.vec){
    #split the pars.vec into beta,theta,cs
    beta<-pars.vec[1:n1]
    theta<-pars.vec[(n1+1):n2]
    cs<-pars.vec[(n2+1):n3]
    theta<-matrix(theta,nrow=K,ncol=K0)
    return(list(beta=beta,theta=theta,cs=cs))
  }

  obj.fun<-function(pars.vec){
    #pars.vec is a vector combining beta theta cs
    pars.split<-pars.vec2split(pars.vec)
    return(obj.fun.raw(pars.split$beta,pars.split$theta,pars.split$cs,fpc = fpc.on))
  }
  gra.fun<-function(pars.vec){
    pars.split<-pars.vec2split(pars.vec)
    return(gra.fun.raw(pars.split$beta,pars.split$theta,pars.split$cs,fpc = fpc.on))
  }

  pars.start<-rnorm(n3,0,1)
  if(1==0){
  u<-c()
  for (i in c(1:n3)) {
    p<-pars.start
    eps<-1e-8
    p[i]<-p[i]-eps
    u[i]<-(obj.fun(pars.start)-obj.fun(p))/eps
  }
  print(u-gra.fun(pars.start)) }
  optim.res<-optim(par=pars.start,fn=obj.fun,gr=gra.fun, method = method.optim)
  for (i in c(2:times)) {
    pars.start<-rnorm(n3,0,i)
    optim.res.new<-optim(par=pars.start,fn=obj.fun,gr=gra.fun, method = method.optim)
    if(optim.res.new$value<optim.res$value){
      optim.res<-optim.res.new
    }
  }
  pars.list<-pars.vec2split(optim.res$par)

  res$predictors.sel<-c()
  res$beta<-matrix(pars.list$beta,nrow=K,ncol=length(res$predictors))
  if(fpc.on){
    res$theta<-pars.list$theta
    res$cs<-pars.list$cs
  }
  else{
    res$theta<-matrix(0,nrow=K,ncol=K0)
    res$cs<-rep(0,n3-n2)
  }
  res$fpc.on<-fpc.on
  res$id.sub<-id.sub
  res$K<-K
  res$K0<-K0
  indexs<-indexs[1:res$n.sub]
  res$subs<-data.sort[indexs,id.sub]
  beta<-res$beta
  n.sel<-0
  for(i in c(1:num.pre)){
    beta.pen.i<-t(beta[((i-1)*K+1):(i*K)])%*%bas.pen.mat%*%beta[((i-1)*K+1):(i*K)]
    if(beta.pen.i>=delta){
      n.sel<-n.sel+1
      res$predictors.sel[n.sel]<-res$predictors[i]
    }
  }
  class(res)='flcmgl'
  return(res)
}



