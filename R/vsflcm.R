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
#' @param lambda Constant to multiply the L1 penalty term
#' @param K0 number of FPCs
#' @param delta Only useful when \code{method.obj}='nonconvex'. When the square of the norm of the coefficient is smaller than \code{delta}, L1 penalty will smoothly turned into L2 penalty
#' @param method.optim the optimization method to be used. For details and other options, please refer to \code{optim}
#' @param times Only useful when \code{method.obj}='nonconvex'. the algorithm randomly chooses \code{times} initial points to optimize and select the best
#' @param fpc.on Logical. Whether fpc will be used
#' @param lam.nuc Only Useful when \code{method.obj}='nuclear'. penalty factor for the nuclear norm
#' @param method.obj the methodology chosen. Should be 'nuclear' or 'nonconvex'. Typically, we recommend 'nuclear'
#' @param spline.fun2d second derivative of \code{spline.fun}, if \code{spline.fun} is default then the second derivative will be computed automaticly and you don't need to provide values for \code{spline.fun2d}. Otherwise you need to provide a function with input and output dimensions in accordance with \code{spline.fun} if you want to model the smoothness.
#' @param lam.smo penalty factor for the second derivative of basis funtion. Only useful when \code{spline.fun} is default or \code{spline.fun2d} is not NULL
#' @param maxit \code{maxit} of \code{optim}
#' @param intercept logical, whether use intercept function
#'
#' @author Hongming Pu \email{phmhappier@@163.com}
#'
#' @importFrom splines bs
#' @importFrom stats terms.formula quantile
#'
#' @export
#'

vsflcm<-function(formula,data=NULL,id.time=NULL, intercept=TRUE,
                 id.sub=NULL,t.min=NULL,t.max=NULL,
                 K=5,spline.fun='B-spline',lambda=0.5,K0=2,lam.nuc=10,method.obj=c('nuclear','nonconvex'),
                 delta=1e-1,method.optim=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                                           "Brent"),times=1,fpc.on=TRUE,spline.fun2d=NULL,lam.smo=0,maxit=100){

  if (is.null(id.sub)) {
    stop("Please specify the subject ID")}
  if(is.null(id.time)){
    stop("Please specify the time ID")}
  if(is.null(data)){
    stop('data missed, illegal input')}

  method.optim<-match.arg(method.optim)
  tf <- terms.formula(formula, specials = NULL)
  res<-list()
  res$predictors <- attr(tf, "term.labels")
  res$intercept=intercept
  if(intercept){
    res$predictors<-c(res$predictors,"intercept")
    data[,'intercept']<-1
  }
  res$method.obj<-match.arg(method.obj)
  variable.L<-as.character(formula)[2]
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
  res$sds<-stan.pre$sds
  num.pre<-length(res$sds)#number of predictors
  pre.mat<-stan.pre$mat
  ts<-as.matrix((data.sort[id.time]-t.min)/(t.max-t.min))
  ys<-as.matrix(data.sort[variable.L])
  res$y.mean<-mean(ys)
  ys<-ys-res$y.mean
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

  #compute basis and design matrix
  knots.bspline<-quantile(ts,probs=seq(0,1,length=K-2))[-c(1,K-2)]
  bspline.fun<-function(x){
    return(bs(c(-0.0001,1.0001,x),knots=knots.bspline,intercept=TRUE,degree=3)[-c(1,2),])
  }
  bspline.fun2d<-function(x){
    eps<-1e-5
    y1<-bspline.fun(x+eps)
    y2<-bspline.fun(x-eps)
    y<-bspline.fun(x)
    return((y1+y2-2*y)/(eps^2))
  }
  if(spline.fun!='B-spline'){
    basis.fun<-spline.fun
    bas.pen.mat<-myInt(basis.fun,K)
    if(! is.null(spline.fun2d)){
      bas.pen2d.mat<-myInt(spline.fun2d,K)
      bas.pen.mat<-bas.pen.mat+lam.smo*bas.pen2d.mat
    }
  }
  else{
    basis.fun<-bspline.fun
    bas.pen.mat<-myInt(bspline.fun,K)+lam.smo*myInt(bspline.fun2d,K)
  }
  res$basis.fun<-basis.fun
  bas.mat<-basis.fun(ts)
  design.mat<-kronecker(pre.mat,t(rep(1,K)))*kronecker(t(rep(1,num.pre)),bas.mat)

  obj.fun.raw1<-function(beta,thetaCs,fpc=TRUE){
    ys.hat<-design.mat%*%beta
    #penalty for beta
    beta.pen<-0
    for(i in c(1:num.pre)){
      beta.pen.i<-t(beta[((i-1)*K+1):(i*K)])%*%bas.pen.mat%*%beta[((i-1)*K+1):(i*K)]
      beta.pen<-beta.pen+sqrt(beta.pen.i)
    }
    beta.pen<-beta.pen*lambda

    #fpc residual
    ys.resi.hat<-rep(0,length(ys))
    for(i in c(1:res$n.sub)){
      ys.resi.hat[indexs[i]:(indexs[i+1]-1)]<-bas.mat[indexs[i]:(indexs[i+1]-1),]%*%thetaCs[,i]
    }
    ys.resi<-ys-ys.hat
    if(fpc){ys.resi<-ys.resi-ys.resi.hat}
    udv<-svd(thetaCs)
    return(sum(ys.resi*ys.resi)+beta.pen+lam.nuc*sum(udv$d))
  }

  gra.fun.raw1<-function(beta,thetaCs,fpc=TRUE){
    #gradient function
    ys.hat<-design.mat%*%beta
    ys.resi<-ys-ys.hat
    if(! fpc){
      der.beta<- t(ys.resi)%*%design.mat}
    der.thetaCs<-matrix(0,nrow=K,ncol=n.sub)

    for(i in c(1:res$n.sub)){
      a<-indexs[i]
      b<-indexs[i+1]-1
      ys.resi.hat<-bas.mat[a:b,]%*%thetaCs[,i]
      ys.resi[a:b]<-ys.resi[a:b]-ys.resi.hat
      der.thetaCs[,i]<-t(ys.resi[a:b])%*%bas.mat[a:b,]
    }

    if(fpc){der.beta<- t(ys.resi)%*%design.mat}
    der.beta<-der.beta*(-2)
    der.thetaCs<-der.thetaCs*(-2)
    udv<-svd(thetaCs)
    der.thetaCs<-der.thetaCs+lam.nuc*udv$u%*%diag(sign(udv$d))%*%t(udv$v)
    der.thetaCs<-as.vector(der.thetaCs)

    for(i in c(1:num.pre)){
      a<-(i-1)*K+1
      b<-i*K
      temp<-bas.pen.mat%*%beta[((i-1)*K+1):(i*K)]
      beta.pen.i<-sum(beta[((i-1)*K+1):(i*K)]*temp)
      if(beta.pen.i!=0){
        der.beta[a:b]<-der.beta[a:b]+lambda*temp/sqrt(beta.pen.i)}
    }
    if(! fpc){
      der.thetaCs<-rep(0,K*n.sub)}
    return(c(der.beta,der.thetaCs))
  }
  obj.fun.raw2<-function(beta,theta,cs,fpc=TRUE){
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
  gra.fun.raw2<-function(beta,theta,cs,fpc=TRUE){
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

  if(res$method.obj=='nuclear'){
    n1<-K*num.pre
    n2<-n1+K*n.sub
    n3=n2}
  else{
    n1<-K*num.pre
    n2<-n1+K*K0
    n3<-n2+K0*res$n.sub
  }
  pars.vec2split1<-function(pars.vec){
    #split the pars.vec
    beta<-pars.vec[1:n1]
    thetaCs<-pars.vec[(n1+1):n2]
    thetaCs<-matrix(thetaCs,nrow=K,ncol=n.sub)
    return(list(beta=beta,thetaCs=thetaCs))
  }
  pars.vec2split2<-function(pars.vec){
    #split the pars.vec
    beta<-pars.vec[1:n1]
    theta<-pars.vec[(n1+1):n2]
    cs<-pars.vec[(n2+1):n3]
    theta<-matrix(theta,nrow=K,ncol=K0)
    return(list(beta=beta,theta=theta,cs=cs))
  }
  if(res$method.obj=='nuclear'){
    pars.vec2split<-pars.vec2split1
    obj.fun.raw<-obj.fun.raw1
    gra.fun.raw<-gra.fun.raw1
  }else{
    pars.vec2split<-pars.vec2split2
    obj.fun.raw<-obj.fun.raw2
    gra.fun.raw<-gra.fun.raw2
  }
  obj.fun<-function(pars.vec){
    #pars.vec is a vector combining beta theta cs or beta thetaCs
    pars.split<-pars.vec2split(pars.vec)
    if(res$method.obj=='nuclear'){
      return(obj.fun.raw(pars.split$beta,pars.split$thetaCs,fpc = fpc.on))}
    else{
      return(obj.fun.raw(pars.split$beta,pars.split$theta,pars.split$cs,fpc = fpc.on))
    }
  }
  gra.fun<-function(pars.vec){
    pars.split<-pars.vec2split(pars.vec)
    if(res$method.obj=='nuclear'){
      return(gra.fun.raw(pars.split$beta,pars.split$thetaCs,fpc = fpc.on))}
    else{
      return(gra.fun.raw(pars.split$beta,pars.split$theta,pars.split$cs,fpc = fpc.on))
    }
  }
  pars.start<-rnorm(n3,0,1)

  cont<-list()
  cont$maxit<-maxit
  optim.res<-optim(par=pars.start,fn=obj.fun,gr=gra.fun, method = method.optim, control=cont)

  if(res$method.obj=='nonconvex'){
    for (i in c(2:times)) {
      pars.start<-rnorm(n3,0,i)
      optim.res.new<-optim(par=pars.start,fn=obj.fun,gr=gra.fun, method = method.optim)
      if(optim.res.new$value<optim.res$value){
        optim.res<-optim.res.new
      }
    }
  }


  thres<-1e-7
  for(i in 1:num.pre){
  temp<-bas.pen.mat%*%optim.res$par[((i-1)*K+1):(i*K)]
  beta.pen.i<-sum(optim.res$par[((i-1)*K+1):(i*K)]*temp)

    par2<-optim.res$par
    par2[((i-1)*K+1):(i*K)]<-0
    val.new<-obj.fun(par2)
    val.old<-obj.fun(optim.res$par)
    if(val.new<val.old || abs(val.new-val.old)/max(abs(val.new),abs(val.old))<thres){
      optim.res$par<-par2
    }
  }
  pars.list<-pars.vec2split(optim.res$par)
  res$beta<-matrix(pars.list$beta,nrow=K,ncol=length(res$predictors))


  if(fpc.on){
    if(res$method.obj=='nuclear'){
      res$thetaCs<-pars.list$thetaCs}
    else{
      theta<-pars.list$theta
      cs<-pars.list$cs
      cs<-matrix(cs,nrow=K0,res$n.sub)
      res$thetaCs<-theta%*%cs
    }
    udv<-svd(res$thetaCs)
    res$thetaD<-udv$d
  }
  else{
    res$thetaCs<-matrix(0,nrow=K,ncol=res$n.sub)
    res$thetaD<-rep(0,K)
  }
  res$fpc.on<-fpc.on
  res$id.sub<-id.sub
  res$K<-K
  res$K0<-K0
  res$lambda<-lambda
  res$lam.nuc<-lam.nuc
  res$lam.smo<-lam.smo
  res$opt<-method.optim
  indexs<-indexs[1:res$n.sub]
  res$subs<-data.sort[indexs,id.sub]
  beta<-res$beta

  class(res)="flcmgl"
  return(res)
}



