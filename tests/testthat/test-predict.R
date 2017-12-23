context('prediction')
library(dplyr)
library(splines)
set.seed(1)
n.sub<-20
n.time=20
n.var<-3
fpc.rate<-0.2
n.total<-n.sub*n.time
train.rate<-0.5
n.train<-floor(n.total*train.rate)
n.test<-n.total-n.train
train.set<-sample(1:n.total,n.train)
test.set<-c(1:n.total)[-train.set]
f1<-function(x){return(sin(2*pi*x))}
f2<-function(x){return(cos(2*pi*x))}
f.fpc<-function(pa,x){return(fpc.rate*sin(2*pi*(x+pa)))}
data.time = runif(n.total,0,1)
data.var<-matrix(rnorm(n.total*n.var),nrow=n.total,ncol=n.var)
res1<-f1(data.time)*data.var[,1]+f2(data.time)*data.var[,2]
subs<-as.vector(rep(1,n.time)%*%t(c(1:n.sub)))
fpc.par<-as.vector(rep(1,n.time)%*%t(runif(n.sub)))
res.fpc<-f.fpc(fpc.par,data.time)
n1<-length(res1)
sd<-0
res.error<-rnorm(n1)*sd
res<-res1+res.fpc+res.error
data.fin<-cbind(res,subs,data.time,data.var)
vars<-paste0('V',1:n.var)
colnames(data.fin)<-c('res','sub','time',paste0('V',1:n.var))
data.fra<-as.data.frame(data.fin)
data.train<-data.fra[train.set,]
data.test<-data.fra[test.set,]
formula = as.formula( paste("res~", paste(vars, collapse = "+")) )


err.mean<-function(x){return(mean(x*x))}
test_that("predict works",{
  gt<-data.test[,'res']
  res<-vsflcm(formula,data = data.train,id.time='time',t.min=0,t.max=1,
              id.sub = 'sub',lambda=0,method.optim='BFGS',method.obj = 'nuclear',
              delta = 0.01,fpc.on = TRUE,lam.nuc = 1)
  pred1<-predict(res,data.test,sub.reg = FALSE)
  pred2<-predict(res,data.test,sub.reg=TRUE)
  expect_equal(err.mean(gt-pred1)<0.1,TRUE)
  expect_equal(err.mean(gt-pred2)<0.1,TRUE)
})
