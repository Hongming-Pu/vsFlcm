context('vsflcm')
library(dplyr)
library(splines)
set.seed(1)
n.sub<-20
n.time<-20
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
sd<-0.05
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

test_that("vsflcm(nuclear) works on toy example",{
  res<-vsflcm(formula,data = data.train,id.time='time',t.min=0,t.max=1,
            id.sub = 'sub',lambda=10,method.optim='BFGS',method.obj = 'nuclear',
            delta = 0.01,fpc.on = TRUE,lam.nuc = 1)
  expect_equal(err.mean(res$beta[,3]),0,tol=0.1)
  expect_equal(res$method.obj,'nuclear')
})
test_that("vsflcm(nonconvex) works on toy example",{
  res<-vsflcm(formula,data=data.train,id.time='time',t.min=0,t.max=1,
              id.sub='sub',lambda=10,method.optim = 'BFGS',method.obj = 'nonconvex',
              delta = 0.01,times=1,fpc.on = TRUE)
  expect_equal(err.mean(res$beta[,3]),0,tol=0.1)
  expect_equal(res$method.obj,'nonconvex')
})

test_that("vsflcm(nuclear) options work",{
  lambda=5.5
  lam.nuc=7.8
  k=8
  m='CG'
  res<-vsflcm(formula,data=data.train,id.time='time',t.min=0,t.max=1,K=k,
              id.sub='sub',lambda=lambda,method.optim = m,method.obj = 'nuclear',
              fpc.on = TRUE,lam.nuc = lam.nuc)
  expect_equal(lambda,res$lambda,tol=0.1)
  expect_equal(lam.nuc,res$lam.nuc,tol=0.1)
  expect_equal(k,res$K)
  expect_equal(m,res$opt)
  res<-vsflcm(formula,data=data.train,id.time='time',t.min=0,t.max=1,K=k,
              id.sub='sub',lambda=lambda,method.optim = m,method.obj = 'nuclear',
              fpc.on =FALSE,lam.nuc = lam.nuc)
  expect_equal(FALSE,res$fpc.on)
})
test_that("vsflcm(nonconvex) options work",{
  lambda=5.5
  k=8
  k0=2
  m='CG'
  res<-vsflcm(formula,data=data.train,id.time='time',t.min=0,t.max=1,K=k,K0=k0,
              id.sub='sub',lambda=lambda,method.optim = m,method.obj = 'nonconvex',
              fpc.on = TRUE)
  expect_equal(lambda,res$lambda,tol=0.1)
  expect_equal(k0,res$K0)
  expect_equal(k,res$K)
  expect_equal(m,res$opt)
  res<-vsflcm(formula,data=data.train,id.time='time',t.min=0,t.max=1,K=k,
              id.sub='sub',lambda=lambda,method.optim = m,method.obj = 'nonconvex',
              fpc.on = FALSE)
  expect_equal(FALSE,res$fpc.on)
})

test_that("check the class of the object vsflcm returns",{
  res<-vsflcm(formula,data=data.train,id.time='time',t.min=0,t.max=1,
              id.sub='sub',lambda=10,method.obj = 'nuclear',
              fpc.on = TRUE,lam.nuc = 5)
  expect_is(res,'flcmgl')
  res<-vsflcm(formula,data=data.train,id.time='time',t.min=0,t.max=1,
              id.sub='sub',lambda=10,method.obj = 'nonconvex',
              fpc.on = FALSE)
  expect_is(res,'flcmgl')
})
