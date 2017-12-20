# vsFlcm： Variable Selection for FLCM(Functional Linear Concurrent Model)#

Author: Hongming Pu (phmhappier@163.com)

Version: 0.2.0

This package implements two statistical methods  for selecting variables in the functional linear concurrent model. These methods are described later in this article.

# Installation#

You can install this from GitHub with [devtools](https://github.com/hadley/devtools):

~~~R
library(devtools)
devtools::install_github("Hongming-Pu/vsFlcm")
~~~

# Example of Use#

The code below simulates a dataset under the functional linear concurrent model. For each of 30 subjects, observations of 30 predictor functions and a response function are observed over times between 0 and 1.  Different subjects have different individual effects.

~~~R
library(vsFlcm)
set.seed(1)

#generate parameters
n.sub<-30
n.time<-30
n.var<-30
fpc.rate<-0.2
n.total<-n.sub*n.time
train.rate<-0.5
n.train<-floor(n.total*train.rate)
n.test<-n.total-n.train
train.set<-sample(1:n.total,n.train)
test.set<-c(1:n.total)[-train.set]

#functions
f1<-function(x){return(sin(2*pi*x))}
f2<-function(x){return(cos(2*pi*x))}
f.fpc<-function(pa,x){return(fpc.rate*sin(2*pi*(x+pa)))}

#generate the data
data.time = runif(n.total,0,1)
data.var<-matrix(rnorm(n.total*n.var),nrow=n.total,ncol=n.var)
res1<-f1(data.time)*data.var[,1]+f2(data.time)*data.var[,2]
subs<-as.vector(rep(1,n.time)%*%t(c(1:n.sub)))
fpc.par<-as.vector(rep(1,n.time)%*%t(runif(n.sub)))
res.fpc<-f.fpc(fpc.par,data.time)
n1<-length(res1)
sd<-0.2
res.error<-rnorm(n1)*sd
res<-res1+res.fpc+res.error
data.fin<-cbind(res,subs,data.time,data.var)
vars<-paste0('V',1:n.var)
colnames(data.fin)<-c('res','sub','time',paste0('V',1:n.var))
data.fra<-as.data.frame(data.fin)
data.train<-data.fra[train.set,]
data.test<-data.fra[test.set,]
formula = as.formula( paste("res~", paste(vars, collapse = "+")) )

#fit the model
res<-vsflcm(formula,data = data.train,id.time='time',t.min=0,t.max=1,
              id.sub = 'sub',lambda=10,method.optim='BFGS',method.obj = 'nuclear',
              delta = 0.01,times = 1,fpc.on = TRUE,lam.nuc = 3)
pred<-predict(res,data.test,sub.reg = FALSE)
error<-pred-data.test$res
print(mean(error*error))
~~~



