#' Cross validation to choose lambda in vsflcm
#'
#' @param lambdas tuning parameter vector for \code{lambda} to choose from
#' @param seed.split seed to randomly split the data set into training set and test set
#' @param num.fold \code{num.fold}-fold cross validation
#'
#' @param formula formula for the regression. should have form \code{ y ~ V1 + V2 + ... + Vk}. Note: don't contain \code{id.time}
#' @param data  data frame
#' @param id.time  the variable that represents time
#' @param id.sub variable giving subject ID
#' @param t.min minimum value to be evaluated on the time domain. if `NULL`, taken to be minium observed value.
#' @param t.max maximum value to be evaluated on the time domain. if `NULL`, taken to be minium observed value.
#' @param K number of spline basis functions for coefficients and fpc
#' @param spline.fun spline basis functions. If not default, it should be a function taking a vector(corresponding to time) belonging to [0,1] as input and returning a \code{length(vector)*K} spline matrix where each row corresponds to a time point in the vector and each column corresponds to a spline basis function
#' @param K0 number of FPCs
#' @param delta Only useful when \code{method.obj}='nonconvex'. When the square of the norm of the coefficient is smaller than \code{delta}, L1 penalty will smoothly turned into L2 penalty
#' @param method.optim the optimization method to be used. For details and other options, please refer to \code{optim}
#' @param times Only useful when \code{method.obj}='nonconvex'. the algorithm randomly chooses \code{times} initial points to optimize and select the best
#' @param fpc.on Logical. Whether fpc will be used
#' @param lam.nuc Only Useful when \code{method.obj}='nuclear'. penalty factor for the nuclear norm
#' @param method.obj the methodology chosen. Should be 'nuclear' or 'nonconvex'. Typically, we recommend 'nuclear'
#'
#' @author Hongming Pu \email{phmhappier@@163.com}
#'
#' @export
#'

cv.vsflcm<-function(formula,data=NULL,id.time=NULL,
                    id.sub=NULL,t.min=NULL,t.max=NULL,
                    K=5,spline.fun='B-spline',lambdas=10^seq(-2,2,1),
                    K0=2,lam.nuc=10,method.obj=c('nuclear','nonconvex'),
                    delta=1e-1,method.optim=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                                              "Brent"),times=1,fpc.on=TRUE,seed.split=1,num.fold=5)
{
  err.mat = matrix(NA, nrow = length(lambdas), ncol = 5)
  id.out = as.character(formula)[2]
  set.seed(seed.split)
  subjs = unique(data[id.sub][,1])
  groups = sample(subjs, length(subjs)) %>%
    split(., ceiling(seq_along(.)/ (length(.)/ num.fold)  ))

  for (fold in 1:num.fold) {
    data.train = filter(data, !(data[id.sub][,1] %in% groups[[fold]]))
    data.test = filter(data, data[id.sub][,1] %in% groups[[fold]])
    for (i in 1:length(lambdas)) {
      res<-vsflcm(formula,data=data.train,id.time=id.time,
                  id.sub=id.sub,t.min=t.min,t.max=t.max,
                  K=K,spline.fun=spline.fun,lambda=lambdas[i],K0=K0,
                  lam.nuc=lam.nuc,method.obj=method.obj,
                  delta=delta,method.optim=method.optim,times=times,fpc.on=fpc.on)
      yHat = predict(res, data = data.test, sub.reg=FALSE)
      err.mat[i, fold] = mean((data.test[id.out]-yHat)^2, na.rm = TRUE)
    }
  }
  return(apply(err.mat, 1, mean))
}
