#' @title Estimation of ATE
#'
#' @description Differentiated Confounder Balancing (DCB) Estimation of ATE
#'
#' @param outcome Y, binary treatment T, and covariates X
#'
#' @return ATE
#'
#' @examples ATE.DCB_V1(Y,T,X)
#'
#' @export ATE.DCB_V1

ATE.DCB_V1 <- function(Y,T,X,ATT=FALSE,lambda=10,delta=0.001,mu=0.001,upsilon=0.001,thold=1e-4,max_iter=100000){
  if(ATT){
    fit <- ATT.DCB(Y,T,X,gp=ATT,lambda=lambda,delta=delta,mu=mu,upsilon=upsilon,thold=thold,max_iter=max_iter)
    wt <- rep(NA,length(T))
    wt[T==1] <- rep(1/sum(T==1),sum(T==1))
    wt[T==0] <- fit$weight
    out <- list("wt"=wt, "E[Y(1)]"=mean(Y[T==1]), "E[Y(0)]"=sum(fit$weight*Y[T==0]), "ATT"=fit$ATT)
  }
  else{
    fit1 <- ATT.DCB(Y,T,X,gp=TRUE,lambda=lambda,delta=delta,mu=mu,upsilon=upsilon,thold=thold,max_iter=max_iter)
    fit2 <- ATT.DCB(Y,1-T,X,gp=TRUE,lambda=lambda,delta=delta,mu=mu,upsilon=upsilon,thold=thold,max_iter=max_iter)
    #wt <- rep(NA,length(T))
    #wt[T==0] <- fit1$weight
    #wt[T==1] <- fit2$weight
    #out <- list("wt"=wt, "E[Y(1)]"=sum(wt[T==1]*Y[T==1]), "E[Y(0)]"=sum(wt[T==0]*Y[T==0]),
    #            "ATE"=sum(wt[T==1]*Y[T==1])-sum(wt[T==0]*Y[T==0]))
    out <- list("ATT"=fit1$ATT,"ATC"=-fit2$ATT,"ATE"=fit1$ATT*sum(T==1)/length(T)-fit2$ATT*sum(T==0)/length(T))
  }
  return(out)
}
