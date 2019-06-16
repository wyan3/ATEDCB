#' @title Estimation of ATE
#'
#' @description Differentiated Confounder Balancing (DCB) Estimation of ATE
#'
#' @param outcome Y, binary treatment T, and covariates X
#'
#' @return ATE
#'
#' @examples ATE.DCB(Y,T,X)
#'
#' @export ATE.DCB

ATE.DCB <- function(Y,T,X,ATT=FALSE,lambda=10,delta=0.001,mu=0.001,upsilon=0.001,thold=1e-4,max_iter=100000){
  
  for(j in 1:ncol(X)){
    if(shapiro.test(X[,j])$p.value<0.05){
      if(sum(X[,j]<=0)==0){
        X[,j] <- bestNormalize::boxcox(X[,j])$x.t
      }
      if(sum(X[,j]<=0)>0){
        X[,j] <- bestNormalize::yeojohnson(X[,j])$x.t
      }
    }
  }
  
  mydata <- as.data.frame(cbind(Y,T,X))
  names(mydata)[1:2] <- c("Y","T")
  names(mydata)[3:ncol(mydata)] <- paste0("X",1:ncol(X))
  fmla <- as.formula(paste("Y ~ ", paste(setdiff(names(mydata),c("Y","T")), collapse= " + ")))
  mod <- lm(formula=fmla, data=mydata[mydata$T==0,])
  gvmodel <- gvlma::gvlma(mod)
  value1 <- gvmodel$GlobalTest$GlobalStat4$pvalue[1,1]
  
  mydata <- as.data.frame(cbind(Y,1-T,X))
  names(mydata)[1:2] <- c("Y","T")
  names(mydata)[3:ncol(mydata)] <- paste0("X",1:ncol(X))
  fmla <- as.formula(paste("Y ~ ", paste(setdiff(names(mydata),c("Y","T")), collapse= " + ")))
  mod <- lm(formula=fmla, data=mydata[mydata$T==0,])
  gvmodel <- gvlma::gvlma(mod)
  value2 <- gvmodel$GlobalTest$GlobalStat4$pvalue[1,1]
  
  if(value1>=0.05 | value2>=0.05){
    M <- X
  }
  
  if(value1<0.05 & value2<0.05){
    M <- X
    for(i in 1:ncol(X)){
      for(j in i:ncol(X)){
        inter <- X[,i]*X[,j]
        if(shapiro.test(inter)$p.value<0.05){
          if(sum(inter<=0)==0){
            M <- cbind(M,bestNormalize::boxcox(inter)$x.t)
          }
          if(sum(inter<=0)>0){
            M <- cbind(M,bestNormalize::yeojohnson(inter)$x.t)
          }
        }
      }
    }
  }
  
  if(ATT){
    fit <- ATT.DCB(Y,T,M,gp=ATT,lambda=lambda,delta=delta,mu=mu,upsilon=upsilon,thold=thold,max_iter=max_iter)
    wt <- rep(NA,length(T))
    wt[T==1] <- rep(1/sum(T==1),sum(T==1))
    wt[T==0] <- fit$weight
    out <- list("wt"=wt, "E[Y(1)]"=mean(Y[T==1]), "E[Y(0)]"=sum(fit$weight*Y[T==0]), "ATT"=fit$ATT)
  }
  else{
    fit1 <- ATT.DCB(Y,T,M,gp=ATT,lambda=lambda,delta=delta,mu=mu,upsilon=upsilon,thold=thold,max_iter=max_iter)
    fit2 <- ATT.DCB(Y,1-T,M,gp=ATT,lambda=lambda,delta=delta,mu=mu,upsilon=upsilon,thold=thold,max_iter=max_iter)
    wt <- rep(NA,length(T))
    wt[T==0] <- fit1$weight
    wt[T==1] <- fit2$weight
    out <- list("wt"=wt, "E[Y(1)]"=sum(wt[T==1]*Y[T==1]), "E[Y(0)]"=sum(wt[T==0]*Y[T==0]),
                "ATE"=sum(wt[T==1]*Y[T==1])-sum(wt[T==0]*Y[T==0]))
  }
  return(out)
}
