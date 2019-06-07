#' @title Estimation of ATT
#'
#' @description Differentiated Confounder Balancing (DCB) Estimation of ATE
#'
#' @param outcome Y, binary treatment T, and covariates X
#'
#' @return ATT
#'
#' @examples ATT.DCB(Y,T,X)
#'
#' @export ATT.DCB


ATT.DCB <- function(Y,T,X,gp=TRUE,lambda=10,delta=0.001,mu=0.001,upsilon=0.001,thold=1e-4,max_iter=100000){
  M <- X
  
  mydata <- as.data.frame(cbind(Y,T,X))
  names(mydata)[1] <- "Y"
  names(mydata)[3:ncol(mydata)] <- paste0("X",1:ncol(X))
  fmla <- as.formula(paste("Y ~ ", paste(setdiff(names(mydata),c("Y","T")), collapse= " + ")))
  mod <- lm(formula=fmla, data=mydata[mydata$T==0,])
  gvmodel <- gvlma::gvlma(mod)
  if(gvmodel$GlobalTest$GlobalStat4$pvalue[1,1]<=0.05){
    for(i in 1:ncol(X)){
      for(j in i:ncol(X)){
        M <- cbind(M,X[,i]*X[,j])
      }
    }
    for(j in 1:ncol(M)){
      ms <- mean(M[,j])
      sd <- sqrt(var(M[,j]))
      M[,j] <- (M[,j] - ms)/sd
    }
  }
  
  n_c <- sum(T==0)
  p <- ncol(M)
  J <- c()
  values <- c()

  ####################################################################
  ##########################Given Variables###########################
  ####################################################################
  Y_c <- Y[T==0] - mean(Y[T==0])
  M_c <- M[T==0,]
  M_t <- M[T==1,]
  M_norm <- M
  for(j in 1:ncol(M)){
    ms <- mean(M_c[,j])
    M_c[,j] <- M_c[,j] - ms
    M_t[,j] <- M_t[,j] - ms
    M_norm[,j] <- M[,j] - ms
  }
  if(gp){
    M_t_bar <- apply(M_t, 2, mean)
  }
  else{
    M_t_bar <- apply(M_norm, 2, mean)
  }

  obj_func <- function(W,beta){
    diff <- M_t_bar - t(M_c) %*% W
    value1 <- (sum(beta*diff))^2 + sum((1+W)*(Y_c - M_c %*% beta)^2) + delta*sum(W^2)
    value2 <- mu*sum(beta^2) + upsilon*sum(abs(beta))
    return(value1 + value2)
  }

  ####################################################################
  #############################Initialize#############################
  ####################################################################
  W <- rep(1/n_c,n_c)
  beta <- rep(1/p,p)

  J[1] <- obj_func(W,beta)

  omega <- sqrt(W)

  fr <- function(omega){
    W <- omega^2
    diff <- M_t_bar - t(M_c) %*% W
    value <- (sum(beta*diff))^2 + sum((1+W)*(Y_c - M_c %*% beta)^2) + delta*sum(W^2)
    return(value)
  }

  grr <- function(omega){

    W <- omega^2
    diff <- M_t_bar - t(M_c) %*% W

    comp1 <- -4 * sum(beta*diff) * matrixcalc::hadamard.prod(M_c %*% beta,omega)
    comp2 <- 4 * delta * matrixcalc::hadamard.prod(W,omega)
    comp3 <- 2 * lambda * matrixcalc::hadamard.prod(omega,(Y_c-M_c %*% beta)^2)
    J_omega <- as.vector(comp1 + comp2 + comp3)
    return(J_omega)
  }

  ####line search
  fit <- rje::armijo(fun=fr, x=omega, dx = -grr(omega))
  eta <- fit$adj


  for(ind in 1:max_iter){
    ####################################################################
    #########################Fix W Update beta##########################
    ####################################################################
    diff <- as.vector(M_t_bar - t(M_c) %*% W)
    Y_prime <- sqrt(lambda*(1+W))*Y_c
    M_prime <- as.vector(sqrt(lambda*(1+W)))*M_c

    ####Elastic net####
    fit <- glmnet::glmnet(x=rbind(M_prime,diff), y=c(Y_prime,0),
                  family="gaussian", alpha=1/(1+2*mu/upsilon), lambda = 2*mu+upsilon)
    beta <- as.vector(fit$beta)

    values <- c(values,obj_func(W,beta))

    ####################################################################
    #########################Fix beta Update W##########################
    ####################################################################
    omega <- sqrt(W)

    fr <- function(omega){
      W <- omega^2
      diff <- M_t_bar - t(M_c) %*% W
      value <- (sum(beta*diff))^2 + sum((1+W)*(Y_c - M_c %*% beta)^2) + delta*sum(W^2)
      return(value)
    }

    grr <- function(omega){

      W <- omega^2
      diff <- M_t_bar - t(M_c) %*% W

      comp1 <- -4 * sum(beta*diff) * matrixcalc::hadamard.prod(M_c %*% beta,omega)
      comp2 <- 4 * delta * matrixcalc::hadamard.prod(W,omega)
      comp3 <- 2 * lambda * matrixcalc::hadamard.prod(omega,(Y_c-M_c %*% beta)^2)
      J_omega <- as.vector(comp1 + comp2 + comp3)
      return(J_omega)
    }

    ####line search
    #fit <- armijo(fun=fr, x=omega, dx = -grr(omega))
    #eta <- fit$adj
    #fit <- linesch_ww(fr,grr,omega,-grr(omega))
    #eta <- fit$alpha

    omega <- omega - eta * grr(omega)
    omega <- omega/sqrt(sum(matrixcalc::hadamard.prod(omega,omega)))
    W <- matrixcalc::hadamard.prod(omega,omega)
    #print(c(obj_func(W,beta),mean(Y[T==1])-sum(W*Y_c)))

    values <- c(values,obj_func(W,beta))

    ####################################################################
    #####################Objective function value#######################
    ####################################################################
    J[ind+1] <- obj_func(W,beta)
    if(abs(J[ind+1]-J[ind])<thold){break}
    #print(c(ind,(J[ind+1]-J[ind]),mean(Y[T==1]) - sum(W*Y[T==0])))
  }

  #plot(J[2:length(J)])

  out <- list("beta"=beta, "weight"=W, "ATT"=mean(Y[T==1]) - sum(W*Y[T==0]), "update"=values)
  return(out)
}

