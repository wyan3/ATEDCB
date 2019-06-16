\name{ATT.DCB}
\alias{ATT.DCB}
\title{ATT Estimation using Differentiated Confounder Balancing (DCB)}
\usage{
ATT.DCB(Y,T,X)
}
\description{
The function provides the estimation of ATT given the outcome Y, binary treatment T, and covariates X.
}
\format{
  The function has the following parameters
  \describe{
    \item{\code{Y}}{Outcome}
    \item{\code{T}}{Binary treatment status indicator: 1=treatment; 0=control.}
    \item{\code{X}}{Covariate matrix}
    \item{\code{lambda}}{Penalty on the loss of outcome regression on covariates in the control group. Default is 10.}
    \item{\code{delta}}{Penalty on the L2-norm of weight for each unit. Default is 0.001}
    \item{\code{mu}}{Penalty on the L2-norm of coefficient for each covariate. Default is 0.001}
    \item{\code{upsilon}}{Penalty on the L1-norm of coefficient for each covariate. Default is 0.001}
    \item{\code{thold}}{Stop criteria of the update process. Default is 1e-4.}
    \item{\code{max_iter}}{Maximum number of iterations. Default is 100000.}
  }
}
\references{
Kuang, K., Cui, P., Li, B., Jiang, M., & Yang, S. (2017). Estimating treatment effect in the wild via differentiated confounder balancing. In Proceedings of the 23rd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (pp. 265-274).
}
\author{
}
\examples{
set.seed(25)
n <- 200
Z <- matrix(rnorm(4*n),ncol=4,nrow=n)
prop <- 1 / (1 + exp(Z[,1] - 0.5 * Z[,2] + 0.25*Z[,3] + 0.1 * Z[,4]))
treat <- rbinom(n, 1, prop)
Y <- 200 + 10*treat+ (1.5*treat-0.5)*(27.4*Z[,1] + 13.7*Z[,2] + 13.7*Z[,3] + 13.7*Z[,4]) + rnorm(n)
X <- cbind(exp(Z[,1])/2,Z[,2]/(1+exp(Z[,1])), (Z[,1]*Z[,3]/25+0.6)^3,(Z[,2]+Z[,4]+20)^2)

fit <- ATT.DCB(Y,T,X)
fit$ATT
}
