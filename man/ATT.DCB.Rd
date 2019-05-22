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
Wangcheng Yan, Wenjun Zhou
}
\examples{
n <- 2000
p <- 50
q <- p*(p-1)/2
tau <- 5
nclust <- 10
beta <- 2 / (1:p) / sqrt(sum(1/(1:p)^2))
clust.ptreat <- rep(c(0.1, 0.9), nclust/2)
cluster.center <- 0.5 * matrix(rnorm(nclust * p), nclust, p)
cluster <- sample.int(nclust, n, replace = TRUE)
X <- cluster.center[cluster,] + matrix(rnorm(n * p), n, p)
T <- rbinom(n, 1, clust.ptreat[cluster])
Y <- 2 + X \%*\% beta + rnorm(n, 0, 1) + tau * T

fit <- ATT.DCB(Y,T,X)
fit$ATT
}