\name{bootLOPR}
\alias{bootLOPR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap Lasso OLS Partial Ridge
}
\description{
Combines residual (or paired) bootstrap Lasso+Partial Ridge and residual (or paired) bootstrap Lasso+OLS, and produces confidence intervals for regression coefficients.
}
\usage{
bootLOPR(x, y, lambda2, B = 500, type.boot = "residual", thres = 0.5, alpha = 0.05, 
         OLS = TRUE, cv.method = "cv", nfolds = 10, foldid, cv.OLS = TRUE, tau = 0, 
         parallel = FALSE, standardize = TRUE, intercept = TRUE, parallel.boot = 
         FALSE, ncores.boot = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

   \item{x}{
Input matrix as in glmnet, of dimension nobs x nvars; each row is an observation vector.
}
  \item{y}{
Response variable.
}
\item{lambda2}{
Tuning parameter in the Partial Ridge. If missing, lambda2 will be set to 1/nobs, where nobs is the number of observations.
}
  \item{B}{
Number of replications in the bootstrap -- default is 500.
}
  \item{type.boot}{
Bootstrap method which can take one of the following two values: "residual" or "paired". The default is residual.
}
\item{thres}{
 A threshold parameter. For the variables/predictors with selection probability (obtained by bootstrap) larger than thres, this function uses bootstrap Lasso+OLS (if OLS=TRUE) or bootstrap Lasso (if OLS=FALSE) to produce confidence intervals; while, for the variables/predictors with selection probability (obtained by bootstrap) smaller than thres, this function uses bootstrap Lasso+Partial Ridge to produce confidence intervals.
}
  \item{alpha}{
Significance level -- default is 0.05.
}
\item{OLS}{
If TRUE, this function uses Lasso+OLS estimator to compute the residuals for residual bootstrap Lasso+Partial Ridge; otherwise, it uses Lasso estimator to compute the residuals for residual bootstrap Lasso+Partial Ridge. The default value is TRUE. This argument can be ignored for paired bootstrap Lasso+Partial Ridge.
}
  \item{cv.method}{
The method used to select lambda in the Lasso -- can be cv, cv1se, and escv; the default is cv.
}
  \item{nfolds, foldid, cv.OLS, tau, parallel}{
Arguments that can be passed to escv.glmnet.
}
\item{standardize}{
Logical flag for x variable standardization, prior to fitting the model. Default is standardize=TRUE.
}
  \item{intercept}{
Should intercept be fitted (default is TRUE) or set to zero (FALSE).
}
\item{parallel.boot}{
If TRUE, use parallel foreach to run the bootstrap replication. Must register parallel before hand, such as doParallel or others. See the example below.
}
\item{ncores.boot}{
Number of cores used in the bootstrap replication.
}
\item{...}{
Other arguments that can be passed to glmnet.
}
}
\details{
The function combines the performance of bootstrap Lasso+Partial Ridge and bootstrap Lasso+OLS (if OLS=TRUE). For "large" regression coefficient in the sense that their selection probability (obtained by bootstrap) is larger than a threshold value (thres), it uses bootstrap Lasso+OLS to produce confidence intervals which can decrease the interval length ; while, for "small" regression coefficients meaning that their selection probability (obtained by bootstrap) is smaller than a threshold value (thres), it uses bootstrap Lasso+Partial Ridge to produce confidence intervals which can guarantee coverage. Note that there are two arguments related to parallel, "parallel" and "parallel.boot": "parallel" is used for parallel foreach in the escv.glmnet; while, "paralle.boot" is used for the parallel foreach in the bootstrap replication precodure.
}
\value{
A list consisting of the following elements is returned.

  \item{lambda.opt}{
The optimal value of lambda selected by cv/cv1se/escv.
}
\item{Beta}{
Lasso+OLS (if OLS=TRUE) or Lasso (if OLS=FALSE) estimate of the regression coefficients.
}
\item{Beta.LPR}{
Lasso+Partial Ridge estimate of the regression coefficients. 
}
\item{interval}{
A 2 by p matrix containing the bootstrap Lasso+OLS (if OLS=TRUE) or bootstrap Lasso (if OLS=FALSE) confidence intervals -- the first row is the lower bounds of the confidence intervals for each of the coefficients and the second row is the upper bounds of the confidence intervals.
}
\item{interval.LPR}{
A 2 by p matrix containing the bootstrap Lasso+Partial Ridge confidence intervals -- the first row is the lower bounds of the confidence intervals for each of the coefficients and the second row is the upper bounds of the confidence intervals. 
}
\item{interval.LOPR}{
A 2 by p matrix containing the combining confidence intervals of bootstrap Lasso+Partial Ridge and bootstrap Lasso+OLS (or bootstrap Lasso if OLS=FALSE) -- the first row is the lower bounds of the confidence intervals for each of the coefficients and the second row is the upper bounds of the confidence intervals.
}

}

\examples{
library("glmnet")
library("mvtnorm") 

## generate the data
set.seed(2015)
n <- 200      # number of obs
p <- 500
s <- 10
beta <- rep(0, p)
beta[1:s] <- runif(s, 1/3, 1)
x <- rmvnorm(n = n, mean = rep(0, p), method = "svd")
signal <- sqrt(mean((x \%*\% beta)^2))
sigma <- as.numeric(signal / sqrt(10))  # SNR=10
y <- x \%*\% beta + rnorm(n)

## residual bootstrap Lasso OLS + Partial Ridge
set.seed(0)
obj <- bootLOPR(x = x, y = y, B = 10)
# confidence interval
obj$interval
sum((obj$interval[1,]<=beta) & (obj$interval[2,]>=beta))

## using parallel in the bootstrap replication
#library("doParallel")
#registerDoParallel(2)
#set.seed(0)
#system.time(obj <- bootLOPR(x = x, y = y))
#system.time(obj <- bootLOPR(x = x, y = y, parallel.boot = TRUE, ncores.boot = 2))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.