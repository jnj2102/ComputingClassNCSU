##Jami Jackson
## HW 8

#Clean up workspace

rm(list=ls())

#
#' Density of Dirichlet-multinomial
#'
#' \code{ddirmult} evaluates the Dirichlet-multinomial density
#'
#' @param x an n-by-d matrix of counts
#' @param alpha an n-by-d or n-by-1 matrix of parameters
#' @param log log-density (TRUE) or density (FALSE)
#' @return an n vector of (log)-densities
#' @seealso 
ddirmult <- function(x, alpha, log = FALSE) {
  
  # check x is count data
  if (any(x < 0))
    stop("Error: x should be nonnegative count data")
  # check positivity of alpha
  if (any(alpha < 0))
    stop("Error: Dirichlet-multinomial parameter alpha should be nonnegative")
  
  # check dimensions of input arguments
  if (is.vector(alpha) && length(alpha) > 1) {
    if (is.vector(x) && length(x) > 1) {
      if (length(x) != length(alpha)) {
        stop("Error: sizes of x and alpha do not match.")
      } else {
        # expand x to a matrix of matching size with alpha
        x <- matrix(x, 1, length(x)) 
      }
    } else if (is.vector(x) && length(x) <= 1) {
      stop("Error: x can not be a scalar")
    }
    # expand alpha to a matrix of matching size with x
    alpha <- matrix(alpha, nrow = nrow(x), ncol = length(alpha), byrow = TRUE)
  }
  if (any(dim(alpha) != dim(x)))
    stop("Error: dimensions of alpha and x do not match")
  
  # compute log-densities
  alphaRowSums <- rowSums(alpha)
  xRowSums <- rowSums(x)
  # lgamma(0 + 0) - lgamma(0) will produce NaNs.
  # This fixes x_{ij} = alpha_{ij} = 0 cases.
  # Is there better way to deal with this?
  alpha[(x == 0) & (alpha == 0)] <- 1
  # assemble log-likelihood
  # lgamma(0) throws a lot of warnings but they are valid
  logl <- suppressWarnings(
    lfactorial(xRowSums) + rowSums(lgamma(x + alpha)) + 
      lgamma(alphaRowSums) - (rowSums(lfactorial(x)) + 
                                rowSums(lgamma(alpha)) + lgamma(alphaRowSums + xRowSums))
  )
  # Deal with alphaRowSums == 0 cases
  # fix the lgamma(0 + 0) - lgamma(0) produces NaN here
  logl[(xRowSums == 0) & (alphaRowSums == 0)] <- 0
  logl[(xRowSums > 0) & (alphaRowSums == 0)] <- - Inf
  
  # output
  if (log)
    return(logl)
  else
    return(exp(logl))
}



#' MLE of Dirichlet-multinomial distribution
#'
#' \code{dirmultfit} finds the MLE of Dirichlet-multinomial
#'
#' @param X an n-by-d matrix of counts
#' @param weights an n vector of observation weights
#' @param alpha0 an n vector of starting point
#' @param algorithm optimizatio algorithm: 'MM'
#' @param tolfun convergence tolerance in log-likelihood values
#' @param maxiters maximum number of iterations
#' @param display verbose mode (TRUE) or not (FALSE)
#' @return a list that contains
#'  alphaHat MLE
#'  se standard error of MLE
#'  maximum log-likelihood at MLE
#'  iterations number of iterations performed 
#'  gradient gradient at MLE
#'  obsinfo observed information at MLE
#'  obsinfoInv inverse of observed information at MLE
#' @seealso
dirmultfit <- function(X, weights = NULL, alpha0 = NULL,
                       tolfun = 1e-6, maxiters = 1000, display = FALSE) {
  
  # check observation weights
  if (!is.null(weights)) {
    weights <- weight[rowSums(X) != 0]
    if (any(weights < 0))
      stop("Error: observation weigths should be positive")
  }
  
  # remove data points with batch size 0
  rsum <- rowSums(X)
  if (any(rsum == 0)) { 
    rmv <- sum(rsum == 0)
    message(paste( "Warning: ", rmv,
                   " rows are removed because the row sums are 0"))
  }
  
  # remove the bins with no observations
  csum <- colSums(X)
  if (any(csum == 0)) {
    rmv <- sum(csum == 0)
    message(paste("Warning: ", rmv,
                  " columns are removed because the column sums are 0"))
  }
  
  # cleaned up data
  data <- X[rsum != 0, csum != 0]
  N <- nrow(data)       ## Sample size
  d <- ncol(data)           ## Number of parameters
  m <- rowSums(data)  ## batch sizes
  
  # set default obs weights to be all 1s
  if (is.null(weights)) 
    weights <- rep(1, N)
  weightsum = sum(weights)
  
  # set starting points
  if (is.null(alpha0)) {
    # method of moment estimate
    rho <- sum(colSums(weights * (data / m)^2) / (colSums(weights * data / m)))
    alpha0 <- as.vector(colSums(weights * data / m) * (d - rho) / (rho - 1) / N)
    alpha0[alpha0 <= 0] = 1e-6
  } else {
    alpha0 <- alpha0[csum != 0]
    if (!is.vector(alpha0) && length(alpha0) != d) {
      stop("Error: dimension of alpha0 does not match data")
    } else if (any(alpha0 <= 0)) {
      # user provided starting values
      stop("Error: starting values should be positive")
    }
  }
  
  # prepare some variables before looping
  alphaHat <- alpha0
  alphaSum <- sum(alphaHat)
  loglIter <- sum(weights * ddirmult(data, alpha0, log = TRUE))
  
  # Print initial log-like if requested 
  if (display)
    print(paste("iteration ", 1, ", logL = ", loglIter, sep=""))
  
  if (maxiters == 1) iter <- 1
  else {

# pre-compute some auxillary variables based on data
# scts: a matrix of size n-by-max(data)
# scts[j,k] is # data points with j-th component > k-1
# rcts: a vector of size max(m)
# rcts[k] is # data points with batch size > k-1
Sj <- function(xj, kvec, weight) 
  Sj <- colSums(weight * outer(xj, kvec, ">"))
scts <- apply(data, 2, Sj, kvec = 0:(max(data) - 1), weight = weights)
rcts <- Sj(m, kvec = 0:(max(m) - 1), weight = weights)

##----------------------------------------##
## MM loop
for (iter in 2:maxiters) {
  
  # MM update
  num <- colSums(scts / (outer(0:(max(data) - 1), alphaHat, "+")))
  den <- sum(rcts / (sum(alphaHat) + (0:(max(m) - 1))))
  alphaHat <- alphaHat * num / den
  
  # log-likelihood at new iterate
  loglOld <- loglIter
  loglIter <- sum(weights * ddirmult(data, alphaHat, log = TRUE))
  
  # Print the iterate log-like if requested
  if(display)
    print(paste("iteration ", iter, ", logL = ", loglIter, sep=""))
  
  # check convergence criterion
  if (abs(loglIter - loglOld) < tolfun * (abs(loglOld) + 1)) 
    break          
  
}
##----------------------------------------##
## End of MM loop

} else {
  stop("unknown algorithm option")
}
}

# score, i.e., gradient
alphaMat <- matrix(alphaHat, nrow = N, ncol = d, byrow = TRUE)
score <- 
  colSums(weights * digamma(data + alphaMat)) - 
  weightsum * (digamma(alphaHat) - digamma(alphaSum)) -
  sum(weights * digamma(alphaSum + m))
# diagonal part of the observed information matrix
obsinfoDvec <- weightsum * trigamma(alphaHat) - 
  colSums(weights * trigamma(data + alphaMat))
obsinfoDvecInv <- 1 / obsinfoDvec
# the constant c in the observed information matrix
obsinfoC <- weightsum * trigamma(alphaSum) - 
  sum(weights * trigamma(m + alphaSum))
# compute standard errors
obsinfo <- diag(obsinfoDvec) - obsinfoC
obsinfoInv <- diag(obsinfoDvecInv) + 
  outer(obsinfoDvecInv, obsinfoDvecInv) / (1 / obsinfoC - sum(obsinfoDvecInv))
se <- sqrt(diag(obsinfoInv))

# restore to original data size
if (any(csum == 0)) {
  colidx <- (csum != 0)
  # parameter estimate
  tmp <- alphaHat
  alphaHat <- rep(0, ncol(X))
  alphaHat[colidx] <- tmp
  # gradient/score vector
  tmp <- score
  score <- rep(0, ncol(X))
  score[colidx] <- tmp
  score[!colidx] <- weightsum * digamma(alphaSum) - 
    sum(weights * digamma(alphaSum + m))
  # obs info matrix
  tmp <- obsinfo
  obsinfo <- matrix(0, ncol(X), ncol(X))
  obsinfo[colidx, colidx] <- tmp
  obsinfo[!colidx, ] <- - obsinfoC
  obsinfo[, !colidx] <- - obsinfoC
  # inverse of obs info matrix
  tmp <- obsinfoInv
  obsinfoInv <- matrix(0, ncol(X), ncol(X))
  obsinfoInv[colidx, colidx] <- tmp
}

# output
list(estimate = alphaHat, se = se, maximum = loglIter, iterations = iter, 
     gradient = score, obsinfo = obsinfo, obsinfoInv = obsinfoInv)
}


############################################################################




#Write a function to generate S sets of simulated data
#this function will take in a fixed I and a fixed d and fixed S
#it will output one matrix of counts I x d and one matrix of parameters I x d

generate <- function(S, I, d) {
  
  dat.param <- NULL
  
  dat.counts <- NULL
  
  for (k in 1 : S) {

    multi.param <- matrix(runif(I * d, min = 0, max = 1), nrow = I, 
                        ncol = d)
  
  #I need to normalize the matrix of parameters 
  
  final.multi.param <- multi.param/rowSums(multi.param)
  
  dat.param <- rbind(dat.param, final.multi.param)
  
  
  matrix.counts <- t(apply(final.multi.param, 1, 
                           function(x) rmultinom(n = 1, size = N,
                                                 prob = x)))
  
  dat.counts <- rbind(dat.counts, matrix.counts)
  

}

#output
out <- list(parameters = dat.param, counts = dat.counts)
return(out)

}


#create the two factors 

#Number of populations

I <- list(10, 20, 50, 100)

#names(I) <- c("10", "20", "50", "100")


#Number of categories

d <- list(2, 5, 10, 20)

#names(d) <- c("2", "5", "10", "20")

#Different combinations for factors

#treat <- expand.grid(I, d) 

param <- list()

count <- list()

for (i in 1 : I[[i]]){
  for (j in d[[1]] : d[[j]]){
    set.seed((i - 1) * I + (j - 1) * d)
    
    param <- generate(S, I, d)$parameters
    
    count <- generate(S, I, d)$counts
    
   output <- list(list(param), list(count)) 
    

    
   # df <- do.call("rbind", param)
    
#    df2 <- do.call("rbind", count)
    
  }
}  
  



set.seed(4)

               
               
               
               
               
               
               
               

sim.matrices <- apply(treat.dataf, 1,
  function(x) treat.grid(I = treat.dataf$X1, d = treat.dataf$X1))
  

  sigma.a <- treat$Var1[[1]]
  P <- as.vector(treat$Var2[[1]]) #Extract sigma.a and npatterns
  sigma.e<-1              # fixed values for errors dist. as standard normal
  g<-rep(1:length(P), P)  # group index
  N<-length(g)            # number of total observations
  
  #Set upvectors for results:
  a<-rep(0,S)
  b<-rep(0,S)
  c<-rep(0,S)
  
  #Simulate data from a multivariate distribution:
  Z<-model.matrix(~factor(g)-1) 
  #Check dim
  mu<-rep(0,N)
  Sigma<-sigma.a*tcrossprod(Z,Z)+sigma.e*diag(N)
  y<-mvrnorm(n=S,mu=mu,Sigma=Sigma)
  
  #Calculate statistics:
  for(i in 1:S){
    datos<-data.frame(cbind(y[i,],g))
    names(datos)<-c("y","g")
    m0<-lm(y ~ 1, data=datos)
    m1<-lme(fixed=y~1,random=~1|g,data=datos,method='ML')
    m2<-lme(fixed=y~1,random=~1|g,data=datos,method='REML')
    # Compute Likelihood Ratio Test
    lrt<-(2*( m1$logLik - logLik(m0)[1] ))
    a[i]=1-pchisq(2*lrt,df=1)
    # Compute Exact Likelihood Ratio Test distribution
    b[i]<-exactLRT(m=m1,m0=m0, seed = 1800)$p
    # Compute Exact Restricted Likelihood Ratio Test distribution
    c[i]<-exactRLRT(m2,seed = 1800)$p
  }
  results<-as.matrix(cbind(a,b,c))
  return(results=results)
}

