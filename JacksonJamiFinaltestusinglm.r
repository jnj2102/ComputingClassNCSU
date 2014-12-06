##Jami Jackson
## HW 8

#Clean up workspace

rm(list=ls())


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
#' @param algorithm optimizatio algorithm: 'Newton' (default) | 'MM'
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
dirmultfit <- function(X, weights = NULL, alpha0 = NULL, algorithm = 'Newton',
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
    ##----------------------------------------##
    ## The Newton loop
    if (algorithm == "Newton") {
      
      # set backtrack max iterations
      backtrackMaxiters <- 10
      
      for (iter in 2:maxiters) {
        
        # score vector
        alphaMat <- matrix(alphaHat, nrow = N, ncol = d, byrow = TRUE)
        score <- colSums(weights * digamma(data + alphaMat)) - 
          weightsum * (digamma(alphaHat) - digamma(alphaSum)) -
          sum(weights * digamma(alphaSum + m))
        # observed info. matrix = diag(obsinfoDvec) - obsinfoC
        obsinfoDvec <- weightsum * trigamma(alphaHat) - 
          colSums(weights * trigamma(data + alphaMat))
        obsinfoDvecInv <- 1 / obsinfoDvec
        obsinfoC <- weightsum * trigamma(alphaSum) - 
          sum(weights * trigamma(m + alphaSum))
        # shrink c if necessary to make obs. info. pos def
        if (obsinfoC * sum(obsinfoDvecInv) >= 1) {
          if (display) print("shrink c")
          obsinfoC <- 0.95 / sum(obsinfoDvecInv)
        }
        
        # compute Newton direction
        newtondir <- obsinfoDvecInv * score
        newtondir <- newtondir + 
          (sum(newtondir) / (1 / obsinfoC - sum(obsinfoDvecInv))) * 
          obsinfoDvecInv
        
        # line search by step halving
        if (any(newtondir < 0)) {
          # make sure Newton iterate always lands within boundary
          stepsize <- min(- alphaHat[newtondir < 0] / newtondir[newtondir < 0])
          stepsize <- min(0.95 * stepsize, 1)
        } else {
          stepsize <- 1
        }
        for (btiter in 1:backtrackMaxiters) {
          alphaNew <- alphaHat + stepsize * newtondir
          loglNew <- sum(weights * ddirmult(data, alphaNew, log = TRUE))
          # line search successful if improving log-L
          if (loglNew > loglIter) break
          else if (btiter == backtrackMaxiters) {
            warning("line search failed")
          } else {
            if (display) print("step halving")
            stepsize <- stepsize / 2
          }
        }
        alphaHat <- alphaNew
        alphaSum <- sum(alphaHat)
        loglOld <- loglIter
        loglIter <- loglNew
        
        # Print the iterate log-like if requested 
        if (display)
          print(paste("iteration ", iter, ", logL = ", loglIter, sep=""))
        
        # check convergence criterion
        if (abs(loglIter - loglOld) < tolfun * (abs(loglOld) + 1)) break
      } 
      ##----------------------------------------##
      ## End of Newton loop
      
    } else if (algorithm == "MM") {
      
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

N <- 20

generate <- function(S, I, d) {
  
  dat.MLE <- NULL
  
  dat.EB <- NULL
  
  for (k in 1 : S) {

    multi.param <- matrix(runif(I * d, min = 0, max = 1), nrow = I, 
                        ncol = d)
  
  #I need to normalize the matrix of parameters.  This is my "TRUTH"
  
  final.multi.param <- multi.param/rowSums(multi.param)
  
  
  matrix.counts <- t(apply(final.multi.param, 1, 
                           function(x) rmultinom(n = 1, size = N,
                                                 prob = x)))
 
 #find the MLE
 
 pi.MLE <- matrix.counts/N
 
 #Estimate the Dich-Mult parameter
 
 dirmult.alpha <- dirmultfit(matrix.counts, algorithm = 'Newton',
                 tolfun = 1e-6, maxiters = 1000, display = FALSE)$estimate
 
 batchsize.alpha <- sum(dirmult.alpha)
 
 alpha.matrix <- t(replicate(dim(matrix.counts)[1], dirmult.alpha))
 
 #Get the empirical Bayes estimate
  
 pi.EB <- (matrix.counts + alpha.matrix) * 1/(N + batchsize.alpha)
 
 #Get the Estimation error for MLE
 
   ERR.MLE <- sum(rowSums(abs(pi.MLE - final.multi.param))) * 1/(2 * I)
 
 #Get the estimation error for EB
 
   ERR.EB <- sum(rowSums(abs(pi.EB - final.multi.param))) * 1/(2 * I)
 
 #Now I can combine each of these estimation errors for each replicate
 #because I want to summarize the estimation errors, not anything else.
 #I'm left with one big column vector after S replicates for ERR MLE and 
 #ERR EB
 
  dat.EB <- rbind(dat.EB, ERR.EB)
  dat.MLE <- rbind(dat.MLE, ERR.MLE)
 
  
print(k)
}


# a quick t test

paired.t  <- t.test(dat.MLE, dat.EB, paired = TRUE)


#find the mean estimation error for EB and MLE

mean.err.EB <- mean(dat.EB)

mean.err.MLE <- mean(dat.MLE)

#find the standard error of the mean estimation error for EB and MLE

se.mean.err.EB <- sqrt(var(dat.EB)/S)

se.mean.err.MLE <- sqrt(var(dat.MLE)/S)


#output the means and standard errors because this is my function
list(mean.EB = mean.err.EB, mean.MLE = mean.err.MLE, 
     se.EB = se.mean.err.EB, se.MLE = se.mean.err.MLE,
	t.test = paired.t)


}


#create the two factors 

#Number of populations

I <- list(10, 20, 50, 100)

names(I) <- c("10", "20", "50", "100")


#Number of categories

d <- list(2, 5, 10, 20)

names(d) <- c("2", "5", "10", "20")


#param <- list()

#count <- list()

#Code to make sure no combo has the same seed number.  There are no 
#duplicates
# 
# for (i in names(I)){
#        for (j in names(d)){
#              print(((I[[i]] - 1) * 3 + (d[[j]] - 1) * 11))
#          }
# }
# 
#  for (i in names(I)){
#     
#               print(I[[i]])
# 
#  }
# 
# 
#   for (j in names(d)){
# 
#      print(d[[j]])
#     }
# 

#biglist.param <- list()

#biglist.count <- list()

biglist.mean.MLE <- list()

biglist.mean.EB <- list()

biglist.ttest <- list()

biglist.se.MLE <- list()

biglist.se.EB <- list()

#MLE <- list()

#MLE.var <- list()

#MLE.se <- list()



#Figure out what number the S should be based on the standard errors of the
#MLE by running the below code with different values of S

# S <- 3
# S <- 100

#Based on all standard errors being <= 0.005, I choose S = 100

S <- 100

#library(matrixStats)

# Start the clock!
ptm <- proc.time()

for (i in names(I)){
  for (j in names(d)){
 #   print(typeof((I[[i]] - 1) * 3 + (d[[j]] - 1) * 11))
    set.seed((I[[i]] - 1) * 3 + (d[[j]] - 1) * 11)
    
    #generate one treatment combo

    one.combo <- generate(S, I[[i]], d[[j]])
    
    meanEB <- list(one.combo$mean.EB)
    
    meanMLE <- list(one.combo$mean.MLE)

   ttest <- list(one.combo$t.test)
    
    name <- paste('item:', (I[[i]] - 1) * 3 + (d[[j]] - 1) * 11, sep = '')
    
    biglist.mean.EB[[name]] <- meanEB
    
    biglist.mean.MLE[[name]] <- meanMLE

 biglist.ttest[[name]] <- ttest
    
    seMLE <- list(one.combo$se.MLE)
    
    biglist.se.MLE[[name]] <- seMLE
    
    seEB <- list(one.combo$se.EB)
    
    biglist.se.EB[[name]] <- seEB
    
  #  df.MLE.count <- as.matrix(data.frame(matrix(unlist(biglist.count[name]), 
                      #            nrow = S * I[[i]], byrow = F)))
    
   # MLE[name] <- list(colMeans(df.MLE.count)/N)
    
    #find the sample variance of the MLEs
    
    #MLE.var <- colVars(df.MLE.count/N)
    
    #MLE.se[name] <- list(sqrt(MLE.var/S))
     
    #output()
  biglist.mean.EB
  biglist.mean.MLE
  biglist.se.EB
  biglist.se.MLE
 biglist.ttest

    print(j)
    

  }

 print(i)
}  
  
# Stop the clock
proc.time() - ptm


#double check that I can reproduce my results
#I can I run it exactly as my code is
#   i <- "10"
# I[[i]]
#  j <- "2"
# d[[j]]
#  set.seed((I[[i]] - 1) * 3 + (d[[j]] - 1) * 11)
# (I[[i]] - 1) * 3 + (d[[j]] - 1) * 11
#  38
# param <- list(generate(S, I[[i]], d[[j]])$parameters)
#  param
#  count <- list(generate(S, I[[i]], d[[j]])$counts)
#  count


#create a (fake) function to calculate the MLE of the multinomial, the mean MLE
#the SE of MLE and the ERR of MLE. It takes in the list, converts it to 
#a dataframe, and then does the calculations. I'm adding this to my
#generate function

# MLE <- function(S, biglist.count["item:38"], I[[i]], N) {
#                
# df.MLE.count <- as.matrix(data.frame(matrix(unlist(biglist.count["item:38"]), 
#                           nrow = S * I[[i]], byrow = F)))
# 
# MLE.d <- colMeans(df.MLE.count)/N

#  a <-  split(df.MLE.count, col(df.MLE.count))
#  
#  n <- dim(df.MLE.count)[2]
# 
# lhs  <- paste("cat.vec", 1 : n, sep = "")
#  
# rhs  <- paste("a[[",1 : n,"]]", sep = "")
# 
# eq   <- paste(paste(lhs, rhs, sep = "<-"), collapse = ";")
# 
# eval(parse(text = eq))
 
#playing around

# vec <- rep(1, 30)
# lh <- paste("cbind(vec")
# rh <- paste("a[[",1:n,"]])", sep = "")
# eqn <- paste(paste(lh, rh, sep = ","), collapse = ";")
# eval(parse(text = eqn))

  
  #MLE.i <- rowSums(df.MLE.count)/N       
  
               
               
#}

#############################################################################
             
#Do Question 3


#create the two factors 

#Sample Size

n <- list(50, 100, 200, 400)

names(n) <- c("50", "100", "200", "400")

library(MASS)
library(glmnet)

#Correlation

rho <- list(0, 0.25, 0.5, 0.9)

names(rho) <- c("0", "0.25", "0.5", "0.9")

#p the number of predictors is fixed 

p <- 50

#My truth beta estimates are fixed since I fix p

# one intercept and p-1 predictors

beta.truth <- rep(1, 50)

#I am keeping my noise factor fixed and constant 

noise <- 1


#Create a function to find the variance-covariance matrix
#that is based on rho

autocorr.mat <- function(n, rho) {
  mat <- diag(n)
  return(rho ^ abs(row(mat) - col(mat)))
  
}

#create a function that finds the OLS, ridge, and lasso estimates and
#calculates the MSE and prediction errors

library(mvtnorm)



#Figure out what number the S should be based on the standard errors of the
#MLE by running the below code with different values of S

# S <- 3
 S <- 500

#Based on all standard errors being <= 0.005, I choose S = 


generate.shrink <- function(S, n, rho) {
  
  dat.ls.mse <- NULL
  
  
#   dat.ls.mse1 <- NULL
  
  dat.ridge.mse <- NULL
  
  dat.lasso.mse <- NULL
  
  dat.ls.pe <- NULL
  
#   dat.ls.pe1 <- NULL
#   
  dat.ridge.pe <- NULL

  dat.lasso.pe <- NULL
  
  for (k in 1 : S) {
    
    #here is how I can generate the X matrix using random numbers
    
    #this version should also be random
     x.train.m <- mvrnorm(n = n, mu = rep(0, p - 1), 
                        Sigma = autocorr.mat(p - 1, rho), empirical = TRUE)
    
#           x.train.m <- rmvnorm(n = n, mean = rep(0, p - 1), 
#                               sigma = autocorr.mat(p - 1, rho))
    
    
#for some reason, I am getting singular matrices when I ask for n 
    #observations using either mvrnorm or rmvnorm.  I'm going to make
    # a loop of size n to get 1 observation n times
    
    x.train.m <- NULL
    
    for (u in 1 : n) {
    x.train.loop <- rmvnorm(n = 1, mean = rep(0, p - 1), 
              sigma = autocorr.mat(p - 1, rho))
    
    #don't know why this doesnt want to work
#     x.train.loop1 <- mvrnorm(n = 1, mu = rep(0, p - 1), 
#                 Sigma = autocorr.mat(p - 1, rho), empirical = TRUE)
    
    x.train.m <- rbind(x.train.m, x.train.loop)
    
  }
    
    x.train <- cbind(rep(1, n), x.train.m)
    
    
#          x.train <- mvrnorm(n = n, mu = rep(0, p), 
#                        Sigma = autocorr.mat(p, rho), empirical = TRUE)
#     

#      x.train <- rmvnorm(n = n, mean = rep(0, p), 
#                          sigma = autocorr.mat(p, rho))

    
    
    #then find the density of y (a vector)
    #based on the x matrix but I need to specify
    #my true betas here and my noise here too
    
    #dmvnorm knows to multiply x.train and my beta.truth together
    #dmvnorm is using chol algorithm to estimate the inverse of the 
    #covariance
    #matrix so it has to be pxp matrix (not the usual NxN matrix)
    
    
    #this y works too

#it's faster to say y.train = X*B + noise * rnorm(x)

#      y.train <- mvrnorm(n = 1, mu = x.train %*% beta.truth, 
#          Sigma = diag(noise, nrow = n))
    
#this is the y that I want
#    y.train <- t(rmvnorm(n = 1, mean = x.train %*% beta.truth, 
#         sigma = diag(noise, nrow = n)))
# 

y.train <- x.train %*% beta.truth + noise * rnorm(1, mean = 0, sd = 1)

# x.test1 <- cbind(rep(1, p), x.test)
# beta.truth1 <- rep(1, p + 1)
# 
# y.train1 <- t(rmvnorm(n = 1, mean = x.train1 %*% beta.truth1, 
#                      sigma = diag(noise, nrow = n)))
# 
#     
#### Model Estimation ####
# estimate models

# Ordinary Least Squares. No intercept
# 
#   fitOLS = lm(y.train ~ 0 + x.train) 

#will fit the intercept

fitOLS = lm(y.train ~  x.train) 
# 
#  fitOLS <- lm(y.train ~  0 + x.train)

#  fitOLS  = glmnet(x.train, y.train, alpha = 0, lambda = 0,
#                    intercept = FALSE)  #this is it!
# # # 
# olsbeta = coef(fitOLS2)[-1]  #this is it!



# glmnet automatically standardizes the predictors
# Ridge Regression. No intercept

fitRidge = glmnet(x.train, y.train, alpha = 0, intercept = TRUE)

# The Lasso. No intercept

fitLasso = glmnet(x.train, y.train, alpha = 1, intercept = TRUE) 

    
    

#### Model Selection ####

# (10-fold) cross validation for the Lasso

cvLasso = cv.glmnet(x.train, y.train, alpha = 1)


# (10-fold) cross validation for Ridge Regression
cvRidge = cv.glmnet(x.train, y.train, alpha = 0)


### Extract Coefficients ###
# OLS coefficient estimates
# 
# betaHatOLS1 = fitOLS1$coefficients

# betaHatOLS = fitOLS$coefficients

# betaHatOLS = coef(fitOLS)[-1] 

betaHatOLS = coef(fitOLS)

# Lasso coefficient estimates 
# s is lambda
#remove the first row because it is the intercept row which is blank

betaHatLasso = as.double(coef(fitLasso, s = cvLasso$lambda.1se))[-1]  

# Ridge  coefficient estimates 
#remove the first row because it is the intercept row which is blank

betaHatRidge = as.double(coef(fitRidge, s = cvRidge$lambda.1se))[-1]

    
#calculate the MSE of OLS, Lasso, and Ridge using the true betas


# MSEOLS1 <- mean((betaHatOLS1 - beta.truth1) ^ 2)

MSEOLS <- sqrt(mean(abs(betaHatOLS - beta.truth)))


MSELasso <- sqrt(mean(abs(betaHatLasso - beta.truth)))
MSERidge <- sqrt(mean(abs(betaHatRidge  - beta.truth)))

dat.ls.mse <- rbind(dat.ls.mse, MSEOLS)
# 
# dat.ls.mse1 <- rbind(dat.ls.mse1, MSEOLS1)

dat.lasso.mse <- rbind(dat.lasso.mse, MSELasso)

dat.ridge.mse <- rbind(dat.ridge.mse, MSERidge)

#calculate the prediction errors of OLS, LAsso and Ridge using a new
#test set for x and a new test set for y


#generate a new x test set using rmvnorm to ensure I get different
#x values

# 
#  x.test <- rmvnorm(n = n, mean = rep(0, p), 
#                     sigma = autocorr.mat(p, rho))

# 
 x.test <- mvrnorm(n = n, mu = rep(0, p), 
                        Sigma = autocorr.mat(p, rho), empirical = TRUE)
                  
#generate a new y from the x test set

# y.test <- t(rmvnorm(n = 1, mean = x.test %*% beta.truth, 
#                      sigma = diag(noise, nrow = n)))

y.test <- x.test %*% beta.truth + noise * rnorm(1, mean = 0, sd = 1)

# calculate predicted values

#predict function gives inaccurate values for OLS
# predOLS =  predict(fitOLS, newdata = as.data.frame(x.test))
# 
# predOLS  <- x.test %*% betaHatOLS

# predOLS  <- x.test %*% betaHatOLS

predOLS = predict(fitOLS, x.test)  #this is it!

# predOLS1  <- x.test1 %*% betaHatOLS1
predLasso <- predict(fitLasso, s = cvLasso$lambda.1se, newx = x.test)
predRidge <- predict(fitRidge, s = cvRidge$lambda.1se, newx = x.test)
    
    

# calculate test set prediction errors
PEOLS = sqrt(mean(abs(predOLS - y.test)))

# 
# PEOLS1 = mean((predOLS1 - y.test) ^ 2)

PELasso = sqrt(mean(abs(predLasso - y.test)))
PERidge = sqrt(mean(abs(predRidge - y.test)))


dat.ls.pe <- rbind(dat.ls.pe, PEOLS)

# dat.ls.pe1 <- rbind(dat.ls.pe1, PEOLS1)

dat.lasso.pe <- rbind(dat.lasso.pe, PELasso)

dat.ridge.pe <- rbind(dat.ridge.pe, PERidge)
    

#standardize y.test to do prediction errors
#This makes the prediction errors much worse than not 
#standardizing y

# y.test.stand <- scale(y.test, center = TRUE, scale = TRUE)
# 
# PELasso.stan = mean((predLasso - y.test.stand) ^ 2)
# PERidge.stan = mean((predRidge - y.test.stand) ^ 2)

print(k)

    
  }
  
  #find the mean MSE for OLS, Ridge, Lasso
  
  mean.MSE.OLS <- mean(dat.ls.mse)

# mean.MSE.OLS1 <- mean(dat.ls.mse1)
  
  mean.MSE.Ridge <- mean(dat.ridge.mse)

 mean.MSE.Lasso <- mean(dat.lasso.mse)

#find the mean prediction error for OLS, Ridge, Lasso

mean.PE.OLS <- mean(dat.ls.pe)

# mean.PE.OLS1 <- mean(dat.ls.pe1)

mean.PE.Ridge <- mean(dat.ridge.pe)

mean.PE.Lasso <- mean(dat.lasso.pe)
  
  #find the standard error of the mean MSE for OLS, Ridge, Lasso
  
  se.mean.MSE.OLS <- sqrt(var(dat.ls.mse)/S)
# 
# se.mean.MSE.OLS1 <- sqrt(var(dat.ls.mse1)/S)
  
  se.mean.MSE.Ridge <- sqrt(var(dat.ridge.mse)/S)

se.mean.MSE.Lasso <- sqrt(var(dat.lasso.mse)/S)

#find the standard error of the mean prediction error for OLS, Ridge, Lasso

se.mean.PE.OLS <- sqrt(var(dat.ls.pe)/S)
# 
# se.mean.PE.OLS1 <- sqrt(var(dat.ls.pe1)/S)

se.mean.PE.Ridge <- sqrt(var(dat.ridge.pe)/S)

se.mean.PE.Lasso <- sqrt(var(dat.lasso.pe)/S)

  
  #output the means and standard errors because this is my function
   list(MeanMSEOLS = mean.MSE.OLS, 
       MeanMSERidge = mean.MSE.Ridge, 
       MeanMSELasso = mean.MSE.Lasso, meanPEOLS = mean.PE.OLS,
        
       meanPERidge = mean.PE.Ridge, meanPELasso = mean.PE.Lasso,
       seMSEOLS = se.mean.MSE.OLS, 
       seMSERidge = se.mean.MSE.Ridge,
       seMSELasso = se.mean.MSE.Lasso, sePEOLS = se.mean.PE.OLS,
    
       sePERidge = se.mean.PE.Ridge, sePELasso = se.mean.PE.Lasso)
  
  
}


#Now find it for all treatment combos


biglist.mean.MSEOLS <- list()
# 
# biglist.mean.MSEOLS1 <- list()

biglist.mean.MSERidge <- list()

biglist.mean.MSELasso <- list()

biglist.mean.PEOLS <- list()

# 
# biglist.mean.PEOLS1 <- list()

biglist.mean.PERidge <- list()

biglist.mean.PELasso <- list()

biglist.se.MSEOLS <- list()
# 
# biglist.se.MSEOLS1 <- list()

biglist.se.MSERidge <- list()

biglist.se.MSELasso <- list()

biglist.se.PEOLS <- list()

# 
# biglist.se.PEOLS1 <- list()

biglist.se.PERidge <- list()

biglist.se.PELasso <- list()


#Code to make sure no combo has the same seed number.  There are no 
#duplicates
#  
#  for (i in names(n)){
#         for (j in names(rho)){
#               print(((n[[i]] - 1) * 3 + (rho[[j]] - 1) * 11))
#           }
#  }
#  
#   for (i in names(n)){
# #     
#                print(n[[i]])
# # 
#   }
# # 
# # 
#    for (j in names(rho)){
# # 
#       print(rho[[j]])
#      }
# # 

# Start the clock!
ptm2 <- proc.time()

for (i in names(n)){
  for (j in names(rho)){
    #   print(typeof((I[[i]] - 1) * 3 + (d[[j]] - 1) * 11))
    set.seed((n[[i]] - 1) * 3 + (rho[[j]] - 1) * 11)
    
    one.combo <- generate.shrink(S, n[[i]], rho[[j]])
    
    MeanMSEOLS <- list(one.combo$MeanMSEOLS)
    
    MeanMSERidge <- list(one.combo$MeanMSERidge)
        
    MeanMSELasso <- list(one.combo$MeanMSELasso)
    
    MeanPELasso <- list(one.combo$meanPELasso)
    
    MeanPERidge <- list(one.combo$meanPERidge)
    
    MeanPEOLS <- list(one.combo$meanPEOLS)
    
    seMSEOLS <- list(one.combo$seMSEOLS)
    
    seMSERidge <- list(one.combo$seMSERidge)
    
    seMSELasso <- list(one.combo$seMSELasso)
        
    sePEOLS <- list(one.combo$sePEOLS)
    
    sePERidge <- list(one.combo$sePERidge)
    
    sePELasso <- list(one.combo$sePELasso)
    
        
   name <- paste('item:', (n[[i]] - 1)*3 + (rho[[j]] - 1)*11, sep = '')
    
    
    biglist.mean.MSEOLS[[name]] <- MeanMSEOLS 
    
    biglist.mean.MSERidge[[name]] <- MeanMSERidge
    
    biglist.mean.MSELasso[[name]] <- MeanMSELasso
    
    biglist.mean.PEOLS[[name]] <- MeanPEOLS
    
    biglist.mean.PERidge[[name]] <- MeanPERidge
    
    biglist.mean.PELasso[[name]] <-  MeanPELasso
    
    biglist.se.MSEOLS[[name]] <- seMSEOLS 
    
    biglist.se.MSERidge[[name]] <- seMSERidge
    
    biglist.se.MSELasso[[name]] <- seMSELasso
    
    biglist.se.PEOLS[[name]] <-  sePEOLS
    
    biglist.se.PERidge[[name]] <- sePERidge
    
    biglist.se.PELasso[[name]] <- sePELasso 
    
    
    
    
    #output
    biglist.mean.MSEOLS 
    
    biglist.mean.MSERidge 
    
    biglist.mean.MSELasso 
    
    biglist.mean.PEOLS
    
    biglist.mean.PERidge 
    
    biglist.mean.PELasso
    
    biglist.se.MSEOLS
    
    biglist.se.MSERidge
    
    biglist.se.MSELasso 
    
    biglist.se.PEOLS 
    
    biglist.se.PERidge 
    
    biglist.se.PELasso 
    
    
    
    print(j)
    
  }
  print(i)
  
}  

# Stop the clock
proc.time() - ptm2

warnings()

totalfile <- save.image(file = "Allproject5.rdata")


#Try using kable function to make a nice table from all of these lists
#Dr. Zhou has an example in HW

#give information on the functions I wrote and clean up my data













