#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# If you use this code, cite as
# Y.W. Park and D. Klabjan, Iteratively reweighted least squares algorithms for l1-norm principal component analysis, 2016 IEEE International Conference on Data Mining (ICDM)
# preprint available at http://arxiv.org/abs/1609.02997
#
# Functions included: wl1pca, awl1pca
# wl1pca and awl1pca correspond to Algorithms 1 and 2, respectively, in the reference paper
# Date: Sep 17, 2016
# Version: 1.0.1
# Comments:
# Most of the arguments, return values, and their definitions follow the package "pcaL1" by Jot and Brooks (2015)
# Function: awl1pca
# Arguments 
#    X: data, must be in matrix or table form
#    projDim: number of dimensions to project data into, must be an integer, default is 1.
#    center: whether to center the data using the mean, default is TRUE
#    tolerance: sets the convergence tolerance (of weights) for the algorithm, default is 0.001
#    iterations: sets the number of iterations to run before returning the result, default is 200
#    beta: algorithm parameter to set up bound for weights
#    gamma: algorithm parameter to determine whether to use approximation formula or prcomp function
awl1pca <- function(X, projDim = 1, center = TRUE, projections="l2", tolerance = 0.001, iterations = 200, beta = 0.99, gamma = 0.1){
  if (inherits(X, "data.frame")) {
    X <- as.matrix(X)
  }

  if(center){
    X <- apply(X,2,function(y) y - mean(y))
  }
  time.begin <- proc.time()
  n <- dim(X)[1]
  m <- dim(X)[2]
  weights.prev <- matrix(2,n,1)
  weights <- matrix(1,n,1)
  obj.best <- .Machine$double.xmax
  my.epsilon <- tolerance
  my.epsilon.b <- 0.1
  U.weights <- matrix(1,n,1)
  num.iteration <- 1
  continue.algorithm <- 1
  my.beta <- beta^num.iteration
  
  while(continue.algorithm){
    weights.prev <- weights
    if(num.iteration > 1){
      my.EV.prev <- my.EV
      my.PC.prev <- as.matrix(my.PC)
    }

    if((num.iteration > 1) && (sum(abs(weight.diff)) < gamma * sum(abs(weights)))){
      X.diff <- diag(as.vector(weight.diff)) %*% X
      X.diff2 <- t(X.diff) %*% X.diff
      my.result <- L2PCA_approx(my.EV.prev, my.PC.prev,projDim,X.diff2)
      my.EV <- my.result$eigenvalues
      my.PC <- my.result$eigenvectors
    }else{
      X.weighted <- diag(as.vector(weights)) %*% X
      my.result <- prcomp(X.weighted)
      my.EV <- my.result$sdev[1:projDim]^2
      my.PC <- as.matrix(my.result$rotation[,1:projDim])
    }
    
    Reconstuction.Errors <- X - X %*% my.PC %*% t(my.PC)
    obj.current <- sum(abs(Reconstuction.Errors))
    if(obj.current < obj.best){
      obj.best <- obj.current
      PC.best <- my.PC
    }
    
    max.U <- 0
    for(i in 1:n){
      L2.norm.obs.i <- sum(Reconstuction.Errors[i,]^2)
      L1.norm.obs.i <- sum(abs(Reconstuction.Errors[i,]))
      if(L2.norm.obs.i > my.epsilon.b){
        U.weights[i,1] <- L1.norm.obs.i/L2.norm.obs.i
        max.U <- max(max.U,U.weights[i,1])
      }else{
        U.weights[i,1] <- -1
      }
    }
    U.weights[U.weights[,1]==(-1),1]<-max.U
    weights[U.weights[,1]<weights[,1]*(1-my.beta),1] <- weights[U.weights[,1]<weights[,1]*(1-my.beta),1]*(1-my.beta)
    weights[U.weights[,1]>weights[,1]*(1+my.beta),1] <- weights[U.weights[,1]>weights[,1]*(1+my.beta),1]*(1+my.beta)
    weights[(weights[,1]*(1-my.beta)<=U.weights[,1]) & (U.weights[,1]<=weights[,1]*(1+my.beta)),1] <- U.weights[(weights[,1]*(1-my.beta)<=U.weights[,1]) & (U.weights[,1]<=weights[,1]*(1+my.beta)),1]
    
    my.beta <- beta^num.iteration
    weight.diff <- weights.prev - weights
    if((num.iteration >= iterations)||(sum(abs(weight.diff))<my.epsilon)){
      continue.algorithm <- 0
    }
    
    num.iteration <- num.iteration + 1
    if(num.iteration > 2){
      if(norm(as.matrix(my.PC) - as.matrix(my.PC.prev),"F") < my.epsilon){
        continue.algorithm <- 0
      }
    }
    cat(".", file=stderr())
  }
  cat("\n", file=stderr())
  time.exe <- as.numeric((proc.time() - time.begin)[3])
  
  Reconstuction.Errors <- X - X %*% PC.best %*% t(PC.best)
  Reconstructions <- X %*% PC.best %*% t(PC.best)
  Projected.Points = X %*% PC.best

  if (projections == "l1") {
    myl1projection <- l1projection(X, as.matrix(PC.best))
    Reconstructions <- myl1projection$projPoints
    Projected.Points <- myl1projection$scores
  }

  solution <-  list(loadings = PC.best, scores = Projected.Points, projPoints = Reconstructions, L1error = obj.best, nIter = num.iteration, ElapsedTime = time.exe)
  class(solution) <- "awl1pca"
  
  cat("# Result summary \n", file=stderr())
  cat("L1 error = ", obj.best, "\n", file=stderr())
  cat("Num iterations = ", num.iteration, "\n", file=stderr())
  cat("Elapsed time = ", time.exe, "seconds \n\n", file=stderr())
  solution
}


# Function: L2PCA_approx
# implementation of Equations (11) and (12) in the reference paper.
L2PCA_approx <- function(ev.prev, pc.prev, projDim, X.diff){
  my.ev <- ev.prev
  my.pc <- pc.prev
  for(k in 1:projDim){
    my.ev[k] <- ev.prev[k] + t(pc.prev[,k]) %*% X.diff %*% pc.prev[,k]
    for(k2 in 1:projDim){
      if(k != k2){
        my.pc[,k] <- my.pc[,k] + (((t(pc.prev[,k2]) %*% X.diff %*% pc.prev[,k])/(ev.prev[k] - ev.prev[k2])) %*% pc.prev[,k2])
      }
    }
  }
  appr_result <- list(eigenvalues = my.ev, eigenvectors = my.pc)
  return (appr_result)
}
