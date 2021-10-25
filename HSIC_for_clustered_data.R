#### HSIC-based independence test for clustered data (using Davies approximation for null distribution) ####
library(CompQuadForm)

## Input arguments:

# m: number of clusters
# expo.K.c: centered exposure kernel matrix (with entries ordered by cluster)
# outcome.K.c: centered outcome kernel matrix (with entries ordered by cluster)
# kernel_type: kernel type (optional). If kernel_type = 'linear', more efficient computation can be used.
# n_expo: the number of exposure variables (specify when kernel_type = 'linear')
# n_outcome: the number of outcome variables (specify when kernel_type = 'linear')
# block_size: parameter used to speed up computation/save memory (optional). 
#             Only used when either expo.K.c or outcome.K.c has n non-zero eigenvalues.
#             Must be a factor of sample size n.

HSIC_cl <- function(m, expo.K.c, outcome.K.c, kernel_type = NULL, n_expo = NULL,
                     n_outcome = NULL, block_size = NULL) {
  # cluster size
  d <- nrow(expo.K.c)/m
  
  # eigendecomposition of kernel matrices
  outcome.K.decomp <- chol.de(outcome.K.c)
  expo.K.decomp <- chol.de(expo.K.c)
  
  if (is.null(kernel_type) | kernel_type != 'linear') {
    R1 <- outcome.K.decomp$R
    R2 <- expo.K.decomp$R 
  } else if (kernel_type == 'linear') {
    R1 <- as.matrix(outcome.K.decomp$R[, 1:min(n_outcome, m*d)] )
    R2 <- as.matrix(expo.K.decomp$R[, 1:min(n_expo, m*d)] )
  }
  
  if (ncol(R1) == m*d | ncol(R2) == m*d) {  # Divide the R1 or R2 matrix into blocks to save memory
    
    if (is.null(block_size))  {
      block_size <- m
    }
    n_block <- m*d/block_size
    
    Sigma_hat <- matrix(0, m, m)
    
    if (ncol(R1) == m*d ) {
      for (s in 1:n_block) {
       V_tilde <- lapply(1:m, function(i){
          get_V_tilde(i, R1[,(block_size*s-block_size+1):(block_size*s)], R2, d)})
       V_tilde <- do.call(cbind, V_tilde)
       Sigma_hat <- Sigma_hat + crossprod(V_tilde, V_tilde)
      }
    } else {
        for (s in 1:n_block) {
          V_tilde <- lapply(1:m, function(i){
            get_V_tilde(i, R1, R2[,(block_size*s-block_size+1):(block_size*s)], d)})
          V_tilde <- do.call(cbind, V_tilde)
          Sigma_hat <- Sigma_hat + crossprod(V_tilde, V_tilde)
        }
    }
    
    rm(V_tilde)
    
  } else {
    
    V_tilde <- lapply(1:m, function(i){get_V_tilde(i, R1, R2, d)})
    V_tilde <- do.call(cbind, V_tilde)
    Sigma_hat <- crossprod(V_tilde, V_tilde)
    
    rm(V_tilde)
    
  }

  
  Sigma_hat.eiv <- eigen(Sigma_hat, symmetric=T, only.values=T)$values  
  Sigma_hat.eiv <- Sigma_hat.eiv[Sigma_hat.eiv > 1e-08]
  
  rm(Sigma_hat)
  
  test_stat <- sum(outcome.K.c * expo.K.c)
  
  eiv.sort <- sort(Sigma_hat.eiv, decreasing=T)
  
  ## derive p-value of HSIC statistic:
  results_score <- davies(test_stat, eiv.sort)
  davies_pvalue <- results_score$Qq
  
  return(list('test_stat'=test_stat,  'pval'=davies_pvalue))
} 



#### Helper functions ####

## eigendecomposition of a psd matrix as X = RR^T
chol.de <- function(x) {
  ei <- eigen(x, symmetric = T)
  eiv.init <- ei$values
  eiv.noneg <- (eiv.init+abs(eiv.init))/2
  R <- ei$vectors %*% diag(sqrt(eiv.noneg))
  return(list("R"=R, "ei.value" = eiv.noneg))
}

## helper function to estimate the null distribution of HSIC statistic
get_V_tilde <- function(i, R1, R2, d) {
  R1_psi <- R1[(d*i-d+1):(d*i),]
  R2_psi <- R2[(d*i-d+1):(d*i),]
  R1_R2_psi_prod <- 0
  for (j in 1:d) {
    R1_R2_psi_prod <- R1_R2_psi_prod + as.vector(R1_psi[j,] %*% t(R2_psi[j,]))
  }
  return(R1_R2_psi_prod)
} 
