#!/usr/local/bin/Rscript

#### Investigate Type I error and Power of HSIC test for clustered data ####

## Load libraries
library(argparse)
library(MASS)
library(CompQuadForm)
source("HSIC_for_clustered_data.R")
source("useful_functions.R")




## load command line arguments
parser <- ArgumentParser()

parser$add_argument("--power", type = "double",
                    help = "indicator for power (1,2) or type I error (0) analysis")
parser$add_argument("--m", type = "double",
                    help = "number of clusters")
parser$add_argument("--d", type = "double",
                    help = "cluster size")
parser$add_argument("--n_expo", type = "double",
                    help = "number of exposures/outcomes")
parser$add_argument("--rho_c", type = "double",
                    help = "parameter for the within-cluster correlation structure")
parser$add_argument("--gamma", type = "double",
                    help = "proportion of outcomes associated with the exposure")
parser$add_argument("--kernel", type = "character",
                    help = "type of kernel (gaussian or linear)")
parser$add_argument("--n_sim", type = "double",
                    help = "number of data sets to be simulated")
parser$add_argument("--seed_num", type = "double",
                    help = "seed number for simulation")
parser$add_argument("--output_file", type = "character",
                    help = "output file name")

args <- parser$parse_args()


## use command line arguments
power_ind <- args$power  # indicator for power (1,2) or type I error (0) analysis
m <- args$m # number of clusters 
d <- args$d # cluster size
n_expo <- args$n_expo # number of exposures/outcomes 
rho_c <- args$rho_c # parameter for the within-cluster correlation structure
gamma <- args$gamma # proportion of outcomes associated with the exposure
kernel_type <- args$kernel  # type of kernel (linear or Gaussian)
n_sim <- args$n_sim # number of simulations
seed_num = args$seed_num # seed for simulation





#### Baseline parameters ####
## within-cluster covariance matrix, with AR(1) structure
cluster_cov_mt <- ar1_cor(d, rho_c)

## covariance matrix for exposure variables
expo_cov_mt <- matrix(0.5, n_expo, n_expo)  
diag(expo_cov_mt) <- 1

## covariance matrix for outcome variables
outcome_cov_mt <- matrix(0.5, n_expo, n_expo)  
diag(outcome_cov_mt) <- 1

## centering matrix (observation level)
I.m=diag(1,m*d)
I.1=rep(1,m*d)
H=I.m-I.1%*%t(I.1)/(m*d)

## centering matrix (cluster level)
I.m.indep=diag(1,m)
I.1.indep=rep(1,m)
H.indep=I.m.indep-I.1.indep%*%t(I.1.indep)/m


#### main simulation codes ####
sim_func <- function(i) {
  ## simulate multivariate normal exposure variables
  expo_list <- replicate(m, sim_mvn_in_cluster(5, n_expo, expo_cov_mt, cluster_cov_mt, d),
                         simplify = F)
  expo_matrix <- do.call(rbind, expo_list)
  
  if (kernel_type == 'gaussian')  {
    expo.K <- kern_g(expo_matrix)    # gaussian kernel
  }  else if (kernel_type ==  'linear') {
    expo.K <- expo_matrix %*% t(expo_matrix)  # linear kernel
  }
  
  ## simulate multivariate normal outcome variables
  outcome_list <- replicate(m, sim_mvn_in_cluster(0, n_expo, outcome_cov_mt, cluster_cov_mt, d),
                          simplify = F)
  outcome_matrix_init <- do.call(rbind, outcome_list)
  
  if (power_ind == 1) {  # Power Scenario 1: linear effect of exposure on outcome
    
    ## randomly select one exposure variable to be the causal exposure
    expo_idx <- sample(1:n_expo, 1)
    aff_expo_vec <- expo_matrix[,expo_idx]
    
    # effect size
    beta <- rep(0, n_expo)
    n_outcome_aff <- gamma*n_expo  # number of outcomes affected
    beta[1:n_outcome_aff] <- runif(n = n_outcome_aff, min = 0, max = sqrt(25/m))
    expo_effect <- lapply(beta, function(x){x*aff_expo_vec})
    expo_effect_mt <- do.call(cbind, expo_effect)
    
    # add the exposure effect to the outcomes
    outcome_matrix <- outcome_matrix_init + expo_effect_mt
    
  }  else if (power_ind == 2) {   # Power Scenario 2: nonlinear effect of exposure on outcome
    ## randomly select one exposure variable to be the causal exposure
    expo_idx <- sample(1:n_expo, 1)
    aff_expo_vec <- expo_matrix[,expo_idx] - 4
    aff_expo_vec <- log(aff_expo_vec^2)
    
    # effect size
    beta <- rep(0, n_expo)
    n_outcome_aff <- gamma*n_expo  # number of outcomes affected
    beta[1:n_outcome_aff] <- runif(n = n_outcome_aff, min = 0, max = sqrt(10/m))
    expo_effect <- lapply(beta, function(x){x*aff_expo_vec})
    expo_effect_mt <- do.call(cbind, expo_effect)
    
    # add the exposure effect to the outcomes
    outcome_matrix <- outcome_matrix_init + expo_effect_mt
    
  } else if (power_ind == 0) {   # independence between exposure and outcome
    outcome_matrix <- outcome_matrix_init
  }
  
  ## construct a kernel based on the outcome matrix
  if (kernel_type == 'gaussian')  {
    outcome.K <- kern_g(outcome_matrix)  # gaussian kernel
  } else if (kernel_type == 'linear') {
    outcome.K <- outcome_matrix %*% t(outcome_matrix)  # linear kernel
  }
  
  # center the kernel matrices
  outcome.K.c <- H%*% outcome.K %*%H
  expo.K.c <- H%*% expo.K %*%H
  
  
  if (power_ind == 0) {  # assess type I error rate
    # Perform HSIC_orig
    p_value_orig <- HSIC_orig(expo.K.c, outcome.K.c)$pval
    
  }  
  

  indep_idx1 <- seq(from=1, to=m*d-d+1, by=d)
  indep_idx2 <- seq(from=2, to=m*d-d+2, by=d)
  indep_idx3 <- seq(from=3, to=m*d, by=d)
  
  #### Perform HSIC_mean: Use the mean across observations from each cluster ####
  expo_matrix.indep <- (expo_matrix[indep_idx1,] + expo_matrix[indep_idx2,] +
                          expo_matrix[indep_idx3,])/3
  outcome_matrix.indep <- (outcome_matrix[indep_idx1,] + outcome_matrix[indep_idx2,] +
                           outcome_matrix[indep_idx3,])/3
  
  if (kernel_type == 'gaussian') {
    expo.K.indep <- kern_g(expo_matrix.indep)
    outcome.K.indep <- kern_g(outcome_matrix.indep)
  } else if (kernel_type == 'linear') {
    expo.K.indep <- expo_matrix.indep %*% t(expo_matrix.indep)
    outcome.K.indep <- outcome_matrix.indep %*% t(outcome_matrix.indep)
    
  }
  
  # center the kernel matrices
  outcome.K.indep.c <- H.indep %*% outcome.K.indep %*%H.indep
  expo.K.indep.c <- H.indep %*% expo.K.indep %*%H.indep
  
  # Get p-value
  p_value_mean <- HSIC_orig(expo.K.indep.c, outcome.K.indep.c)$pval
  
  
  #### Perform HSIC_cat: Concatenate all observations at cluster level ####
  expo_matrix.indep <- cbind(expo_matrix[indep_idx1,], expo_matrix[indep_idx2,],
                             expo_matrix[indep_idx3,])  
  outcome_matrix.indep <- cbind(outcome_matrix[indep_idx1,], outcome_matrix[indep_idx2,],
                              outcome_matrix[indep_idx3,])   
  
  if (kernel_type == 'gaussian') {
    expo.K.indep <- kern_g(expo_matrix.indep)
    outcome.K.indep <- kern_g(outcome_matrix.indep)
  } else if (kernel_type == 'linear') {
    expo.K.indep <- expo_matrix.indep %*% t(expo_matrix.indep)
    outcome.K.indep <- outcome_matrix.indep %*% t(outcome_matrix.indep)
  }
  
  # center the kernel matrices
  outcome.K.indep.c <- H.indep %*% outcome.K.indep %*%H.indep
  expo.K.indep.c <- H.indep %*% expo.K.indep %*%H.indep
  
  # Get p-value
  p_value_cat <- HSIC_orig(expo.K.indep.c, outcome.K.indep.c)$pval
  

  #### Perform HSIC_cl ####
  p_value_cl <- HSIC_cl(m, expo.K.c, outcome.K.c, kernel_type, n_expo,
                           n_expo, block_size = m)$pval
    
  
  
  if (power_ind  == 0) {
      return(c(p_value_orig, p_value_cl))
    } else {
      return(c(p_value_cl, p_value_cat, p_value_mean))
    }
  
}

if (power_ind  == 0) {
   p_value_df <- data.frame('HSIC_orig' = rep(NA, n_sim),
                              'HSIC_cl' = rep(NA, n_sim))
} else {
   p_value_df <- data.frame('HSIC_cl' = rep(NA, n_sim),
                              'HSIC_cat' = rep(NA, n_sim),
                              'HSIC_mean' = rep(NA, n_sim))
}


set.seed(1000*seed_num)
for (i in 1:n_sim) {
  p_value_df[i,] <- sim_func(i)
  gc()
}


write.csv(p_value_df, file=args$output_file)

warnings()








