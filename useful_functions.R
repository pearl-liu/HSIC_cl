#### Functions used in simulations and demo ####
library(MASS)
library(CompQuadForm)
library(reshape2)
library(ggplot2)

## generate AR(1) correlation structure
ar1_cor <- function(d, rho) {
  exponent <- abs(matrix(1:d - 1, nrow = d, ncol = d, byrow = TRUE) - 
                    (1:d - 1))
  return(rho^exponent)
}

## simulate multivariate normal exposure/outcome variables within a cluster
sim_mvn_in_cluster <- function(mu_val, n_expo, var_cov_mt, cluster_cov_mt, d) {
  Y <- mvrnorm(mu = rep(0, d*n_expo) + mu_val, 
                     Sigma = var_cov_mt %x% cluster_cov_mt )
  Y_mt <- matrix(Y, d, n_expo)   
  # rows are observations within each cluster
  # columns are variables
  return(Y_mt)
}

## construct a gaussian kernel matrix
kern_g <- function(input_mt){
  n <- nrow(input_mt)
  D <- matrix(NA, n, n)  # pairwise distance matrix
  for(i in 1:n){
    for(j in i:n){
      D[i,j] <- sum((input_mt[i,]-input_mt[j,])^2)
      if (i != j) {
        D[j,i] <- D[i,j]
      }
    }}	
  D_vec <- c(D)
  tau <- median(D_vec[D_vec>0])   # use the median distance as the bandwidth
  K <- exp(-D/tau)
  return(K)
}

## Helper function for HSIC_orig (using Davies approximation).
## This function is adapted from https://github.com/epstein-software/GAMuT
TestGAMuT <- function(Yc, lambda_Y, Xc, lambda_X) {
  
  ## test statistic:
  m = nrow(Yc) # number of subjects in study
  GAMuT = (1/m) * sum(sum(t(Yc) * Xc))  
  
  ## populate vector of all pairwise combination of eigenvalues
  ## from the phenotype and genotype similarity matrices:
  Z <- (as.matrix(lambda_Y)) %*% t(as.matrix(lambda_X))
  Zsort <- sort(Z, decreasing=T)
  
  ## derive p-value of GAMuT statistic:
  scoredavies = GAMuT*m^2
  results_score <- davies(scoredavies, Zsort)
  davies_pvalue <- (results_score$Qq)
  
  return(list('test_stat'=scoredavies,  'pval'=davies_pvalue))
} 

## HSIC_orig (using Davies approximation)
HSIC_orig <- function(expo.K.c, outcome.K.c) {
  # Obtain eigenvalues of expo.K.c and outcome.K.c
  expo.K.eiv <- eigen(expo.K.c, symmetric=T, only.values=T)$values  
  expo.K.eiv <- expo.K.eiv[expo.K.eiv > 1e-08]
  
  outcome.K.eiv <- eigen(outcome.K.c, symmetric=T, only.values=T)$values  
  outcome.K.eiv <- outcome.K.eiv[outcome.K.eiv > 1e-08]
  
  # Get test statistic and p-value
  return( TestGAMuT(expo.K.c, expo.K.eiv, outcome.K.c, outcome.K.eiv) )
}


## Output empirical type I error rates or powers based on p-value results
# power_ind: 0 = type I error, 1 = power
# sig_level: significance level
result_eval <- function(power_ind, file_name, sig_level) {
  p_val_list <-  read.csv(file_name)
  if (power_ind == 0) { # evaluate type I error
    hsic_orig_type_I_err <- sum(p_val_list$HSIC_orig < sig_level)/nrow(p_val_list)
    hsic_cl_type_I_err <- sum(p_val_list$HSIC_cl < sig_level)/nrow(p_val_list)
    return(c(hsic_orig_type_I_err, hsic_cl_type_I_err))
  } else {  # evaluate power
    hsic_cl_power <- sum(p_val_list$HSIC_cl < sig_level)/nrow(p_val_list)
    hsic_cat_power <- sum(p_val_list$HSIC_cat < sig_level)/nrow(p_val_list)
    hsic_mean_power <- sum(p_val_list$HSIC_mean < sig_level)/nrow(p_val_list)
    return(c(hsic_cl_power, hsic_cat_power, hsic_mean_power))
  }
}

## Make a plot comparing power between different methods
power_plot <- function(power_result, kernel_type, rho_c,
                       effect_size_input, effect_size_output) {
  power_result_m <- melt(power_result, id.vars = 'Method')
  colnames(power_result_m)[2:3] <- c('Number', 'Power') 
  power_result_m$Number <- factor(power_result_m$Number, 
                                 levels = effect_size_input,
                                 labels =  effect_size_output)
  power_result_m$Method <- factor(power_result_m$Method,
                                  levels = c("HSIC_mean", "HSIC_cat", "HSIC_cl"))
  output_plot <- ggplot(power_result_m, aes(Number, Power, fill = Method)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1)) +
    scale_fill_manual(values = c('green','orange','blue'),
              labels = c(expression(HSIC[mean]), expression(HSIC[cat]), 
                         expression(HSIC[cl])))+
    labs(title = bquote(rho[c]~"="~.(rho_c)~","~.(kernel_type)), 
         x='Proportion of outcomes affected') +
    scale_x_discrete(labels = effect_size_output)
  
  return(output_plot)
}


