X_cov_source <- function(p){
  cov_ma = sapply(1:p, function(i){
    0.6^abs(i-(1:p))
  })
  diag(cov_ma) = diag(cov_ma) + rnorm(p, mean = 0.1, sd = 0.01)
  # cov_ma = diag(x = rep(1, p))
  return(cov_ma)
}

X_cov_target <- function(p){
  cov_ma = sapply(1:p, function(i){
    0.7^abs(i-(1:p))
  })
  # cov_ma = diag(x = rep(1, p))
  return(cov_ma)
}



gen_coeff_highD_sparseT <- function(p, s0, s1, L){
  Beta_matrix = matrix(rep(c(0.3, 0.1, 0.5, -0.2, -0.7, 0, 0, rep(0.15, 23), rep(-0.15, 20)), L), ncol = L)
  other_null = matrix(0, ncol = L, nrow = p - s0)
  other_notnull = sapply(1:L, function(i){
    s1 * rep(1/s0, s0) * 
      sample(c(1, -1), size = s0, replace = TRUE, prob = c(.5, .5))
  })
  other = rbind(other_notnull, other_null)
  Beta_matrix = rbind(Beta_matrix, other)
  target_beta = c(0.3, 0.1, 0.5, -0.2, -0.7, 0, 0, rep(0.15, 23), rep(-0.15, 20), rep(0, p)) 
  return(list(Beta_matrix, target_beta, gamma_0 = rep(1/L, L)))
}




##### Simulation
simu_highD_sparseT <- function(p, s0, s1, L, ratio_target_n, s_Q=1, s_P=1){
  coeff_matrix = gen_coeff_highD_sparseT(p, s0, s1, L)
  X_lst = list()
  y_lst = list()
  N_vec = rep(20000, L)
  N_target = round(mean(N_vec)*ratio_target_n)
  for(l in 1:L){
    N = N_vec[l]
    beta = coeff_matrix[[1]][,l]
    p_all = length(beta)
    X_l = MASS::mvrnorm(n = N, mu = rnorm(p_all, mean = rexp(1), sd = 0.01),
                        Sigma = X_cov_source(p_all))
    epsilon = matrix(rnorm(N, mean = 0, sd = s_P), N, 1)
    y_l = X_l %*% beta + epsilon
    X_lst = c(X_lst, list(X_l))
    y_lst = c(y_lst, list(y_l)) 
  }
  X_target = MASS::mvrnorm(n = N_target, mu = rep(0, p_all),
                           Sigma = X_cov_target(p_all))
  epsilon = matrix(rnorm(N_target, mean = 0, sd = s_Q), N_target, 1)
  
  beta_target = coeff_matrix[[2]]
  y_target = X_target %*% beta_target + epsilon
  
  N_valid = 5000
  X_target_valid = MASS::mvrnorm(n = N_valid, mu = rep(0, p_all),
                                 Sigma = X_cov_target(p_all))
  epsilon_valid = matrix(rnorm(N_valid, mean = 0, sd = s_Q), N_valid, 1)
  beta_valid = coeff_matrix[[2]]
  y_target_valid = X_target_valid %*% beta_valid + epsilon_valid
  
  target_data = list(`X` = X_target, `y` = y_target, `beta` = beta_target)
  source_data = list(`X` = X_lst, `y` = y_lst, `beta` = coeff_matrix[[1]])
  valid_data = list(`X` = X_target_valid, `y` = y_target_valid, `beta` = beta_valid)
  return(list(`target` = target_data, `source` = source_data,
              `valid` = valid_data,
              `gamma_0` = coeff_matrix[[3]],
              `beta_target` = beta_target,
              `beta_valid` = beta_valid,
              `beta_source` = coeff_matrix[[1]],
              `s_Q` = s_Q, `s_P` = s_P,
              `N_source` = N_vec, `N_target` = length(y_target)))
}

