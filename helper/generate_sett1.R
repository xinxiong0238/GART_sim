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



gen_coeff_eqSource <- function(sd_source, p=30, s0, s1, L){
  target_part = c(0.3, 0.1, 0.5, -0.2, -0.3, 0.2, 0.1, 0.15, 0.03, 0.01, -0.15, -0.03, -0.01)
  target_beta = c(target_part, rep(0.1, 17))
  other_null = sapply(1:(L-1), function(i){
    target_part + rnorm(13, mean = 0, sd = sd_source)
  })
  other_notnull = sapply(1:(L-1), function(i){
    # junk = rnorm(17, mean = 0.1, sd = sd_source)
    if(i %% 2 == 0){
      sign_i = 1
    }else{
      sign_i = -1
    }
    junk = rep(0.1, 17) + sign_i * sd_source * ceiling(i/2)/(L-1)
    return(junk)
  })
  other = rbind(other_null, other_notnull)
  Beta_matrix = cbind(
    # target_beta + rnorm(30, sd = 0.05),
    other,
    -target_beta + rnorm(30, sd = 0.05))
  return(list(Beta_matrix, target_beta, gamma_0 = rep(1/L, L)))
}



##### Simulation
simu_eqSource <- function(sd_source, L, ratio_target_n, s_Q=1, s_P=1){
  coeff_matrix = gen_coeff_eqSource(sd_source, p, s0, s1, L)
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
  
  Sigma_Q = X_cov_target(p_all)
  # set.seed(N_target)
  X_target = MASS::mvrnorm(n = N_target, mu = rep(0, p_all),
                           Sigma = Sigma_Q)
  # set.seed(seed)
  epsilon = matrix(rnorm(N_target, mean = 0, sd = s_Q), N_target, 1)
  
  ## beta contains noise..?
  beta_target = coeff_matrix[[2]]
  ## target =d source
  y_target = X_target %*% beta_target + epsilon
  
  N_valid = 5000
  X_target_valid = MASS::mvrnorm(n = N_valid, mu = rep(0, p_all),
                                 Sigma = Sigma_Q)
  epsilon_valid = matrix(rnorm(N_valid, mean = 0, sd = s_Q), N_valid, 1)
  
  beta_valid_target_weight = 1
  par_multinom = MCMCpack::rdirichlet(1, alpha = rep(1, L))
  par_multinom = c(beta_valid_target_weight,
                   par_multinom * (1 - beta_valid_target_weight))
  beta_valid = cbind(coeff_matrix[[2]],
                     coeff_matrix[[1]]) %*% par_multinom
 
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

