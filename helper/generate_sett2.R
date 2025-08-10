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



gen_coeff_bv <- function(type, dim, alpha = 0.01){
  Beta_matrix = c()
  if(dim == 'low'){
    Beta_matrix = matrix(c(0.3, 0.1, 0.5, -0.2, -0.7, 0, 0, 
                           0.2, 0.05, 0.4, -0.15, -0.6,0, 0,
                           0.3, 0.2, 0.6, -0.3, -0.6,0, 0, 
                           -0.2, -0.1, -0.5, 0.2, 0.6,0, 0), ncol = 4)
  }else{
    Beta_matrix =  matrix(c(0.3, 0.1, 0.5, -0.2, -0.7, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0, 1, 0, 0, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0, 1, 0, 0, 
                            0.2, 0.05, 0.4, -0.15, -0.6, 0, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0, -1, 0,0, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0, -1, 0,
                            0.3, 0.2, 0.6, -0.3, -0.6, 0, 0, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0, 1,0, 0, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0, 1,
                            -0.2, -0.1, -0.5, 0.2, 0.6, 0, 0, 0, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0, 0, 0, 0, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0), ncol = 4)
  }
  if(type == 'p1'){
    target_beta = rowMeans(Beta_matrix[,-1]) + alpha * rep(c(-1, 1, -1, -1, 1), 7)
    gamma_0 = c(1/3, 1/3, 1/3, 0)
  }else{
    target_beta = c(0.2, 0.5, -0.3, -0.2, 0.4)
    gamma_0 = c(0, 0, 0, 0)
  }
  return(list(Beta_matrix, target_beta, gamma_0))
}




##### Simulation
simu_bv <- function(alpha, ratio_target_n, s_Q=1, s_P=1,
                    valid_type = 0, beta_valid_diff = 0, beta_valid_target_weight = 0){
  coeff_matrix = gen_coeff_bv(type = 'p1', dim = 'high', alpha)
  X_lst = list()
  y_lst = list()
  N_vec = c(20000, 20000, 20000, 20000)
  L = length(N_vec)
  N_target = round(mean(N_vec)*ratio_target_n)
  for(l in 1:4){
    N = N_vec[l]
    beta = coeff_matrix[[1]][,l]
    p = length(beta)
    X_l = MASS::mvrnorm(n = N, mu = rnorm(p, mean = rexp(1), sd = 0.01),
                        Sigma = X_cov_source(p))
    epsilon = matrix(rnorm(N, mean = 0, sd = s_P), N, 1)
    y_l = X_l %*% beta + epsilon
    X_lst = c(X_lst, list(X_l))
    y_lst = c(y_lst, list(y_l)) 
  }
  Sigma_Q = X_cov_target(p)
  X_target = MASS::mvrnorm(n = N_target, mu = rep(0, p),
                           Sigma = Sigma_Q)
  epsilon = matrix(rnorm(N_target, mean = 0, sd = s_Q), N_target, 1)
  
  ## beta contains noise..?
  beta_target = coeff_matrix[[2]]
  ## target =d source
  y_target = X_target %*% beta_target + epsilon
  
  N_valid = 5000
  X_target_valid = MASS::mvrnorm(n = N_valid, mu = rep(0, p),
                                 Sigma = Sigma_Q)
  epsilon_valid = matrix(rnorm(N_valid, mean = 0, sd = s_Q), N_valid, 1)
  if(valid_type == 0){
    beta_valid = coeff_matrix[[2]]
    y_target_valid = X_target_valid %*% beta_valid + epsilon_valid
  }
  if(valid_type == 1){
    # beta_noise_ma = beta_target + matrix(rnorm(p*10, sd = 0.2), ncol = 10)
    gamma_source = Variable(L)
    model = Problem(Minimize(quad_form(beta_target - coeff_matrix[[1]] %*% gamma_source, Sigma_Q)),
                    constraints = list(gamma_source >= 0, sum(gamma_source) == 1))
    result <- solve(model)
    gamma_valid = result$getValue(gamma_source)
    beta_noise_ma = coeff_matrix[[1]] %*% gamma_valid
    
    beta_valid = apply(beta_noise_ma, 2, function(beta_noise){
      beta = Variable(p)
      ben_dist = t(beta_noise - beta_target) %*% Sigma_Q %*% (beta_noise - beta_target)
      model = Problem(Minimize(quad_form(beta - beta_noise, Sigma_Q)),
                      constraints = list(quad_form(beta - beta_target, Sigma_Q) <= ben_dist*beta_valid_target_weight))
      result <- solve(model)
      beta_valid = result$getValue(beta)
      return(beta_valid)
    })
    y_target_valid = apply(beta_valid, 2, function(beta_noise){
      epsilon_valid = matrix(rnorm(N_valid, mean = 0, sd = s_Q), N_valid, 1)
      X_target_valid %*% beta_noise + epsilon_valid
    })
  }
  
  if(valid_type == 2){
    beta_noise_ma = c()
    for(n_valid_beta in 1){
      B_0 = cbind(coeff_matrix[[2]], coeff_matrix[[1]])
      gamma_source = Variable(L)
      model = Problem(Minimize(quad_form(beta_target - coeff_matrix[[1]] %*% gamma_source, Sigma_Q)),
                      constraints = list(gamma_source >= 0, sum(gamma_source) == 1))
      result <- solve(model)
      ben_dist = result$value
      judge = T
      while(judge>=1){
        par_multinom = as.vector(MCMCpack::rdirichlet(1, alpha = c(5, rep(1, L))))
        beta_valid = B_0 %*% par_multinom
        judge = t(beta_valid - beta_target) %*% Sigma_Q %*% (beta_valid - beta_target) <= 
          max(c(beta_valid_target_weight*ben_dist, 0.0001))
        judge = 1 - as.numeric(judge)
      }
      # print(par_multinom)
      beta_noise_ma = cbind(beta_noise_ma, beta_valid)
    }
    beta_valid = beta_noise_ma
    y_target_valid = apply(beta_valid, 2, function(beta_noise){
      epsilon_valid = matrix(rnorm(N_valid, mean = 0, sd = s_Q), N_valid, 1)
      X_target_valid %*% beta_noise + epsilon_valid
    })
  }
  
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





gen_coeff_eqSource <- function(sd_source, p, s0, s1, L){
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
simu_bv_new <- function(sd_source, p, s0, s1, L, ratio_target_n, s_Q=1, s_P=1, beta_valid_target_weight = 0){
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
  X_target = MASS::mvrnorm(n = N_target, mu = rep(0, p_all),
                           Sigma = X_cov_target(p_all))
  epsilon = matrix(rnorm(N_target, mean = 0, sd = s_Q), N_target, 1)
  
  ## beta contains noise..?
  beta_target = coeff_matrix[[2]]
  ## target =d source
  y_target = X_target %*% beta_target + epsilon
  
  N_valid = 5000
  X_target_valid = MASS::mvrnorm(n = N_valid, mu = rep(0, p_all),
                                 Sigma = X_cov_target(p_all))
  epsilon_valid = matrix(rnorm(N_valid, mean = 0, sd = s_Q), N_valid, 1)
  
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


