# Equality constraints
eval_g_eq <- function(x){
  return (sum(x) - 1)
}

GART_function_more_ini <- function(data, tau){
  X_target = data$target$X; y_target = data$target$y
  X_source = data$source$X; y_source = data$source$y
  X_target_valid = data$valid$X
  y_target_valid = data$valid$y
  N_source = data$N_source; N_target = data$N_target
  L = length(X_source)
  p = ncol(X_target)
  ### estimate b(l), Sigma_Q, Gamma_beta
  # message("estimate b(l), Sigma_Q, Gamma_beta")
  X_source_matrix = Reduce(rbind, X_source)
  y_source_vec = unlist(y_source)
  B_source_hat = sapply(1:L, function(l){
    est = cv.glmnet(X_source[[l]], y_source[[l]])
    b_l = est$glmnet.fit$beta
    b_l = b_l[, which(est$lambda == est$lambda.min)]
    return(b_l)
  })
  
  est = cv.glmnet(X_target, y_target)
  b_l = est$glmnet.fit$beta
  beta_target_only = b_l[, which(est$lambda == est$lambda.min)]
  
  base_names = c('weighted', 'mixST', 'mixT0', 'source', 'target', 'zero', 'convex')
  beta_est = gamma_est = c()
  ini_bound = 1
  for(which_beta_to_use in c('target_only')){
    base_est = maximin_base_est = beta_target_hat_lst = c()
    for(fold in c(0, 1)){
      # if(fold == 0){
      #   message('only step 1, first half')
      # }else{
      #   message('only step 1, second half')
      # }
      test_rid = (1:N_target)[((1:N_target) %% 2) == fold]
      train_rid = setdiff(1:N_target, test_rid)
      est = cv.glmnet(X_target[train_rid, ], y_target[train_rid])
      b_l = est$glmnet.fit$beta
      beta_target_hat = b_l[, which(est$lambda == est$lambda.min)]
      
      beta_target_hat_lst = c(beta_target_hat_lst, list(beta_target_hat))
      
      gammaHat <- Variable(ncol(B_source_hat)+1)
      B_hat_1 = cbind(beta_target_hat, B_source_hat)
      objective <- Minimize(mean((X_target[test_rid, ] %*% B_hat_1 %*% gammaHat - 
                                    y_target[test_rid])^2))
      maximin_step1 = step1 = c()
      for(base in base_names){
        # message(paste0('Baseline: ', base))
        problem = switch(base, 
                         convex = Problem(objective, constraints =
                                            list(gammaHat >= 0, gammaHat <= 1, sum(gammaHat)==1)),
                         source = Problem(objective, constraints =
                                            list(gammaHat >= 0, gammaHat <= 1, sum(gammaHat)==1, gammaHat[1] == 0)),
                         mixST = Problem(objective, constraints =
                                           list(gammaHat >= 0, gammaHat <= 1, sum(gammaHat)==1, gammaHat[1] == 0)),
                         mixT0 = Problem(objective, constraints =
                                           list(gammaHat[1]>=0, gammaHat[2:(L+1)] >= -ini_bound, gammaHat <= ini_bound)),
                         target = Problem(objective, constraints =
                                            list(gammaHat >= 0, gammaHat <= 1, sum(gammaHat)==1, gammaHat[1] == 1)),
                         weighted = Problem(objective, constraints =
                                              list(gammaHat[1]>=0, gammaHat[2:(L+1)] >= -ini_bound, gammaHat <= ini_bound)),
                         zero = Problem(objective, constraints =
                                          list(gammaHat == 0)))
        result1 <- solve(problem)
        if(result1$status %in% c('optimal', 'optimal_inaccurate')){
          gamma_step1_1 = as.vector(result1$getValue(gammaHat))
        }else{
          print('Failed!')
          gamma_step1_1 = c(1, rep(0, L))
        }
        beta_hat_step1_1 = B_hat_1 %*% gamma_step1_1
        if(base == 'mixT0'){
          gamma_weight = gamma_step1_1
          
          sample_mixT = sample(1:nrow(X_target), size = round((1-tau)*N_target))
          sample_0 = sample(1:nrow(X_target), size = round((tau)*N_target))
          X_ini = rbind(X_target[sample_mixT, ], X_target[sample_0, ])
          Y_ini = c(X_target[sample_mixT, ] %*% B_hat_1 %*% gamma_weight,
                    X_target[sample_0, ] %*% rep(0, p))
          
          beta = Variable(p)
          model = Problem(Minimize(mean((Y_ini - X_ini %*% beta)^2)))
          result <- solve(model)
          beta_hat_step1_1 = result$getValue(beta)
          mean((beta_hat_step1_1 - data$beta_valid)^2)
        }
        if(base == 'mixST'){
          gamma_source = gamma_step1_1
          
          gammaHat <- Variable(ncol(B_source_hat)+1)
          B_hat_1 = cbind(beta_target_hat, B_source_hat)
          objective <- Minimize(mean((X_target[test_rid, ] %*% B_hat_1 %*% gammaHat - 
                                        y_target[test_rid])^2))
          problem = Problem(objective, constraints =
                              list(gammaHat[1]>=0, gammaHat[2:(L+1)] >= -ini_bound, gammaHat <= ini_bound))
          result = solve(problem)
          gamma_weight = result$getValue(gammaHat)
          
          sample_mixST = sample(1:nrow(X_target), size = round((1-tau)*N_target))
          sample_mixsource = sample(1:nrow(X_target), size = round((tau)*N_target))
          X_ini = rbind(X_target[sample_mixST, ], X_target[sample_mixsource, ])
          Y_ini = c(X_target[sample_mixST, ] %*% B_hat_1 %*% gamma_weight,
                    X_target[sample_mixsource, ] %*% B_source_hat %*% gamma_source[-1])
          
          beta = Variable(p)
          model = Problem(Minimize(mean((Y_ini - X_ini %*% beta)^2)))
          result <- solve(model)
          beta_hat_step1_1 = result$getValue(beta)
          mean((beta_hat_step1_1 - data$beta_valid)^2)
        }
        if(base == 'target'){
          beta_hat_step1_1 = beta_target_only
        }
        
        result = list(`gamma` = gamma_step1_1, `beta` = beta_hat_step1_1)
        step1 = c(step1, list(result))
        result = Trans_label_est_split(fold = fold, tau,  X_source_matrix, base, B_source_hat, beta_target_hat,
                                       X_target, y_target, N_target, N_source, 
                                       beta_hat_step1_1, beta_target_only, which_beta_to_use)
        maximin_step1 = c(maximin_step1, list(result))
        # if(base == base_names[1]){
        #   print(paste0('s2_Q in target only:', round(result$ini_split$ini_target$s2_Q, 3)))
        #   print(paste0('s2_Q in linear comb:', round(result$ini_split$ini_linear_comb$s2_Q, 3)))
        #   print(paste0('s2_Q in robust linear comb:', round(result$ini_split$ini_robust_linear_comb$s2_Q, 3)))
        # }
      }
      names(maximin_step1) = names(step1) = base_names
      maximin_base_est = c(maximin_base_est, list(maximin_step1))
      base_est = c(base_est, list(step1))
    }
    
    maximin_base_ave_later_lst = base_lst = c()
    for(base in base_names){
      maximin_base_avg = base_avg = list()
      maximin_base_avg$gamma_hat = (maximin_base_est[[1]][[base]]$min_s2Q_robust_A_base_split$gamma_hat +
                                      maximin_base_est[[2]][[base]]$min_s2Q_robust_A_base_split$gamma_hat) / 2
      maximin_base_avg$beta_hat = (maximin_base_est[[1]][[base]]$min_s2Q_robust_A_base_split$beta_hat +
                                     maximin_base_est[[2]][[base]]$min_s2Q_robust_A_base_split$beta_hat) / 2
      base_avg$gamma_hat = (base_est[[1]][[base]]$gamma + base_est[[2]][[base]]$gamma) / 2
      base_avg$beta_hat = (base_est[[1]][[base]]$beta + base_est[[2]][[base]]$beta) / 2
      maximin_base_ave_later_lst = c(maximin_base_ave_later_lst, list(maximin_base_avg))
      base_lst = c(base_lst, list(base_avg))
    }
    names(base_lst) = paste0('base_', base_names)
    names(maximin_base_ave_later_lst) = paste0('GART_later_', base_names, '_', which_beta_to_use)
    
    maximin_base_ave_first_lst = c()
    beta_target_hat_1 = beta_target_hat_lst[[1]]
    beta_target_hat_2 = beta_target_hat_lst[[2]]
    for(base in names(base_lst)){
      beta_hat_step1_1 = base_lst[[base]]$beta
      base0 = str_remove(base, 'base\\_')
      result = Trans_label_est_split_ave(X_source_matrix, tau, B_source_hat, base0, beta_target_hat_1, beta_target_hat_2,
                                         X_target, y_target, N_target, N_source,
                                         beta_hat_step1_1, beta_target_only, which_beta_to_use)
      junk = list(`beta_hat` = result$min_s2Q_robust_A_robust_split$beta_hat,
                  `gamma_hat` = result$min_s2Q_robust_A_robust_split$gamma_hat)
      maximin_base_ave_first_lst = c(maximin_base_ave_first_lst, list(junk))
    }
    names(maximin_base_ave_first_lst) = paste0('GART_', base_names)
    
    beta_est = cbind(beta_est, sapply(maximin_base_ave_first_lst, function(x) x$beta_hat))
    gamma_est = cbind(gamma_est, sapply(maximin_base_ave_first_lst, function(x) x$gamma_hat)) 
  }
  
  
  
  ### evaluate
  message('Start evaluating...')
  coef_compare_est = compare_est(B_source_hat, X_source, y_source, X_target, y_target)
  beta_est = cbind(`true_beta` = data$beta_target,
                   beta_est,
                   coef_compare_est$beta
  )
  gamma_est = cbind(gamma_est, rbind(c(0, 0), coef_compare_est$gamma))
  beta_valid = data$valid$beta
  eval = evaluate_fun(beta_valid, beta_est, 
                              X_target_valid, y_target_valid)
  message('R2:')
  print(eval$R2)
  return(list(`beta_est` = beta_est, `gamma_est` = gamma_est,
              `eval` = eval))
}


evaluate_fun <- function(beta_valid, beta_est, 
                                 X_target_valid, y_target_valid){
  
  if(ncol(y_target_valid) == 1){
    epsilon = X_target_valid %*% (as.matrix(beta_est) - as.vector(beta_valid))
    mse = apply(epsilon, 2, function(x){
      mean(x^2)
    })
    predict_mse0 = as.vector(y_target_valid) - X_target_valid %*% as.matrix(beta_est)
    predict_mse = apply(predict_mse0, 2, function(x){
      mean(x^2)
    })
    R2 = apply(predict_mse0, 2, function(x){
      1 - mean(x^2) / var(y_target_valid)
    })
    beta_diff2_norm = sqrt(apply(beta_est[,-1] - beta_est[, 1], 2, function(x){
      sum(x^2)
    }))
  }else{
    mse = apply(beta_valid, 2, function(beta){
      epsilon = X_target_valid %*% (as.matrix(beta_est) - as.vector(beta))
      apply(epsilon, 2, function(x){
        mean(x^2)
      })
    })
    mse_mean = rowMeans(mse); mse_low = apply(mse, 1, min); mse_high = apply(mse, 1, max)
    mse = data.frame(`method` = names(mse_mean), `mean` = as.vector(mse_mean),
                     `min` = as.vector(mse_low), `max` = as.vector(mse_high))
    predict_mse = NULL
    R2 = apply(beta_est, 2, function(beta){
      predict_error = (y_target_valid - as.vector(X_target_valid %*% beta))^2
      predict_error = apply(predict_error, 2, mean)
      var_y = apply(y_target_valid, 2, var)
      mean(1 - predict_error / var_y)
    })
    
    beta_diff2_norm = apply(beta_est, 2, function(beta_hat){
      mean(colMeans((beta_valid - as.vector(beta_hat))^2))
    })
    
  }
  return(list(`mse` = mse,
              `R2` = R2,
              `predict_mse` = predict_mse,
              `beta_diff2_norm` = beta_diff2_norm)
  )
}




est_ini_w_target <- function(X_target_train, y_target_train,
                             X_target_test, y_target_test, beta_target_hat, L){
  n_train = length(y_target_train)
  n_test = length(y_target_test)
  s2_Q_hat = mean((y_target_test - X_target_test %*% beta_target_hat)^2)
  gamma_ini = rep(1/L, L)
  return(list(`beta` = beta_target_hat, `s2_Q` = s2_Q_hat, `gamma` = gamma_ini))
}

est_ini_w_linear_comb <- function(B_source_hat, beta_target_hat,
                                  X_target_train, y_target_train,
                                  X_target_test, y_target_test, L){
  n_train = length(y_target_train)
  n_test = length(y_target_test)
  ### Use gamma to estimate beta_hat_ini!
  gammaHat <- Variable(ncol(B_source_hat))
  objective <- Minimize(mean((X_target_train %*% B_source_hat %*% gammaHat - 
                                y_target_train)^2))
  problem <- Problem(objective, constraints =
                       list(gammaHat >= 1e-10, gammaHat <= 1, sum(gammaHat)==1))
  result <- solve(problem)
  if(result$status %in% c('optimal', 'optimal_inaccurate')){
    gamma_ini = as.vector(result$getValue(gammaHat))
    beta_hat_ini = B_source_hat %*% gamma_ini
    s2_Q_hat = mean((y_target_test - X_target_test %*% beta_hat_ini)^2)
  }else{
    gamma_ini = rep(1/L, L)
    beta_hat_ini = B_source_hat %*% gamma_ini
    s2_Q_hat = mean((y_target_test - X_target_test %*% beta_hat_ini)^2) 
  }
  return(list(`beta` = beta_hat_ini, `s2_Q` = s2_Q_hat, `gamma` = gamma_ini))
}


est_ini_w_robust_linear_comb <- function(B_source_hat, beta_target_hat,
                                         X_target_train, y_target_train,
                                         X_target_test, y_target_test, L){
  n_train = length(y_target_train)
  n_test = length(y_target_test)
  
  
  ### Use gamma to estimate beta_hat_ini!
  gammaHat <- Variable(ncol(B_source_hat)+1)
  B_hat = cbind(beta_target_hat, B_source_hat)
  
  objective <- Minimize(mean((X_target_test %*% B_hat %*% gammaHat - 
                                y_target_test)^2))
  problem <- Problem(objective, constraints =
                       list(gammaHat >= 1e-10, gammaHat <= 1, sum(gammaHat)==1))
  result <- solve(problem)
  
  if(result$status == 'optimal'){
    gamma_ini = as.vector(result$getValue(gammaHat))
    beta_hat_ini = B_hat %*% gamma_ini
    s2_Q_hat = mean((y_target_train - X_target_train %*% beta_hat_ini)^2)
  }else{
    gamma_ini = rep(1/L, L)
    beta_hat_ini = B_hat %*% gamma_ini
    s2_Q_hat = mean((y_target_train - X_target_train %*% beta_hat_ini)^2)
  }
  return(list(`beta` = beta_hat_ini, `s2_Q` = s2_Q_hat, `gamma` = gamma_ini))
}



est_ini <- function(fold=1, sample = c('split', 'cv'), nfold, X_target, y_target, 
                    Sigma_Q_hat, Gamma_beta_hat, B_source_hat, beta_target_hat, L){
  N_target = nrow(X_target)
  
  test_rid = (1:N_target)[((1:N_target) %% 2) == fold]
  # test_rid = 1:40
  train_rid = setdiff(1:N_target, test_rid)
  X_target_train = X_target[train_rid, ]
  y_target_train = y_target[train_rid]
  X_target_test = X_target[test_rid, ]
  y_target_test = y_target[test_rid]
  ini_target = est_ini_w_target(X_target_train, y_target_train,
                                X_target_test, y_target_test, beta_target_hat, L)
  ini_linear_comb = est_ini_w_linear_comb(B_source_hat, beta_target_hat, 
                                          X_target_train, y_target_train,
                                          X_target_test, y_target_test, L)
  ini_robust_linear_comb = est_ini_w_robust_linear_comb(B_source_hat, beta_target_hat,
                                                        X_target_train, y_target_train,
                                                        X_target_test, y_target_test, L)
  return(list(`test_rid` = test_rid,
              `ini_target` = ini_target,
              `ini_linear_comb` = ini_linear_comb,
              `ini_robust_linear_comb` = ini_robust_linear_comb
  )
  )
}


Trans_label_est_split <- function(fold, tau,  X_source_matrix, base, B_source_hat, beta_target_hat,
                                  X_target, y_target, N_target, N_source,
                                  beta_hat_step1, beta_target_only, which_beta_to_use){
  L = ncol(B_source_hat)
  ### estimate b(l), Sigma_Q, Gamma_beta
  Sigma_Q_hat = t(X_target) %*% X_target / nrow(X_target)
  Gamma_beta_hat = t(B_source_hat) %*% Sigma_Q_hat %*% B_source_hat
  
  ### estimate initialization of gamma_l, s_Q, beta
  ini_all = est_ini(fold = fold, sample = 'split', nfold = NA, X_target, y_target, 
                    Sigma_Q_hat, Gamma_beta_hat, B_source_hat, beta_target_hat, L)
  
  s2_Q_hat_ini_target = ini_all$ini_target$s2_Q
  gamma_hat_ini_target = ini_all$ini_target$gamma
  
  s2_Q_hat_ini_linear = ini_all$ini_linear_comb$s2_Q
  gamma_hat_ini_linear = ini_all$ini_linear_comb$gamma
  
  s2_Q_hat_ini_linear_robust = ini_all$ini_robust_linear_comb$s2_Q
  gamma_hat_ini_linear_robust = ini_all$ini_robust_linear_comb$gamma
  
  s2_Q_hat_min = min(s2_Q_hat_ini_linear, s2_Q_hat_ini_target) + tau
  
  
  if(which_beta_to_use == 'target_only'){
    beta_hat_ini_target = beta_target_only 
  }else{
    beta_hat_ini_target = ini_all$ini_target$beta
  }
  
  
  test_rid = ini_all$test_rid
  X_target_test = X_target[test_rid, ]
  y_target_test = y_target[test_rid]
  beta_hat_ini_robust = beta_hat_step1
  
  ### estiamte beta
  coef_min_s2Q_robust_1 = 
    Trans_label_robust_A_robust(X_target_test, y_target_test, base,
                                B_source_hat, Sigma_Q_hat, ini_all$ini_target$beta,
                                beta_hat_ini_target,  s2_Q_hat_min, beta_hat_ini_robust)
  
  coef_min_s2Q_robust_A_base = list(`gamma_hat` = coef_min_s2Q_robust_1$gamma_hat,
                                    `beta_hat` = coef_min_s2Q_robust_1$beta_hat,
                                    `success` = coef_min_s2Q_robust_1$success)
  result = list(
    `min_s2Q_robust_A_base_split` = coef_min_s2Q_robust_A_base,
    `ini_split` = ini_all
  )
  return(result)
}



Trans_label_est_split_ave <-  function(X_source_matrix, tau,  B_source_hat, base,
                                       beta_target_hat_1, beta_target_hat_2,
                                       X_target, y_target, N_target, N_source,
                                       beta_hat_step1, beta_target_only, which_beta_to_use){
  L = ncol(B_source_hat)
  ### estimate b(l), Sigma_Q, Gamma_beta
  Sigma_Q_hat = t(X_target) %*% X_target / nrow(X_target)
  Gamma_beta_hat = t(B_source_hat) %*% Sigma_Q_hat %*% B_source_hat
  
  ### estimate initialization of gamma_l, s_Q, beta
  ini_all = est_ini(fold = 1, sample = 'split', nfold = NA, X_target, y_target,
                    Sigma_Q_hat, Gamma_beta_hat, B_source_hat, beta_target_hat_2, L)
  
  s2_Q_hat_ini_target = ini_all$ini_target$s2_Q
  gamma_hat_ini_target = ini_all$ini_target$gamma
  
  s2_Q_hat_ini_linear = ini_all$ini_linear_comb$s2_Q
  gamma_hat_ini_linear = ini_all$ini_linear_comb$gamma
  
  s2_Q_hat_ini_linear_robust = ini_all$ini_robust_linear_comb$s2_Q
  gamma_hat_ini_linear_robust = ini_all$ini_robust_linear_comb$gamma
  
  s2_Q_hat_min = min(s2_Q_hat_ini_linear, s2_Q_hat_ini_target) + tau
  # s2_Q_hat_min = s2_Q_hat_ini_linear_robust + tau
  
  if(which_beta_to_use == 'target_only'){
    beta_hat_ini_target = beta_target_only
  }else{
    beta_hat_ini_target = ini_all$ini_target$beta
  }
  
  test_rid = ini_all$test_rid
  X_target_test = X_target[test_rid, ]
  y_target_test = y_target[test_rid]
  beta_hat_ini_robust = beta_hat_step1
  
  ### estiamte beta
  coef_min_s2Q_robust_1 =
    Trans_label_robust_A_robust(X_target_test, y_target_test, base,
                                B_source_hat, Sigma_Q_hat, ini_all$ini_target$beta,
                                beta_hat_ini_target,  s2_Q_hat_min, beta_hat_ini_robust)
  
  ### use the second half!
  ### estimate initialization of gamma_l, s_Q, beta
  ini_all = est_ini(fold = 0, sample = 'split', nfold = NA, X_target, y_target,
                    Sigma_Q_hat, Gamma_beta_hat, B_source_hat, beta_target_hat_1, L)
  
  s2_Q_hat_ini_target = ini_all$ini_target$s2_Q
  gamma_hat_ini_target = ini_all$ini_target$gamma
  
  s2_Q_hat_ini_linear = ini_all$ini_linear_comb$s2_Q
  gamma_hat_ini_linear = ini_all$ini_linear_comb$gamma
  
  s2_Q_hat_ini_linear_robust = ini_all$ini_robust_linear_comb$s2_Q
  gamma_hat_ini_linear_robust = ini_all$ini_robust_linear_comb$gamma
  
  s2_Q_hat_min = min(s2_Q_hat_ini_linear, s2_Q_hat_ini_target) + tau
  
  if(which_beta_to_use == 'target_only'){
    beta_hat_ini_target = beta_target_only
  }else{
    beta_hat_ini_target = ini_all$ini_target$beta
  }
  test_rid = ini_all$test_rid
  X_target_test = X_target[test_rid, ]
  y_target_test = y_target[test_rid]
  
  ### estiamte beta
  coef_min_s2Q_robust_2 =
    Trans_label_robust_A_robust(X_target_test, y_target_test, base,
                                B_source_hat, Sigma_Q_hat, ini_all$ini_target$beta,
                                beta_hat_ini_target,  s2_Q_hat_min, beta_hat_ini_robust)
  
  
  coef_min_s2Q_robust_A_robust = list(`gamma_hat` = (coef_min_s2Q_robust_1$gamma_hat + coef_min_s2Q_robust_2$gamma_hat)/2,
                                      `beta_hat` = (coef_min_s2Q_robust_1$beta_hat + coef_min_s2Q_robust_2$beta_hat)/2,
                                      `success` = (coef_min_s2Q_robust_1$success + coef_min_s2Q_robust_2$success)/2)
  result = list(
    `min_s2Q_robust_A_robust_split` = coef_min_s2Q_robust_A_robust,
    `ini_split` = ini_all
  )
  return(result)
}


Trans_label_robust_A_robust <- function(X_target_test, y_target_test, base,
                                        B_source_hat, Sigma_Q_hat, beta_hat_ini_old,
                                        beta_hat_ini, s2_Q_hat, beta_hat_ini_robust){
  B_hat_old = as.matrix(cbind(beta_hat_ini_old, B_source_hat))
  B_hat_new = as.matrix(cbind(beta_hat_ini, B_source_hat))
  B_hat_new_fun = as.matrix(cbind(beta_hat_ini, 
                                  B_source_hat)) - as.vector(beta_hat_ini_robust)
  
  Gamma_beta_hat_new = t(B_hat_new_fun) %*% Sigma_Q_hat %*% B_hat_new_fun
  
  gammaHat <- Variable(ncol(B_hat_new))
  objective <- Minimize(quad_form(gammaHat,
                                  Gamma_beta_hat_new)
                        # + sum(gammaHat/c(N_target, N_source))
  )
  constraints0 = mean((X_target_test %*% B_hat_old %*% gammaHat -
                         y_target_test)^2) - s2_Q_hat <= 0
  problem = Problem(objective, constraints =
                      list(constraints0, gammaHat >= 0, gammaHat <= 1, sum(gammaHat)==1))
  
  
  result <- solve(problem)
  if(result$status %in% c('optimal', 'optimal_inaccurate')){
    gamma_hat = as.vector(result$getValue(gammaHat))
    beta_hat = B_hat_new %*% gamma_hat
    success = 1
  }else{
    gamma_hat = beta_hat = NA
    success = 0
    print('Failed!')
  }
  return(list(`gamma_hat` = gamma_hat,
              `beta_hat` = beta_hat,
              `success` = success))
}




compare_est <- function(B_source_hat, X_source, y_source, X_target, y_target){
  L = length(X_source)
  X_source_matrix = do.call(rbind, X_source)
  y_source_vec = do.call(rbind, y_source)
  
  Sigma_Q_hat = t(X_target) %*% X_target / nrow(X_target)
  Gamma_beta_hat = t(B_source_hat) %*% Sigma_Q_hat %*% B_source_hat
  n_target = length(y_target)
  
  
  gammaHat <- Variable(ncol(B_source_hat))
  objective <- Minimize(mean((X_target %*% B_source_hat %*% gammaHat -
                                y_target)^2))
  problem <- Problem(objective, constraints =
                       list(gammaHat >= 1e-10, sum(gammaHat)==1))
  result <- solve(problem)
  
  gamma_linear_source = as.vector(result$getValue(gammaHat))
  beta_linear_source = B_source_hat %*% gamma_linear_source
  
  
  gammaHat <- Variable(ncol(B_source_hat))
  objective <- Minimize(quad_form(gammaHat, Gamma_beta_hat))
  problem <- Problem(objective, constraints =
                       list(gammaHat >= 1e-10, sum(gammaHat)==1))
  result <- solve(problem)
  gamma_maximin = as.vector(result$getValue(gammaHat))
  beta_maximin = B_source_hat %*% gamma_maximin
  
  out = cv.glmnet(X_target, y_target)
  out = glmnet(X_target, y_target,
               lambda = out$lambda.min)
  beta_target = out$beta
  
  source_data = sapply(1:L, function(i){
    list(`x` = X_source[[i]], `y` = y_source[[i]])
  }, simplify = FALSE)
  target_data = list(`x` = X_target, `y` = y_target)
  out_transglm = glmtrans(target = target_data,  source = source_data, family = "gaussian",
                          intercept = FALSE, detection.info = FALSE)
  beta_transglm = out_transglm$beta[-1]
  
  n.vec = c(n_target, sapply(X_source, FUN = nrow))
  prop.re1 <- tryCatch(Trans.lasso(rbind(X_target, X_source_matrix),
                                   c(y_target, y_source_vec), n.vec,
                                   I.til = 1:round(n_target/2), l1 = T), error = function(e) NULL)
  
  prop.re2 <- tryCatch(Trans.lasso(rbind(X_target, X_source_matrix),
                                   c(y_target, y_source_vec), n.vec,
                                   I.til = (round(n_target/2)+1):n_target, l1=T), error = function(e) NULL)
  if(!is.null(prop.re1) & !is.null(prop.re2)){
    beta.prop <- (prop.re1$beta.hat + prop.re2$beta.hat) / 2
    gamma.prop = (prop.re1$theta.hat + prop.re2$theta.hat) / 2
  }else{
    beta.prop = rep(NA, length(beta_target))
    gamma.prop = NA
  }
  result = list(`beta` = data.frame(
    'target' = as.vector(beta_target),
    'linear_source_beta' = as.vector(beta_linear_source),
    'trans_lasso' = as.vector(beta.prop),
    'trans_glm' = as.vector(beta_transglm),
    'maximin' = as.vector(beta_maximin)
  ),
  `gamma` = data.frame('maximin' = gamma_maximin,
                       'linear_source' = gamma_linear_source)
  )
  return(result)
}




