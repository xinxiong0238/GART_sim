rm(list = ls())
library(LaplacesDemon)
library(glmnet)
library(ggplot2)
library(dplyr)
library(CVXR)
library(glmtrans)
library(stringr)

source('helper/function_GART_more_ini.R')
source('helper/function_translasso.R')
s_Q = 1; s_P = 0.5
ratio_target_n = 0.01

if (!dir.exists('out')) {
  dir.create('out')
}
if (!dir.exists('out/setting4.1')) {
  dir.create('out/setting4.1')
}
if (!dir.exists('out/setting4.2')) {
  dir.create('out/setting4.2')
}

for(jj in 1:100){
  eval_model_1 =  eval_model_2 = c()
  source('helper/generate_sett2.R')
  alpha_vec = rev(c(0.1, 0.2, 0.4))
  for(alpha in alpha_vec){
    set.seed(jj)
    data0 = simu_bv(alpha = alpha, ratio_target_n = ratio_target_n, s_Q=s_Q, s_P=s_Q, 
                    valid_type = 0, beta_valid_target_weight = 1)
    for(tau in c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)){
      message(paste0('rep:', jj, ' alpha:', alpha, ' tau:', tau))
      result = GART_function_more_ini(data0, tau)
      
      ### setting 4.1
      for(beta_valid_target_weight in c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)){
        set.seed(jj * 5000)
        data = simu_bv(alpha = alpha, ratio_target_n = ratio_target_n,
                       s_Q=1, s_P=.5, valid_type = 1,
                       beta_valid_target_weight = beta_valid_target_weight)
        eval = evaluate_fun(beta_valid = data$valid$beta, 
                                    beta_est = result$beta_est,
                                    X_target_valid = data$valid$X,
                                    y_target_valid = data$valid$y)
        method_nm = names(eval$R2)
        eval_model_1 = rbind(eval_model_1,
                                 data.frame(`par` = alpha,
                                            `par_nm` = 'alpha',
                                            `tau` = tau,
                                            `setting` = 'setting4.1',
                                            `beta_valid_target_weight` = beta_valid_target_weight,
                                            `method` = method_nm,
                                            `eval` = c(eval$mse, eval$R2),
                                            `eval_type` = c(rep('mse', length(method_nm)), 
                                                            rep('R2', length(method_nm))),
                                            `rep` = jj))
      }
      saveRDS(eval_model_1, file = paste0('out/setting4.1/part_',jj,'.RDS'))
      
      ### setting 4.2
      for(beta_valid_target_weight in c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)){
        set.seed(jj * 5000)
        data = simu_bv(alpha = alpha, ratio_target_n = ratio_target_n,
                       s_Q=1, s_P=.5, valid_type = 2,
                       beta_valid_target_weight = beta_valid_target_weight)
        eval = evaluate_fun(beta_valid = data$valid$beta, 
                            beta_est = result$beta_est,
                            X_target_valid = data$valid$X,
                            y_target_valid = data$valid$y)
        method_nm = names(eval$R2)
        eval_model_2 = rbind(eval_model_2,
                                 data.frame(`par` = alpha,
                                            `par_nm` = 'alpha',
                                            `tau` = tau,
                                            `setting` = 'setting4.2',
                                            `beta_valid_target_weight` = beta_valid_target_weight,
                                            `method` = method_nm,
                                            `eval` = c(eval$mse, eval$R2),
                                            `eval_type` = c(rep('mse', length(method_nm)),
                                                            rep('R2', length(method_nm))),
                                            `rep` = jj))
      }
      saveRDS(eval_model_2, file = paste0('out/setting4.2/part_',jj,'.RDS'))
    }
  }
}
