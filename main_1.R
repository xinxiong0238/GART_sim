rm(list = ls())
library(LaplacesDemon)
library(glmnet)
library(ggplot2)
library(dplyr)
library(CVXR)
library(glmtrans)
library(stringr)

source('helper/function_GART.R')
source('helper/function_translasso.R')
s_Q = 1; s_P = 0.5
eval_model = c()
tau = 1/200
ratio_target_n = 0.01

if (!dir.exists('out')) {
  dir.create('out')
}
if (!dir.exists('out/setting1')) {
  dir.create('out/setting1')
}
if (!dir.exists('out/setting2')) {
  dir.create('out/setting2')
}
if (!dir.exists('out/setting3')) {
  dir.create('out/setting3')
}

for(sett in c(1,2,3)){
  for(jj in 1:100){
    if(sett == 1){
      ## Setting 1
      sd_source_vec = c(0, 0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1, 1.3, 1.6, 1.8, 2, 2.5, 3)
      eval_model = data.frame()
      source('helper/generate_sett1.R')
      eval_model_sd = c()
      for(sd_source in sd_source_vec){
        set.seed(jj)
        message(paste0('Setting ', sett, ' rep:', jj, ' par:', sd_source))
        data0 = simu_eqSource(sd_source, L = 8, ratio_target_n, s_Q, s_P)
        result = GART_fun(data0, tau = tau)
        eval = evaluate_fun(beta_valid = data0$valid$beta,
                            beta_est = result$beta_est,
                            X_target_valid = data0$valid$X,
                            y_target_valid = data0$valid$y)
        eval_model_sd = rbind(eval_model_sd,
                              data.frame(`par` = sd_source,
                                         `par_nm` = 'kappa',
                                         `tau` = tau,
                                         `setting` = 'setting1',
                                         `method` = names(eval$mse),
                                         `eval` = as.vector(eval$mse),
                                         `eval_type` = 'mse',
                                         `rep` = jj))
        saveRDS(eval_model_sd, file = paste0('out/setting1/part_',jj,'.RDS'))
      }
    }
    
    if(sett == 2){
      ## Change alpha
      source('helper/generate_sett2.R')
      eval_model_alpha = c()
      alpha_vec = c(0.001, 0.01, 0.02, 0.05, 0.075, 0.1, 0.12, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75)
      for(alpha in alpha_vec){
        message(paste0('Setting ', sett, ' rep:', jj, ' par:', alpha))
        set.seed(jj)
        data0 = simu_bv(alpha = alpha, ratio_target_n = ratio_target_n,
                        s_Q=s_Q, s_P=s_P)
        result = GART_fun(data0, tau = tau)
        eval = evaluate_fun(beta_valid = data0$valid$beta, 
                            beta_est = result$beta_est,
                            X_target_valid = data0$valid$X,
                            y_target_valid = data0$valid$y)
        eval_model_alpha = rbind(eval_model_alpha,
                                 data.frame(`par` = alpha,
                                            `par_nm` = 'alpha',
                                            `tau` = tau,
                                            `setting` = 'setting2',
                                            `method` = names(eval$mse),
                                            `eval` = as.vector(eval$mse),
                                            `eval_type` = 'mse',
                                            `rep` = jj))
        
        saveRDS(eval_model_alpha, file = paste0('out/setting2/part_',jj,'.RDS'))
      }
    }
    
    if(sett == 3){
      ### change sparsity
      source('helper/generate_sett3.R')
      # high D: fix l0, change l1 norm
      L = 10; p = 150; s0 = 50
      s1_vec = c(1, 1.5, 2, 3, 4, 7, 10, 15, 20, 22.5, 25, 27.5, 30, 35)
      eval_model_s1 = c()
      for(s1 in s1_vec){
        set.seed(jj)
        message(paste0('Setting ', sett, ' rep:', jj, ' par:', s1))
        data0 = simu_highD_sparseT(p, s0, s1, L, ratio_target_n, s_Q=1, s_P=.5)
        result = GART_fun(data0, tau = tau)
        eval = evaluate_fun(beta_valid = data0$valid$beta,
                            beta_est = result$beta_est,
                            X_target_valid = data0$valid$X,
                            y_target_valid = data0$valid$y)
        eval_model_s1 = rbind(eval_model_s1,
                              data.frame(`par` = s1,
                                         `par_nm` = 's1',
                                         `tau` = tau,
                                         `setting` = 'setting3',
                                         `method` = names(eval$mse),
                                         `eval` = as.vector(eval$mse),
                                         `eval_type` = 'mse',
                                         `rep` = jj))
        saveRDS(eval_model_s1, file = paste0('out/setting3/part_',jj,'.RDS'))
        
      }
    }
  }
}
