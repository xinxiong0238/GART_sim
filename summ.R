rm(list = ls())
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)
library(RColorBrewer)
scaleFUN <- function(x) paste0(sprintf("%.1f", x*100),'%')
manual_color = c(brewer.pal(8, 'Set1')[1], 'violet', 'blue', 'cyan1', 'forestgreen', 'lawngreen',  brewer.pal(8, 'Set1')[4:5], 'lightblue', brewer.pal(8, 'Set1')[7:8])
manual_value = c('GART weighted', 'GART convex', 'GART mixST', 'GART source', 'GART mixT0', 'GART zero',
                 'Source mixture', 'Target only', 'TransGLM', 'TransLasso', 'Maximin')
if (!dir.exists('plot')) {
  dir.create('plot')
}
eval_model = c()
folders = list.files('out/')
for(folder in folders){
  print(folder)
  for(file in list.files(paste0('out/', folder))){
    sub_eval = readRDS(paste0('out/', folder, '/', file))
    if(folder %in% c('setting1', 'setting2', 'setting3')){
      sub_eval$beta_valid_target_weight = NA
    }
    eval_model = rbind(eval_model, sub_eval)
  }
}


eval_model_summary = eval_model %>%
  group_by(par, par_nm, tau, setting, beta_valid_target_weight, method, eval_type) %>%
  summarise('mean_eval_est' = mean(eval, na.omit = TRUE),
            'sd_eval_est' = sd(eval),
            'upper_eval' = quantile(eval, 0.975),
            'lower_eval' = quantile(eval, 0.025),
            'rep_N' = n())
table(eval_model_summary$method)


eval_model_summary$method = factor(eval_model_summary$method, 
                                     levels =  c('GART_weighted', 'GART_convex', 'GART_source', 'GART_mixST', 'GART_mixT0' , 'GART_zero',
                                                 'linear_source_beta', 'target', 'trans_glm',
                                                 'trans_lasso', 'maximin'),
                                   labels = c('GART weighted', 'GART convex', 'GART source', 'GART mixST', 'GART mixT0', 'GART zero', 
                                              'Source mixture', 'Target only', 'TransGLM', 'TransLasso', 'Maximin'))
scaleFUN2 <- function(x) paste0(sprintf("%.3f", x))

eval_model_summary$par[eval_model_summary$setting == 'setting2'] = 
  eval_model_summary$par[eval_model_summary$setting == 'setting2'] * sqrt(35)

table(eval_model_summary$method)
p1 = ggplot(eval_model_summary %>%
              filter(setting %in% c('setting1'),
                     method %in% c('GART weighted', 'Source mixture',
                                   'Target only', 'TransGLM', 'TransLasso', 'Maximin'), 
                     par > 0.01) %>%
              filter(round(par, 2) %in% c(0.05, 0.1, 0.3, 1, 1.8, 3)),
            aes(x = par, y = log(mean_eval_est), color = method, group = paste(method, par))) +
  geom_point(position = position_dodge(width = .1)) +
  geom_errorbar(aes(ymax = log(upper_eval), ymin = log(lower_eval)),
                position = position_dodge(width = .1), 
                width = .2, size = 0.8) +
  scale_fill_manual(breaks = manual_value, values = manual_color) +
  scale_color_manual(breaks = manual_value, values = manual_color) +
  xlab(expression(kappa)) +
  ylab(expression(logMSE)) +
  scale_x_log10() +
  facet_wrap(setting~., scales = 'free', nrow = 1) +
  theme_bw()+
  theme(legend.position = 'right') +
  theme(legend.key.width = unit(1.2,"cm")) +
  theme(text = element_text(size=16))
p1

p2 = ggplot(eval_model_summary %>% filter(setting %in% c('setting2'),
                                          method %in% c('GART weighted', 'Source mixture',
                                                        'Target only', 'TransGLM', 'TransLasso', 'Maximin'), 
                                          par > 0.01) %>%
              filter(round(par,1) %in% c(0, 0.1, 0.3, 0.7, 1, 1.5, 2.5, 3)),
            aes(x = par, y = log(mean_eval_est), color = method, group = paste(method, par))) +
  geom_point(position = position_dodge(width = .1)) +
  geom_errorbar(aes(ymax = log(upper_eval), ymin = log(lower_eval)),
                position = position_dodge(width = .1),
                width = .2, size = 0.8) +
  # geom_smooth(se = F) +
  scale_fill_manual(breaks = manual_value, values = manual_color) +
  scale_color_manual(breaks = manual_value, values = manual_color) +
  xlab(expression(kappa)) +
  ylab(expression(logMSE)) +
  scale_x_log10() +
  facet_wrap(setting~., scales = 'free', nrow = 1) +
  theme_bw()+
  theme(legend.position = 'right') +
  theme(legend.key.width = unit(1.2,"cm")) +
  theme(text = element_text(size=16))
p2

p12 = ggarrange(p1, p2, common.legend = T, legend = 'right')

p3 = ggplot(eval_model_summary %>% filter(setting %in% c('setting3'),
                                          method %in% c('GART weighted', 'Source mixture',
                                                        'Target only', 'TransGLM', 'TransLasso', 'Maximin'), 
                                          par %in% c(1.5, 7, 15, 20, 25, 30, 35)),
            aes(x = par, y = log(mean_eval_est), color = method, group = paste(method, par))) +
  geom_point(position = position_dodge(width = 3)) +
  geom_errorbar(aes(ymax = log(upper_eval), ymin = log(lower_eval)),
                position = position_dodge(width = 3),
                width = 3, size = 0.8) +
  scale_fill_manual(breaks = manual_value, values = manual_color) +
  scale_color_manual(breaks = manual_value, values = manual_color) +
  xlab(expression(kappa)) +
  ylab(expression(logMSE)) +
  facet_wrap(setting~., scales = 'free', nrow = 1) +
  theme_bw()+
  theme(legend.position = 'right') +
  theme(legend.key.width = unit(1.2,"cm")) +
  theme(text = element_text(size=16))
p3




p_ini = ggplot(eval_model_summary %>% filter(setting %in% c('setting1','setting2'),
                                             method %in% c('GART weighted', 'GART convex', 'GART zero', 'GART source'),
                                             par > 0.01) %>%
                 filter(round(par,1) %in% c(0, 0.1, 0.3, 0.7, 1, 1.5, 3)),
               aes(x = par, y = log(mean_eval_est), color = method, group = paste(method, par))) +
  geom_point(position = position_dodge(width = .1)) +
  geom_errorbar(aes(ymax = log(upper_eval), ymin = log(lower_eval)),
                position = position_dodge(width = .1),
                width = .2, size = 0.8) +
  scale_color_manual(breaks = manual_value, values = manual_color) +
  xlab(expression(kappa)) +
  ylab(expression(logMSE)) +
  scale_x_log10() +
  facet_wrap(setting~., scales = 'free', nrow = 1) +
  theme_bw()+
  theme(legend.position = 'right') +
  theme(legend.key.width = unit(1.2,"cm")) +
  theme(text = element_text(size=16))
p_ini


pdf('plot/figure3.pdf', width = 12, height = 5)
print(p12)
dev.off()

pdf('plot/figure1_supp.pdf', width = 7, height = 5)
print(p3)
dev.off()

pdf('plot/figure3_supp.pdf', width = 11, height = 5)
print(p_ini)
dev.off()


eval_model_alpha_summ = eval_model_summary %>%
  filter(setting %in% c('setting4.1', 'setting4.2'),
         !is.na(method))

eval_model_alpha_summ$method_type = as.character(eval_model_alpha_summ$method)
eval_model_alpha_summ$method_type[!str_detect(eval_model_alpha_summ$method, 'GART')] = 'other'
eval_model_alpha_summ$method_type[str_detect(eval_model_alpha_summ$method, 'GART weighted')] = 'other'
eval_model_alpha_summ$tau[eval_model_alpha_summ$method_type == 'other' & eval_model_alpha_summ$method != 'GART_weighted'] = 0.01
eval_model_alpha_summ = eval_model_alpha_summ %>% 
  group_by(par, par_nm, tau, setting, beta_valid_target_weight, method, eval_type, method_type) %>%
  summarise(mse = mean(mean_eval_est),
            upper_eval = mean(upper_eval),
            lower_eval = mean(lower_eval))
eval_model_alpha_filter = eval_model_alpha_summ %>% 
  filter(! (str_detect(method, 'GART weighted') & tau != 0.01))

eval_model_alpha_filter$method_type[eval_model_alpha_filter$method_type != 'other'] = 
  paste0(str_replace(as.character(eval_model_alpha_filter$method[eval_model_alpha_filter$method_type != 'other']), ' ', '['), ']')

eval_model_alpha_filter$par = factor(eval_model_alpha_filter$par, levels = c(0.01, 0.1, 0.2, 0.4),
                                     labels = c(expression(kappa==0.01), expression(kappa==0.1),
                                                expression(kappa==0.2), expression(kappa==0.4)))

method_filter = 'GART source|GART zero|GART convex'
p4_1 = ggplot(eval_model_alpha_filter %>% 
                filter(tau %in% c(0.01, 0.1, 0.3, 0.5, 0.7),
                       !str_detect(method, method_filter),
                       str_detect(par, '0\\.2|0\\.4'),
                       setting == 'setting4.1',
                       eval_type == 'R2',
                       beta_valid_target_weight %in% c(0.1, 0.3, 0.5)) %>%
                mutate(tau = as.character(tau)), 
              aes(x = as.character(beta_valid_target_weight), y = mse, group = paste(method, tau))) +
    geom_bar(stat = 'identity', width = 0.8,
             aes(fill = method, alpha = tau), col = 'black',
             position = position_dodge()) +
    scale_alpha_manual(breaks = c('0.01', '0.1', '0.3', '0.5', '0.7'),
                       values = c(0.99, 0.9, 0.7, 0.5, 0.3)) +
    geom_point(position = position_dodge(width = 0.8), size = 1) +
    geom_errorbar(aes(ymax = upper_eval, ymin = lower_eval),
                  position = position_dodge(width = 0.8),
                  width = 0.3, size = 0.8, alpha = 0.8) +
    facet_grid(par~method_type, scales = 'free', labeller = label_parsed) +
    scale_color_manual(breaks = manual_value, values = manual_color) +
    scale_fill_manual(breaks = manual_value, values = manual_color) +
    xlab(expression(nu)) +
    ylab(expression(R^2)) +
    labs(shape = expression(tau)) +
    # ggtitle(title_nm) +
    # scale_y_log10() +
    theme_bw() +
    theme(legend.position = 'right') +
    theme(legend.key.width = unit(1.2,"cm")) +
    theme(text = element_text(size=16))

p4_2 = ggplot(eval_model_alpha_filter %>% 
                filter(tau %in% c(0.01, 0.1, 0.3, 0.5, 0.7),
                       !str_detect(method, method_filter),
                       str_detect(par, '0\\.2|0\\.4'),
                       setting == 'setting4.2',
                       eval_type == 'R2',
                       beta_valid_target_weight %in% c(0.1, 0.3, 0.5)) %>%
                mutate(tau = as.character(tau)), 
              aes(x = as.character(beta_valid_target_weight), y = mse, group = paste(method, tau))) +
  geom_bar(stat = 'identity', width = 0.8,
           aes(fill = method, alpha = tau), col = 'black',
           position = position_dodge()) +
  scale_alpha_manual(breaks = c('0.01', '0.1', '0.3', '0.5', '0.7'),
                     values = c(0.99, 0.9, 0.7, 0.5, 0.3)) +
  geom_point(position = position_dodge(width = 0.8), size = 1) +
  geom_errorbar(aes(ymax = upper_eval, ymin = lower_eval),
                position = position_dodge(width = 0.8),
                width = 0.3, size = 0.8, alpha = 0.8) +
  facet_grid(par~method_type, scales = 'free', labeller = label_parsed) +
  scale_color_manual(breaks = manual_value, values = manual_color) +
  scale_fill_manual(breaks = manual_value, values = manual_color) +
  xlab(expression(nu)) +
  ylab(expression(R^2)) +
  labs(shape = expression(tau)) +
  # ggtitle(title_nm) +
  # scale_y_log10() +
  theme_bw() +
  theme(legend.position = 'right') +
  theme(legend.key.width = unit(1.2,"cm")) +
  theme(text = element_text(size=16))


eval_model_alpha_summ = eval_model_summary %>%
  filter(setting %in% c('setting4.1', 'setting4.2'),
         !is.na(method))
eval_model_alpha_filter = eval_model_alpha_summ
eval_model_alpha_filter$method_type = 
  paste0(str_replace(as.character(eval_model_alpha_filter$method), ' ', '['), ']')
eval_model_alpha_filter$method_type = 
  factor(eval_model_alpha_filter$method_type,
         levels = c('GART[weighted]', 'GART[convex]', 
                    'GART[mixST]', 'GART[source]', 
                    'GART[mixT0]', 'GART[zero]'))

scaleFUN <- function(x) sprintf("%.2g", x)
p4_1_ini = ggplot(eval_model_alpha_filter %>% 
              filter(tau %in% c(0.01, 0.1, 0.3, 0.5, 0.7),
                     str_detect(par, '0\\.2|0\\.4'),
                     setting == 'setting4.1',
                     str_detect(method, 'GART'),
                     beta_valid_target_weight %in% c(0.1, 0.3, 0.5),
                     eval_type == 'R2') %>%
              mutate(tau = as.character(tau)), 
            aes(x = as.character(beta_valid_target_weight), y = mean_eval_est, group = paste(method, tau))) +
  geom_bar(stat = 'identity', width = 0.8,
           aes(fill = method, alpha = tau), col = 'black',
           position = position_dodge()) +
  scale_alpha_manual(breaks = c('0.01', '0.1', '0.3', '0.5', '0.7'),
                     values = c(0.99, 0.9, 0.7, 0.5, 0.3)) +
  geom_point(position = position_dodge(width = 0.8), size = 1) +
  geom_errorbar(aes(ymax = upper_eval, ymin = lower_eval),
                position = position_dodge(width = 0.8),
                width = 0.3, size = 0.8, alpha = 0.8) +
  facet_grid(par~method_type, scales = 'free', labeller = label_parsed) +
  scale_color_manual(breaks = manual_value, values = manual_color) +
  scale_fill_manual(breaks = manual_value, values = manual_color) +
  xlab(expression(nu)) +
  ylab(expression(R^2)) +
  labs(shape = expression(tau)) +
  theme_bw() +
  theme(legend.position = 'right') +
  theme(legend.key.width = unit(1.2,"cm")) +
  theme(text = element_text(size=16))
p4_1_ini


p4_2_ini = ggplot(eval_model_alpha_filter %>% 
                    filter(tau %in% c(0.01, 0.1, 0.3, 0.5, 0.7),
                           str_detect(par, '0\\.2|0\\.4'),
                           setting == 'setting4.2',
                           str_detect(method, 'GART'),
                           beta_valid_target_weight %in% c(0.1, 0.3, 0.5),
                           eval_type == 'R2') %>%
                    mutate(tau = as.character(tau)), 
                  aes(x = as.character(beta_valid_target_weight), y = mean_eval_est, group = paste(method, tau))) +
  geom_bar(stat = 'identity', width = 0.8,
           aes(fill = method, alpha = tau), col = 'black',
           position = position_dodge()) +
  scale_alpha_manual(breaks = c('0.01', '0.1', '0.3', '0.5', '0.7'),
                     values = c(0.99, 0.9, 0.7, 0.5, 0.3)) +
  geom_point(position = position_dodge(width = 0.8), size = 1) +
  geom_errorbar(aes(ymax = upper_eval, ymin = lower_eval),
                position = position_dodge(width = 0.8),
                width = 0.3, size = 0.8, alpha = 0.8) +
  facet_grid(par~method_type, scales = 'free', labeller = label_parsed) +
  scale_color_manual(breaks = manual_value, values = manual_color) +
  scale_fill_manual(breaks = manual_value, values = manual_color) +
  xlab(expression(nu)) +
  ylab(expression(R^2)) +
  labs(shape = expression(tau)) +
  theme_bw() +
  theme(legend.position = 'right') +
  theme(legend.key.width = unit(1.2,"cm")) +
  theme(text = element_text(size=16))
p4_2_ini



pdf('plot/figure4.pdf', width = 14, height = 5.5)
print(p4_1)
dev.off()


pdf('plot/figure2_supp.pdf', width = 14, height = 5.5)
print(p4_2)
dev.off()

pdf('plot/figure4_supp.pdf', width = 14, height = 5.5)
print(p4_1_ini)
print(p4_2_ini)
dev.off()


