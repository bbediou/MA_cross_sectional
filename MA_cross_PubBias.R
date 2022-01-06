######################################################################################
# Publication Bias analyses
######################################################################################
## Notes: Pub Bias using regression methods tests for presence of Small Study Effect
##   - the Significance of the slope coeff. is the small-study effect detection test.   
##   - the Intercept should be included and gives the "adjusted" effect size.   
######################################################################################
library(tidyverse) # for data-cleaning
library(kableExtra) # for making formatted table
library(robumeta) # for original implementation of RVE
library(metafor) # for multivariate meta-analytic working models 
library(clubSandwich) # for RVE with metafor; requires version 0.5.1 or higher
library(DT) # for datatable and formatRound
# library(patchwork) # to combine ggplots  
library(puniform)# p-uniform
library(weightr)# for 3-PSM

setwd("/Users/bediou/GoogleDrive/Meta-Analysis/DATA/ANALYSIS/2021/MA_cross-sectional/")
load("MA_data_cross.RData")
source("/Users/bediou/GoogleDrive/Meta-Analysis/DATA/ANALYSIS/2021/MA_cross-sectional/tidy_functions.R")
source("/Users/bediou/GoogleDrive/Meta-Analysis/DATA/ANALYSIS/2021/MA_cross-sectional/bootstrapping_functs.R")

rho = 0.8

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Egger's Sandwich test using modified precision covariate (Pustejovsky & Rodgers, 2019) 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#--- NULL MODELS ---------------------------------------------------------------
## RVE version
egg_sand <- robu(formula = es_win ~ 1 + sda, # Modified covariate for Egger's Regression test. 
                 var.eff.size = g_var, # ES variance
                 studynum = Paper,# clustering
                 modelweights = "CORR",# Correlated weights
                 # small = TRUE,# small sample correction ON => check df's if <4 then use smaller alpha to test significance. 
                 data = es_data_win)

## MULTILEVEL model
# variance-covariance matrix
vcov_mat <- impute_covariance_matrix(es_data_win$g_var, 
                                     cluster = es_data_win$Paper,                            
                                     r = rho) 
# multilevel model
egg_mlma <- rma.mv(es_win ~ 1  + sda, # inclusion of modified covariate (Egger's regression)
                   V = vcov_mat, # variance covariance matrix defined above
                   random = ~ 1 | Paper / ES_ID, # random factors (nesting)
                   sparse = TRUE, # 
                   data = es_data_win,
                   test = "t" # small sample correction isn't needed if doing CHE model with RVE
                   )
# Cluster robust CIs and DFs are obtained using a the following command 
# clubSandwich::coef_test(egg_mlma, vcov = "CR2")
# embedded in the function tidy_CHE 


#--- FULL MODELS ---------------------------------------------------------------
## RVE model with correlated weights  
egg_sand_full <- robu(formula = es_win ~ 1 + sda +# Modified covariate for Egger's Regression test. 
                        Cognitive_domain + # PRIMARY moderator
                        DV_type + # control moderator
                        Effect + # control moderator
                        Recruitment,  # control moderator, 
                      var.eff.size = g_var, # ES variance
                      studynum = Paper,# clustering
                      modelweights = "CORR",# Correlated weights
                      # small = TRUE,# small sample correction ON => check df's if <4 then use smaller alpha to test significance. 
                      data = es_data_win)
## MULTILEVEL
egg_mlma_full <- rma.mv(es_win ~ 1  + 
                          sda + # New precision estimate
                          Cognitive_domain + # PRIMARY moderator
                          DV_type + # control moderator
                          Effect + # control moderator
                          Recruitment, # control moderator, 
                        V = vcov_mat, 
                        random = ~ 1 | Paper / ES_ID, # inclusion of modified covariate (Egger's regression)
                        sparse = TRUE,
                        data = es_data_win,
                        test = "t" # This should be put on for small sample correction (but isn't needed if doing CHE model with robust standard errors (RVE))
)



# compare slopes
egger_tab <- bind_rows(cbind("method" = "Egger Sandwich (RVE)","model" = "NULL", tidy_robu(egg_sand)),
          cbind("method" = "Egger CHE (MLMA + RVE)", "model" = "NULL", tidy_CHE(egg_mlma)),
          cbind("method" = "Egger Sandwich (RVE)","model" = "FULL", tidy_robu(egg_sand_full)),
          cbind("method" = "Egger CHE (MLMA + RVE)", "model" = "FULL", tidy_CHE(egg_mlma_full))) %>%
  filter(term == "sda" | term == "intrcpt") %>%# keep only the slope & intercept
  arrange(method, model, term)
# egger_tab %>% 
#   kable() %>%
#   kable_styling(bootstrap_options = c("condensed", "striped", "hover"), full_width = FALSE, position = "left") %>%
#   row_spec(which(egger_tab$p<.05), bold = TRUE) %>%
#   row_spec(which(egger_tab$method == "Egger CHE (MLMA + RVE)"), background = "yellow") %>%
#   collapse_rows(columns = 1:2, valign = "top")



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# PET-PEESE 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#--- Multilevel model ----------------------------------------------------------
pet_mlma <- rma.mv(es_win ~ 1 + sqrt(g_var), # INTERCEPT MUST BE INCLUDED
                   V = vcov_mat,
                   random = ~ 1 | Paper / ES_ID, # inclusion of modified covariate (Egger's regression)
                   sparse = TRUE,
                   data = es_data_win, 
                   test = "t")
pet_mlma_full <- rma.mv(es_win ~ 1 + sqrt(g_var) +
                          Cognitive_domain + # control moderator
                          DV_type + # control moderator
                          Effect + # control moderator
                          Recruitment,  # control moderator, 
                        V = vcov_mat,
                        random = ~ 1 | Paper / ES_ID, # inclusion of modified covariate (Egger's regression)
                        sparse = TRUE,
                        data = es_data_win,
                        test = "t")

#--- RVE models with correlated or hierarchical weights ------------------------
## version used in Bediou et al 2018: hierarchical weights 
## https://osf.io/3gd8v/ 
pet_rve.hier <- robu(formula = es_win ~ 1 + sqrt(g_var), # INTERCEPT MUST BE INCLUDED
                     var.eff.size = g_var,
                     studynum = Paper,
                     modelweights = "HIER",
                     small = TRUE,
                     digits = 3,
                     data = es_data_win)
pet_rve.corr <- robu(formula = es_win ~ 1 + sqrt(g_var), # INTERCEPT MUST BE INCLUDED
                     var.eff.size = g_var,
                     studynum = Paper,
                     modelweights = "CORR",
                     small = TRUE,
                     digits = 3,
                     data = es_data_win)
pet_rve.hier_full <- robu(formula = es_win ~ 1 + sqrt(g_var)+
                            Cognitive_domain + # control moderator
                            DV_type + # control moderator
                            Effect + # control moderator
                            Recruitment , # control moderator, 
                          var.eff.size = g_var,
                          studynum = Paper,
                          modelweights = "HIER",
                          small = TRUE,
                          digits = 3,
                          data = es_data_win)
pet_rve.corr_full <- robu(formula = es_win ~ 1 + sqrt(g_var)+
                            Cognitive_domain + # control moderator
                            DV_type + # control moderator
                            Effect + # control moderator
                            Recruitment , # control moderator, 
                          var.eff.size = g_var,
                          studynum = Paper,
                          modelweights = "CORR",
                          small = TRUE,
                          digits = 3,
                          data = es_data_win)
peese_mlma <- rma.mv(es_win ~ 1 + g_var, # INTERCEPT MUST BE INCLUDED
                     V = vcov_mat,
                     random = ~ 1 | Paper / ES_ID, # inclusion of modified covariate (Egger's regression)
                     sparse = TRUE,
                     data = es_data_win, 
                     test = "t")
peese_mlma_full <- rma.mv(es_win ~ 1 + g_var+
                            Cognitive_domain + # control moderator
                            DV_type + # control moderator
                            Effect + # control moderator
                            Recruitment , # control moderator, 
                          V = vcov_mat,
                          random = ~ 1 | Paper / ES_ID, # inclusion of modified covariate (Egger's regression)
                          sparse = TRUE,
                          data = es_data_win, 
                          test = "t")

#--- version used in Bediou et al 2018 -----------------------------------------
peese_rve.hier <- robu(formula = es_win ~ 1 + g_var, # ERROR: INTERCEPT MUST BE REMOVED
                       var.eff.size = g_var, 
                       studynum = Paper,
                       modelweights = "HIER", 
                       small = TRUE, 
                       digits = 3,
                       data = es_data_win, 
                       test = "t")
peese_rve.corr <- robu(formula = es_win ~ 1 + g_var, # ERROR: INTERCEPT MUST BE REMOVED
                       var.eff.size = g_var, 
                       studynum = Paper,
                       modelweights = "CORR", 
                       small = TRUE, 
                       digits = 3,
                       data = es_data_win, 
                       test = "t")
peese_rve.hier_full <- robu(formula = es_win ~ 1 + g_var+
                              Cognitive_domain + # control moderator
                              DV_type + # control moderator
                              Effect + # control moderator
                              Recruitment, # control moderator, 
                            var.eff.size = g_var, 
                            studynum = Paper,
                            modelweights = "HIER", 
                            small = TRUE, 
                            digits = 3,
                            data = es_data_win)
peese_rve.corr_full <- robu(formula = es_win ~ 1 + g_var+
                              Cognitive_domain + # control moderator
                              DV_type + # control moderator
                              Effect + # control moderator
                              Recruitment,  # control moderator, 
                            var.eff.size = g_var, 
                            studynum = Paper,
                            modelweights = "CORR", 
                            small = TRUE, 
                            digits = 3,
                            data = es_data_win)

# format results
petpeese_tab <- bind_rows(cbind("method" = "PET (RVE corr)", "model" = "NULL", tidy_robu(pet_rve.corr)),
                          cbind("method" = "PET (RVE hier)", "model" = "NULL", tidy_robu(pet_rve.hier)),
                          cbind("method" = "PET (RVE corr)", "model" = "FULL", tidy_robu(pet_rve.corr_full)),
                          cbind("method" = "PET (RVE hier)", "model" = "FULL", tidy_robu(pet_rve.hier_full)),
                          cbind("method" = "PET (CHE)", "model" = "NULL", tidy_CHE(pet_mlma)),
                          cbind("method" = "PET (CHE)", "model" = "FULL", tidy_CHE(pet_mlma_full)),
                          cbind("method" = "PEESE (RVE corr)", "model" = "NULL", tidy_robu(peese_rve.corr)),
                          cbind("method" = "PEESE (RVE hier)", "model" = "NULL", tidy_robu(peese_rve.hier)),
                          cbind("method" = "PEESE (RVE corr)", "model" = "FULL", tidy_robu(peese_rve.corr_full)),
                          cbind("method" = "PEESE (RVE hier)", "model" = "FULL", tidy_robu(peese_rve.hier_full)),
                          cbind("method" = "PEESE (CHE)", "model" = "NULL", tidy_CHE(peese_mlma)),
                          cbind("method" = "PEESE (CHE)", "model" = "FULL", tidy_CHE(peese_mlma_full))
                          ) %>%
  mutate(term = str_replace_all(term, "sqrt.g_var.", "slope")) %>%
  mutate(term = str_replace_all(term, "sqrt(g_var)", "slope")) %>%
  mutate(term = str_replace_all(term, "g_var", "slope")) %>%
  filter(term == "slope" | term == "intrcpt") %>%# keep only the slope & intercept
  arrange(desc(method), term)
# # visualize table
# petpeese_tab %>% 
#   kable() %>%
#   kable_styling(bootstrap_options = c("condensed", "striped", "hover"), full_width = FALSE, position = "left") %>%
#   row_spec(which(petpeese_tab$p<.05), bold = TRUE) %>%
#   row_spec(which(petpeese_tab$method == "Egger CHE (MLMA + RVE)"), background = "yellow") %>%
#   collapse_rows(columns = 1:2, valign = "top")
# # plot results 
# petpeese_tab %>%
#   ggplot(., aes(x=method, y = b, ymin = b-se, ymax = b+se, shape = model, lty = model))+
#   geom_point(position = position_dodge(width =.2), cex= 3)+
#   geom_linerange(position = position_dodge(width =.2))+
#   geom_hline(yintercept = 0)+
#   coord_flip()+
#   scale_shape_manual(values = c(19, 1))+
#   facet_grid(~term, scales = "free_x")+
#   labs(x = "estimate", y = "")+
#   theme_apa()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Trim and Fill -- Bootstrappedpu
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
nStrap <- 1000
bootstrapped_tf <- bootstrap_tf(es_data_win, Paper)
# tidy_boot(boostrapped_tf$b)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# P-UNIFORM and 3-PSM 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#--- P-Uniform (van Assen, van Aert & Wicherts, 2015) --------------------------
bootstrapped_punif <- bootstrap_punif(es_data_win, Paper)
# tidy_boot(bootstrapped_punif$est)

#--- 3-Parameter Selection Model (3PSM, Hedges & Vevea 1996) -------------------      
# Selection models are a general class of models that attempt to model the process 
# by which the studies included in a meta-analysis may have been influenced by 
# some form of publication bias. If a particular selection model is an adequate approximation 
# for the underlying selection process, then the model provides estimates of the 
# parameters of interest (e.g., the average true outcome and the amount of heterogeneity in the true outcomes) 
# that are ‘corrected’ for this selection process (i.e., they are estimates of the parameters 
# in the population of studies before any selection has taken place). 

# first run univariate model on random sampled data (1 es per Paper)
mod_uni <- es_data_win %>% 
  group_by(Paper) %>%
  sample_n(1, replace = FALSE, set.seed = 123) %>%
  # ungroup() %>% 
  rma(es_win ~ 1, vi = g_var, dat = ., test = "t") 

# fit step model
sel <- selmodel(x = mod_uni, type="stepfun", alternative="greater", steps=c(.025,.10,.50,1))

# fit half-normal, negative-exponential, logistic, and power selection models
sel1 <- selmodel(mod_uni, type="halfnorm", alternative="two.sided")
sel2 <- selmodel(mod_uni, type="negexp",   alternative="two.sided")
sel3 <- selmodel(mod_uni, type="logistic", alternative="two.sided")
sel4 <- selmodel(mod_uni, type="power",    alternative="two.sided")

# plot selection functions
# par(mfrow=c(1,1))
# plot(sel)
# plot(sel1, add=TRUE, col="magenta")
# plot(sel2, add=TRUE, col="blue")
# plot(sel3, add=TRUE, col="red")
# plot(sel4, add=TRUE, col="green")
# legend(.2, .8, legend=c("stepfun", "halfnorm", "negexp", "logistic", "power"),
#        col=c("black", "magenta", "blue", "red", "green"), lty=1, lwd=2, cex=0.8)



#--- 3-Parameter Selection Model (3PSM, Hedges & Vevea 1996) -------------------      
# Here we focus on the 3-PSM approach, which models publication bias with 3-parameters. 
# The first one represents how much less likely a non-significant result is to be published 
# than a significant result. The other two parameters represent the estimated bias-adjusted 
# mean effect and the estimated heterogeneity of the effects.   
# run 3PSM with different step parameters
# bootstrapped_3psm <- bind_rows(bootstrap_3psm(es_data_win, Paper, .6),
#                               bootstrap_3psm(es_data_win, Paper, .2),
#                               bootstrap_3psm(es_data_win, Paper, .1),
#                               bootstrap_3psm(es_data_win, Paper, .05),
#                               # bootstrap_3psm(es_data_win, Paper, .025),
#                               bootstrap_3psm(es_data_win, Paper, .01),
#                               bootstrap_3psm(es_data_win, Paper, .001))
# save(bootstrapped_3psm, file = 'bootstrapped_3psm.RData')
# bootstrapped_3psm_cogdom <- bind_rows(bootstrap_3psm_withmoderator(es_data_win, Paper, .6),
#                                      bootstrap_3psm_withmoderator(es_data_win, Paper, .2),
#                                      bootstrap_3psm_withmoderator(es_data_win, Paper, .1),
#                                      bootstrap_3psm_withmoderator(es_data_win, Paper, .05),
#                                      # bootstrap_3psm_withmoderator(es_data_win, Paper, .025),
#                                      bootstrap_3psm_withmoderator(es_data_win, Paper, .01),
#                                      bootstrap_3psm_withmoderator(es_data_win, Paper, .001))
# save(bootstrapped_3psm_cogdom, file = 'bootstrapped_3psm_cogdom.RData')
#
# load('bootstrapped_3psm.RData')
# load('bootstrapped_3psm_cogdom.RData')
# bootstrapped_3psm_merged <- bind_rows(cbind("method" = "3-PSM", "model" = "NULL", bootstrapped_3psm),
#                                       cbind("method" = "3-PSM", "model" = "FULL", bootstrapped_3psm_cogdom))
# save(bootstrapped_3psm_merged, file = 'bootstrapped_3psm_merged.RData')
load('bootstrapped_3psm_merged.RData')

# table(bootstrapped_3psm_merged$model)
# install.packages("beepr")
# library(beepr)
# beep(sound = 1)
# tidy_boot(bootstrapped_3psm$prob)
# tidy_boot(bootstrapped_3psm_cogdom$prob)


### Publication-bias-corrected estimates obtained from regression methods   
# compare intercepts & slopes 
## CHECK P-VALUES 
## FOR SLOPES I TOOK THE PT-ONE-TAILED (FROM MR CODE) 
## FPR INTERCEPT, THIS SEEMS TO GIVE WRONG P-VALUES 
## SO I USED THE P VALUES FROM THE MODEL ! 

pb_results <- bind_rows(
  data.frame("method" = "3-PSM", "model" = "NULL",tidy_boot(bootstrapped_3psm$prob)),
  data.frame("method" = "3-PSM", "model" = "FULL",tidy_boot(bootstrapped_3psm_cogdom$prob)),
  data.frame("method" = "trim and fill", "model" = "NULL",tidy_boot(bootstrapped_tf$tf)),
  data.frame("method" = "p-uniform", "model" = "NULL",tidy_boot(bootstrapped_punif$est)),
  egger_tab,
  petpeese_tab,
  ) %>% 
  arrange(model, method, desc(term)) %>%
  # mutate(est = as.numeric(est), se = as.numeric(se), p_nomissing = as.numeric(p_nomissing)) %>%
  mutate(pSignif = symnum(p, corr = FALSE, na = FALSE, 
                           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                           symbols = c("<.001", "<.01", "<.05", ".", " "))) 
# pb_results %>%
#   arrange(b) %>%
#   kable(., digits = 3) %>%
#   kable_styling(bootstrap_options = c("condensed", "striped", "hover"), full_width = FALSE, position = "left") %>%
#   row_spec(which(pb_intercepts$model =="full"), color = "gray15", background = "#85C1E9") %>%
#   row_spec(which(pb_intercepts$model =="null"), color = "gray15", background = "#ABEBC6") %>%
#   row_spec(which(pb_intercepts$p_nomissing < .05), bold = TRUE, italic = TRUE, color = "black") %>%
#   collapse_rows(columns = 1:2, valign = "top")

save.image('MA_cross_pubBias.RData')
