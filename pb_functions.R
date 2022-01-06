egg_sand <- function(es_data_win){
  robu(formula = es_win ~ 1 + sda, # Modified covariate for Egger's Regression test. 
       var.eff.size = g_var, # ES variance
       studynum = Paper,# clustering
       modelweights = "CORR",# Correlated weights
       # small = TRUE,# small sample correction ON => check df's if <4 then use smaller alpha to test significance. 
       data = es_data_win) 
}
egg_sand_mod <- function(es_data_win){
  robu(formula = es_win ~ 1 + sda +# Modified covariate for Egger's Regression test. 
         Cognitive_domain + # PRIMARY moderator
         DV_type + # control moderator
         Effect + # control moderator
         Recruitment,  # control moderator, , # Modified covariate for Egger's Regression test. 
       var.eff.size = g_var, # ES variance
       studynum = Paper,# clustering
       modelweights = "CORR",# Correlated weights
       # small = TRUE,# small sample correction ON => check df's if <4 then use smaller alpha to test significance. 
       data = es_data_win) 
}
egg_mlma <- function(es_data_win){
  rma.mv(es_win ~ 1  + sda, # inclusion of modified covariate (Egger's regression)
         V = vcov_mat, # variance covariance matrix defined above
         random = ~ 1 | Paper / ES_ID, # random factors (nesting)
         sparse = TRUE, # 
         # test = "t", # small sample correction isn't needed if doing CHE model with RVE
         data = es_data_win)
}
egg_mlma_mod <- function(es_data_win){
  rma.mv(es_win ~ 1  + sda +# Modified covariate for Egger's Regression test. 
           Cognitive_domain + # PRIMARY moderator
           DV_type + # control moderator
           Effect + # control moderator
           Recruitment,  # control moderator, 
         V = vcov_mat, # variance covariance matrix defined above
         random = ~ 1 | Paper / ES_ID, # random factors (nesting)
         sparse = TRUE, # 
         # test = "t", # small sample correction isn't needed if doing CHE model with RVE
         data = es_data_win)
}
pet_mlma <- function(es_data_win){
  rma.mv(es_win ~ 1 + sqrt(g_var), # INTERCEPT MUST BE INCLUDED
         V = vcov_mat,
         random = ~ 1 | Paper / ES_ID, # inclusion of modified covariate (Egger's regression)
         sparse = TRUE,
         data = es_data_win)
}
pet_mlma_mod <- function(es_data_win){
  rma.mv(es_win ~ 1 + sqrt(g_var) +
           Cognitive_domain + # control moderator
           DV_type + # control moderator
           Effect + # control moderator
           Recruitment,  # control moderator, 
         V = vcov_mat,
         random = ~ 1 | Paper / ES_ID, # inclusion of modified covariate (Egger's regression)
         sparse = TRUE,
         data = es_data_win)
}
# version used in Bediou et al 2018
pet_rve_hier <- function(es_data_win){
  robu(formula = es_win ~ 1 + sqrt(g_var), 
       var.eff.size = g_var,
       studynum = Paper,
       modelweights = "HIER",
       small = TRUE,
       digits = 3,
       data = es_data_win)
}
pet_rve_hier_mod <- function(es_data_win){
  robu(formula = es_win ~ 1 + sqrt(g_var)+
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
}
pet_rve_corr <- function(es_data_win){
  robu(formula = es_win ~ 1 + sqrt(g_var), 
       var.eff.size = g_var,
       studynum = Paper,
       modelweights = "CORR",
       small = TRUE,
       digits = 3,
       data = es_data_win)
}
pet_rve_corr_mod <- function(es_data_win){
  robu(formula = es_win ~ 1 + sqrt(g_var)+
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
}
peese_mlma <- function(es_data_win){
  rma.mv(es_win ~ 1 + g_var, 
         V = vcov_mat,
         random = ~ 1 | Paper / ES_ID, 
         sparse = TRUE,
         data = es_data_win)
}
peese_mlma_mod <- function(es_data_win){
  rma.mv(es_win ~ 1 + g_var + #
           Cognitive_domain + # control moderator
           DV_type + # control moderator
           Effect + # control moderator
           Recruitment , # control moderator, 
         V = vcov_mat,
         random = ~ 1 | Paper / ES_ID,
         sparse = TRUE,
         data = es_data_win)
}
# version used in Bediou et al 2018
peese_rve_hier <- function(es_data_win){
  robu(formula = es_win ~ 1 + g_var, 
       var.eff.size = g_var, 
       studynum = Paper,
       modelweights = "HIER", 
       small = TRUE, 
       digits = 3,
       data = es_data_win)
}
peese_rve_hier_mod <- function(es_data_win){
  robu(formula = es_win ~ 1 + g_var +
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
}
peese_rve_corr <- function(es_data_win){
  robu(formula = es_win ~ 1 + g_var, 
       var.eff.size = g_var, 
       studynum = Paper,
       modelweights = "CORR", 
       small = TRUE, 
       digits = 3,
       data = es_data_win)
}
peese_rve_corr_mod <- function(es_data_win){
  robu(formula = es_win ~ 1 + g_var+
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
}

pb_results <- tibble(data = list(es_data_win)) %>%
  mutate(
    # Egger's Sandwich tests
    # NULL MODELS
    egg_sand_null = map(data, egg_sand),
    egg_mlma_null = map(data, egg_mlma),
    # egg_CHE_null = coef_test(egg_mlma_null, vcov = "CR2"),
    # FULL MODELS
    egg_sand_full = map(data, egg_sand_mod),
    egg_mlma_full = map(data, egg_mlma_mod),
    # egg_CHE_full = coef_test(egg_mlma_mod, vcov = "CR2"),
    # PET-PEESE - CHE version
    pet_mlma_null = map(data, pet_mlma),
    # pet_CHE_null = coef_test(pet_mlma, vcov = "CR2"),
    pet_mlma_full = map(data, pet_mlma_mod),
    # pet_CHE_full = coef_test(pet_mlma_full, vcov = "CR2"),
    peese_mlma_null = map(data, peese_mlma),
    # peese_CHE_null = coef_test(peese_mlma, vcov = "CR2"),
    peese_mlma_full = map(data, peese_mlma_mod),
    # peese_CHE_full = coef_test(peese_mlma_full, vcov = "CR2"),
    # PET-PEESE - RVE versions
    pet_rve.hier_null = map(data, pet_rve_hier),
    pet_rve.hier_full = map(data, pet_rve_hier_mod),
    pet_rve.corr_null = map(data, pet_rve_corr),
    pet_rve.corr_full = map(data, pet_rve_corr_mod),
    peese_rve.hier = map(data, peese_rve_hier),
    peese_rve.hier_full = map(data, peese_rve_hier_mod),
    peese_rve.corr = map(data, peese_rve_corr),
    peese_rve.corr_full = map(data, peese_rve_corr_mod)
  )

