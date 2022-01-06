library(tidyverse) # for data-cleaning
library(kableExtra) # for making formatted table
library(robumeta) # for original implementation of RVE
library(metafor) # for multivariate meta-analytic working models 
library(clubSandwich) # for RVE with metafor; requires version 0.5.1 or higher
library(Hmisc)# for %nin%

load("MA_data_cross.RData")

# default correlation
rho = 0.8

# list of moderators
moderators <- c("Cognitive_domain", "DV_type", "Effect", "Recruitment")

# get relative frequencies of each moderator level
mod_freq <- list()
for (i in 1:length(moderators)){
  mod_name <- moderators[i]
  mod_levels <- as.character(unique(es_data_win[[mod_name]]))
  for (j in 1:length(mod_levels)){
    level_name <- mod_levels[j]
    mod_freq[[level_name]] <- sum(es_data_win[[mod_name]] == level_name)/nrow(es_data_win)
  }
}

# moderator analysis
seq <- list(c(1,2,3,4),
            c(2,1,3,4), 
            c(3,1,2,4), 
            c(4,1,2,3))

Moderator_estimates <- data.frame()
Wald_results <- data.frame()
# i=1
for (i in 1:length(moderators)){
  # loop across sequences of primary/secondary moderators
  modseq <- seq[[i]]
  
  # moderators
  mod1 <- moderators[[modseq[1]]]
  mod2 <- moderators[[modseq[2]]]
  mod3 <- moderators[[modseq[3]]]
  mod4 <- moderators[[modseq[4]]]
  mods <- c(mod1, mod2, mod3, mod4)
  
  # impute the variance-covariance matrix
  vcov_mat <- impute_covariance_matrix(es_data_win$g_var, 
                                       cluster = es_data$Paper,                            
                                       r = rho, # 0.8 is default in RVE
                                       smooth_vi = TRUE) # average Vi within each cluster

  # run model
  model_formula = as.formula(paste("es_win ~ 0", paste(mods, collapse=" + "), sep=" + "))
  mod_mlma_nointercept <- metafor::rma.mv(yi = model_formula,
                                          V = vcov_mat,
                                          random = ~ 1 | Paper / ES_ID,# clustering
                                          sparse = TRUE,
                                          # intercept = TRUE,# ignored when using formulas
                                          data = es_data_win,
                                          test="t")
  # Wald test of moderator effect
  primary_moderator_levels <- levels(es_data_win[[mod1]])
  nlevels_primary <- length(primary_moderator_levels)
  eff_wald <- Wald_test(mod_mlma_nointercept,
                        cluster = es_data$Paper,
                        constraints = constrain_equal(1:nlevels_primary), 
                        vcov="CR2") %>% 
    mutate(moderator = mod1)
  # Append
  Wald_results <- rbind(Wald_results, eff_wald) 
  
  # Add RVE CIs and df
  mod_CHE_nointercept <- clubSandwich::coef_test(mod_mlma_nointercept, vcov = "CR2")
  # extract the coefficients from the model
  coef_est <- mod_CHE_nointercept$beta      
  mod_coef <- mod_mlma_nointercept$beta
  
  # get info about primary/secondary moderators and frequencies
  secondary_moderators <- rownames(mod_coef) %>%
    str_replace_all(., "Cognitive_domain", "") %>%
    str_replace_all(., "DV_type", "") %>%
    str_replace_all(., "Effect", "") %>%
    str_replace_all(., "Recruitment", "") 
  # get secondary moderators control levels
  secondary_levels <- secondary_moderators[secondary_moderators %nin% primary_moderator_levels]
  # extract corresponding frequencies
  tmp_names <- names(mod_freq[secondary_levels])
  tmp_freqs <- unlist(mod_freq[secondary_levels]) %>% as.vector()
  
  # loop across primary levels
  for (j in 1:nlevels_primary){
    # vector to test a given level of the primary moderator     
    test_vec <- t(c(rep(0,nlevels_primary),# moderator of interest has nlevels_primary levels!
                    tmp_freqs))# control for secondary moderators frequencies
    
    # add 1 for level to be tested, and zero for other levels of same moderator 
    test_vec[j] <- 1# add 1 to level to test
    g_hat <- sum(as.vector(test_vec)*coef_est)# correct estimate
    
    # test level vs. zero 
    Test_result <- Wald_test(mod_mlma_nointercept, 
                             constraints = test_vec, 
                             cluster = es_data_win$Paper,
                             vcov="CR2")
    df <- Test_result$df_denom# CHECK IF CORRECT? 
    Fvalue <- Test_result$F
    SE_Yhat <- g_hat/sqrt(Fvalue)
    pval <- Test_result$p
    tcrit <- qt(.975, df)
    LB <- round((g_hat - tcrit*SE_Yhat),digits=3)
    UB <- round((g_hat+ tcrit*SE_Yhat),digits=3)
    CI <- paste(LB,UB,sep=",")
    
    # add relevant info
    moderator <- mod1
    level <- levels(es_data_win[[mod1]])[j]
    k <- sum(es_data_win[[mod1]]==level)
    m <- length(unique(es_data_win$Paper[es_data_win[[mod1]]==level]))
    
    # save to dataframe     
    tmp <- data.frame("moderator" = moderator,"level"=level,"k"=k,"m"=m,
                      "g"=g_hat,"CI"=CI,"df"=df,"pval"=pval,"LB"=LB,"UB"=UB)
    
    # Append
    Moderator_estimates <- rbind(Moderator_estimates, tmp)
    
  }
}

Wald_results %>%
  mutate(moderator = str_replace_all(moderator, "\\_", " ")) %>%
  kable() %>% kable_paper() 

Moderator_estimates %>%
  mutate(moderator = str_replace_all(moderator, "\\_", " ")) %>%
  kable() %>% kable_paper() %>%
  collapse_rows(columns = 1, col_names = "moderator")

save(list = c("Wald_results", "Moderator_estimates"), file= "MA_cross_moderator-effects.RData")
