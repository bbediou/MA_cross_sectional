library(readxl) # for reading excel files
library(tidyverse) # for data-cleaning
library(kableExtra) # for making formatted table
library(robumeta) # for original implementation of RVE
library(metafor) # for multivariate meta-analytic working models 
library(clubSandwich) # for RVE with metafor; requires version 0.5.1 or higher
library(DT) # for datatable and formatRound
library(DescTools) # for Winsorize

setwd('/Users/bediou/GoogleDrive/Meta-Analysis/DATA/ANALYSIS/2021/')

# load included data
data <- readxl::read_xlsx('/Users/bediou/Dropbox/Bavelier Lab Team Folder/MA files/ES_cross-sectional_2022_01.xlsx',
# data <- readxl::read_xlsx('/Users/bediou/Dropbox/Bavelier Lab Team Folder/MA files/ES_cross-sectional_2021-07.xlsx',
                          range = "A3:AM461", sheet = "data") %>%
  filter(!is.na(Paper)) %>%
  replace_na(list(SubSample = "")) %>%
  mutate_at(vars(nExp:MaleRatioDiff), as.numeric) %>%
  mutate_at(vars(`Cohen's d`:g_var), as.numeric) %>%
  mutate(Recruitment = case_when(  # recode problematic recruitment levels
    Recruitment == "Overt" ~ "Overt",
    Recruitment == "Covert" ~ "Covert",
    Recruitment == "Covert but explicit goals" ~ "Overt",
    Recruitment == "OVERT + COVERT" ~ "Overt",
    Recruitment == "overt?" ~ "Overt",
    Recruitment == "Overt?" ~ "Overt",
    TRUE ~"Recruitment!?"
    )
  )

# load excluded data
data_excluded <- readxl::read_xlsx('/Users/bediou/Dropbox/Bavelier Lab Team Folder/MA files/ES_cross-sectional_2022_01.xlsx',
                          # data <- readxl::read_xlsx('/Users/bediou/Dropbox/Bavelier Lab Team Folder/MA files/ES_cross-sectional_2021-07.xlsx',
                          range = "A3:AM251", sheet = "excluded") %>%
  # filter(!is.na(Paper)) %>%
  replace_na(list(SubSample = "")) %>%
  mutate_at(vars(nExp:MaleRatioDiff), as.numeric) %>%
  mutate_at(vars(`Cohen's d`:g_var), as.numeric) %>%
  mutate_at(vars(SS:SS_overlap), as.character) %>% # convert to character
  mutate(`RT/accuracy`  = as.character(`RT/accuracy`)) %>% # convert to character
  mutate(Recruitment = case_when(  # recode problematic recruitment levels
    Recruitment == "Overt" ~ "Overt",
    Recruitment == "Covert" ~ "Covert",
    Recruitment == "Covert but explicit goals" ~ "Overt",
    Recruitment == "OVERT + COVERT" ~ "Overt",
    Recruitment == "overt?" ~ "Overt",
    Recruitment == "Overt?" ~ "Overt",
    TRUE ~"Recruitment!?"
  )
  )

# define moderators to recode as factor
moderators <- c("Cognitive_domain", "DV_type", "Effect", "Recruitment")
domains <- c("perception", 
             "motor control", 
             "bottom-up attention", 
             "top-down attention", 
             "inhibition", 
             "spatial cognition", 
             "multi-tasking", 
             "verbal cognition", 
             "problem solving")



# filter data
es_data <- data %>% #filter(is.na(keep)) %>% View()
  filter(!is.na(keep) & keep != "IQ" & keep!="MISSING") %>% 
  mutate_at(vars(Paper:SS_overlap), factor) %>%# recode cluster info as factor
  mutate_at(vars(moderators), factor) %>%# moderators as factor
  mutate(Cognitive_domain = factor(Cognitive_domain, 
                                   levels = domains, 
                                   ordered = FALSE)) %>%
  mutate(Recruitment = droplevels(Recruitment)) %>% # removed unused levels
  mutate_at(vars(moderators), ~droplevels(.)) %>%# remove unused levels
  mutate(Paper = as.character(Paper)) %>%# recode paper as character (or wald_test _0.53 crashes)
  mutate(g_se = sqrt(g_var),# add significance
         z = g / sqrt(g_var),
         p.2t = 2*pnorm(z, lower.tail = F),
         p.1t = pnorm(z, lower.tail = F)) %>%
  mutate(pSignif.2t = if_else(p.2t < .05, "significant", "non-significant")) %>%
  # Create new variable for Egger Sandwich/Egger MLMA analysis as well as ES ID variable for MLMA
  mutate(Va = 4 / (nExp + nCtl), 
         sda = sqrt(Va),
         ES_ID = factor(1:nrow(.))) %>%
  filter(!is.na(g)) # There are missing values for g so I dropped them to test the code fully. 

# winsorize
es_data <- es_data %>% 
  mutate_at(moderators, droplevels) %>%
  mutate(es_win = Winsorize(es_data$g, na.rm = TRUE), 
         flg_win = ifelse(g!=es_win, 1, 0),
         g_final = es_win) 

es_data_win <- es_data
save.image(file = "MA_data_cross.RData")






