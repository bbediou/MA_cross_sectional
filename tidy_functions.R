# tidy model outputs to compare results

#--- Robumeta output -----------------------------------------------------------
# .$regtable 
# colnames: labels    b.r    SE     t  dfs     prob   CI.L  CI.U sig
tidy_robu <- function(x){
  out <- x$reg_table %>% 
    select(term = labels, b = b.r, se = SE, t, dfs, p = prob,
           ci.lb = CI.L, ci.ub = CI.U) %>%
    mutate(term = str_replace_all(term, "X.Intercept.", "intrcpt"))
  return(out)
}
# tidy_robu(egg_sand)

#--- RMA model -----------------------------------------------------------------
# colnames: b beta se zval pval ci.lb ci.ub vb........   
tidy_rma <- function(x) {
  out <- with(x, data.frame(b = beta,se, t = zval, dfs, p = pval,
                            ci.lb, ci.ub))
  out <- data.frame(term = row.names(out), out) 
  row.names(out) <- NULL
  return(out)
}

#--- CHE model -----------------------------------------------------------------
# colnames: beta SE tstat df p_Satt
# NOTE: STATS BASED ON T-TEST INSTEAD OF Z-TEST!
tidy_CHE <- function(x){
  out <- full_join(clubSandwich::coef_test(x, vcov = "CR2"),
                   clubSandwich::conf_int(x, vcov = "CR2"))
  out <- data.frame(term = rownames(x$beta), out) %>%
    rename(b = beta, se = SE, t = tstat, dfs = df, p = p_Satt, ci.lb = CI_L, ci.ub = CI_U)
  row.names(out) <- NULL
  return(out)
}
# conf_int(egg_mlma, vcov = "CR2")

#--- 3-PSM ---------------------------------------------------------------------
tidy_psm <- function(x) {
  b <- as.numeric(x$adj_est)
  se.psm <- as.numeric(x$adj_se)
  mods <- colnames(x$XX)[-1]# moderator names
  steps <- NULL
  for (i in 2:length(x$steps)) {# extract selection steps
    steps <- c(steps, paste0(x$steps[i-1], "<p<", x$steps[i]))
  }
  term <- c("variance", "intrcpt", mods, steps)  # make term column
  if(x$fe == T) term <- term[-1]# if it's a fixed-effects model, forget "variance"
  out <- data.frame(term, b = b, se = se.psm) %>% 
    mutate(z = b/se,
           p = 2*pnorm(abs(z), lower.tail = F),
           ci.lb = b - 1.96*se,
           ci.ub = b + 1.96*se)
  return(out)
}

#--- p-uniform -----------------------------------------------------------------
tidy_punif <- function(x) {
  out <- with(x, data.frame(term = "intrcpt",
                            b = est.fe,
                            se = se.fe,
                            z = zval.fe,
                            ci.lb = ci.lb.fe,
                            ci.ub = ci.ub.fe))
  out <- mutate(out, p = 2*pnorm(abs(z), lower.tail = F))
  return(out)
}

#--- bootstrap -----------------------------------------------------------------
tidy_boot <- function(x) {
  out <- data.frame(term = "intrcpt", b = mean(x), se = sd(x), 
                    ci.lb = quantile(x, .025), ci.ub = quantile(x, .975))
  row.names(out) <- NULL
  return(out)
}


