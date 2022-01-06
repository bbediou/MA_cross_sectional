# bootstrapping functions
# Note the use of enquo & !! to pass arguments to dplyr

#--- Trim and fill -------------------------------------------------------------
# this function requires independent effects so we sample 1 ES per cluster
bootstrap_tf <- function(data, cluster, nStrap = 1000) {
  out <- NULL
  cluster_var <- enquo(cluster)# cluster for evaluation
  for (i in 1:nStrap) {
    run <- data %>% 
      group_by(!!cluster_var) %>% 
      sample_n(1) %>% 
      rma(es_win ~ 1, vi = g_var, data = .)
    tf <- trimfill(run)# trim & fill
    tmp <- data.frame(run = i, 
                      b = run$beta,k = run$k,  
                      tf = tf$beta,k_added = tf$k0)
    out <- rbind(out, tmp)
  }
  return(out)
}

#--- 3-PSM ---------------------------------------------------------------------
# note the use of likelihood ratio test
bootstrap_3psm <- function(data, cluster, step = .075, nStrap = 1000) {
  out <- NULL
  cluster <- enquo(cluster)
  print(paste0("step",step))  
  for (i in 1:nStrap) {
    if(i %in% seq(100, nStrap, 100)){print(paste0("step:", step, " - sample: ", i))}
    iter <- 1
    err <- 0
    while(TRUE){# repeat until model converges
      # display counter
      if(iter %in% seq(10,100,10)){print(paste0("step:", step, " - sample: ", i, "- iteration: ",iter))}
      #run code
      psm_fit <- data %>% 
        group_by(!!cluster) %>% 
        sample_n(1) %>% 
        with(., try(weightfunct(es_win, g_var, steps = c(step, 1)),
                    silent = TRUE), verbose = TRUE)
      # check if error
      if(!is(psm_fit, 'try-error')) break# repeat until model converges
      iter <- iter + 1
      err <- err + 1
    }
    # save results
    lrchisq <- 2 * (abs(psm_fit[[1]]$value - psm_fit[[2]]$value))# likelihood test ratio
    pval <- 1 - pchisq(lrchisq, 1L)
    psm_results <- with(psm_fit, data.frame(
      sample = i, 
      step = step,
      # fit = fit,
      niter = iter,
      prob = psm_fit[[2]]$par[length(par)],
      lrchisq = lrchisq,
      p_nomissing = pval))
    
    out <- rbind(out, psm_results)
    
  }
  
  # display warning if more than 5 errors
  if(err>5) {print(paste0("TOO MANY ERRORS (",err,"), CHECK CODE!!!)"))}
  return(out)
}
# bootstrap_3psm(es_data_filt, SubSample, step = 0.2)

bootstrap_3psm_withmoderator <- function(data, cluster, step = .075, nStrap = 1000) {
  out <- NULL
  cluster <- enquo(cluster)
  print(paste0("step: ",step))
  for (i in 1:nStrap) {
    if(i %in% seq(100, nStrap, 100)){print(paste0("step:", step, " - sample: ", i))}
    iter <- 1
    err <- 0
    while(TRUE){# repeat until model converges
      if(iter %in% seq(10, 100,10)){print(paste0("step:", step, " - sample: ", i, "- iteration: ",iter))}
      psm_fit <- data %>%
        group_by(!!cluster) %>%
        sample_n(1) %>%
        with(., try(weightfunct(es_win, g_var,# repeat until model converges
                                mods = ~ Cognitive_domain,
                                steps = c(step, 1)), silent = TRUE), verbose = TRUE)
      if(!is(psm_fit, 'try-error')) break# repeat until model converges
      iter <- iter + 1
      # err <- err + 1
      if (iter>100) {# stop loop if too many iterations
        break
        print(paste0("TOO MANY ERRORS (",err,"), CHECK CODE!!!)"))
        err <- err + 1
      }
    }
    
    lrchisq <- 2 * (abs(psm_fit[[1]]$value - psm_fit[[2]]$value))# likelihood test ratio
    pval <- 1 - pchisq(lrchisq, 1L)
    psm_results <- with(psm_fit, data.frame(
      sample = i, 
      step = step,
      # fit = fit,
      niter = iter,
      prob = psm_fit[[2]]$par[length(par)],
      lrchisq = lrchisq,
      p_nomissing = pval))
    out <- rbind(out, psm_results)
  }
  
  # if (iter>10000) {
  #   break
  #   print(paste0("TOO MANY ERRORS (",err,"), CHECK CODE!!!)"))
  #   err <- err + 1
  #   }
  
  if (err>5) {print(paste0("TOO MANY ERRORS (",err,"), CHECK CODE!!!)"))}
  return(out)
}

bootstrap_punif <- function(data, cluster, nStrap = 1000) {
  out <- NULL
  cluster <- enquo(cluster)
  for (i in 1:nStrap) {
    run <- data %>% 
      group_by(!!cluster) %>%  
      sample_n(1) %>% 
      with(., puniform(yi = es_win, vi = g_var, side = "right", method = "ML", plot = FALSE))
    tmp <- with(run, data.frame(est, ci.lb, ci.ub, L.0, pval.0, ksig, L.pb, pval.pb))
    out <- rbind(out, tmp)
  }
  # Export result after nStrap cycles
  return(out)
}

