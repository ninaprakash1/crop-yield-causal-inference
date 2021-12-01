
# Randomly permute the difference column for FRT,
# then compute the estimator for that permutation.
# Helper function for compute_p_value_mc
generate_diff_test_stat <- function(diff_vec, Z = NULL) {
  if (is.null(Z)) {
    Z <- rbinom(length(diff_vec),1,.5)
  }
  tau_hat <- mean((Z*2-1)*diff_vec)
  return (tau_hat)
}

# Run an FRT given a permutation function and column that represents
# difference vector to compute the estimator.
compute_p_value_mc <- function(f, diff_vec, n_draws = 1000, seed = 42) {
  set.seed(seed)
  T_sampled <- replicate(n = n_draws, expr = f(diff_vec), simplify = TRUE)
  T_obs <- f(diff_vec, Z = rep(1, length(diff_vec)))
  if (T_obs < mean(T_sampled)) {
    p_val <- (sum(T_sampled == T_obs)/2 + sum(T_sampled < T_obs))/n_draws
  } else {
    p_val <- (sum(T_sampled == T_obs)/2 + sum(T_sampled > T_obs))/n_draws
  }
  return(p_val)
}

# Run FRT on given data by stratifying on the given covariate.

runFRT_each_strata <- function(data, covariate, n_draw, seed) {
  # @desc               Run FRT using diff-in-means and log yield ratio
  #                     test statistics.
  # @param  data        Full data to stratify on
  # @param  covariate   Name of covariate variable
  # @param  n_draw      Number of iterations for FRT
  # @param  seed        Seed for randomization
  # @return             Data frame of p-values for each covariate value for each
  #                     test statistic.
  
  cov_vals = sort(unique(data[[covariate]]))
  C = length(cov_vals)
  
  P.Value.Ln_ratio = vector(length = C)
  
  for (i in 1:C) {
    # Get strata
    mask = data[[covariate]] == cov_vals[i]
    
    # Compute p-value for log ratio
    p_val_lr = compute_p_value_mc(f = generate_diff_test_stat,
                                  diff_vec = data$Ln_ratio[mask],
                                  n_draw = n_draw,
                                  seed = seed)
    P.Value.Ln_ratio[i] = p_val_lr
  }
  
  return (data.frame(cov_vals, P.Value.Ln_ratio))
  
}

# Compute p-value using finite-sample CLT.

compute_CLT_estimator <- function(data, covariate, seed) {
  # @desc                 Get estimator using finite-sample CLT (t-test)
  # @param    data        Dataset to subset from
  # @param    covariate   Name of covariate variable to stratify on
  # @param    seed        Random seed
  # @return               Vector of p_vals for each covariate value
  
  cov_vals = sort(unique(data[[covariate]]))
  C = length(cov_vals)
  
  P.Value.CLT = vector(length=C)
  
  for (i in 1:C) {
    mask = data[[covariate]] == cov_vals[i]
    results = t.test(x = data$Ln_ratio[mask], mu=0, paired=FALSE)
    P.Value.CLT[i] = results[[3]]
  }
  
  return (data.frame(cov_vals, P.Value.CLT))
  
}

# Compute rank statistic estimator
compute_rank_statistic <- function(data, covariate, n_draw=10000,seed=100) {
  
  cov_vals = sort(unique(data[[covariate]]))
  C = length(cov_vals)
  
  Tau.Rank = vector(length=C)
  P.Value.Rank = vector(length=C)

  # For each covariate value
  for (i in 1:C) {
    
    # Filter the data for that value
    data_subset = data[data[[covariate]] == cov_vals[i],]
    
    # Get the rank statistics after subsetting
    rank_all = order(c(unlist(data_subset$Yield.of.CT), unlist(data_subset$Yield.of.NT)))
    data_subset$CT.Rank = rank_all[1:(length(rank_all)/2)]
    data_subset$NT.Rank = rank_all[(length(rank_all)/2 + 1) : (length(rank_all))]
    data_subset$Rank.Diff = data_subset$NT.Rank - data_subset$CT.Rank
    
    # Compute the estimate T_rank
    Tau.Rank[i] = mean(data_subset$Rank.Diff)
    
    # Compute the p-value for T_rank
    P.Value.Rank[i] = compute_p_value_mc(f = generate_diff_test_stat,
                                         diff_vec = data_subset$Rank.Diff,
                                         n_draw = n_draw,
                                         seed = seed)
  }
  return (data.frame(cov_vals, Tau.Rank, P.Value.Rank))
}

get_basic_causal_results <- function(data, crop, covariate, n_draw=10000, seed=100) {
  # @param  crop       Crop name
  # @param  covariate  Name of covariate variable to stratify on
  # @param  n_draw     Number of iterations for FRT
  # @param  seed       Random seed for FRT
  
  # 1. Get subset of data for crop
  data_subset = subset(data, Crop==crop)
  data_subset <- data_subset[!is.na(data_subset[[covariate]]),]
  
  
  # 2. Get unique covariate values
  cov_vals = sort(unique(data_subset[[covariate]]))
  
  # 3. Get within-strata results, for each covariate value
  #         Log yield ratio and its Neymanian variance and CI
  Within.Strata.Results <- data_subset %>% group_by(.data[[covariate]]) %>%
    dplyr::summarise(
      Mean.Ln_ratio = mean(Ln_ratio),
      Var.Ln_ratio = var(Ln_ratio)/(n()-1),
      n = n()) %>%
    
    mutate(CI.L.Ratio  = Mean.Ln_ratio - qnorm(.975)*sqrt(Var.Ln_ratio), 
           CI.U.Ratio  = Mean.Ln_ratio + qnorm(.975)*sqrt(Var.Ln_ratio))

  # 4. Get p-value using FRT
  FRT_results = runFRT_each_strata(data_subset,covariate,n_draw,seed)
  Within.Strata.Results['P.Value.Ln_ratio'] = FRT_results['P.Value.Ln_ratio']
  
  # 5. Get p-value using CLT
  CLT_results = compute_CLT_estimator(data_subset,covariate,seed)
  Within.Strata.Results['P.Value.CLT'] = CLT_results['P.Value.CLT']
  
  # 6. Get rank statistic and its p-value
  Rank_results = compute_rank_statistic(data_subset,covariate,n_draw,seed)
  Within.Strata.Results['Tau.Rank'] = Rank_results['Tau.Rank']
  Within.Strata.Results['P.Value.Rank'] = Rank_results['P.Value.Rank']

  # 7. Get estimator across strata (still within crop subset)
  #         ATE, variance of ATE, and 95% CI for ATE
  n_total <- sum(Within.Strata.Results$n)
  Tau.ATE.Log.Ratio <- sum(Within.Strata.Results$Mean.Ln_ratio* Within.Strata.Results$n/n_total)
  Var.ATE.Log.Ratio <- sum(Within.Strata.Results$Var.Ln_ratio* (Within.Strata.Results$n/n_total)**2)
  
  CI.ATE.Log.Ratio <- c(Tau.ATE.Log.Ratio - qnorm(.975)*Var.ATE.Log.Ratio, Tau.ATE.Log.Ratio + qnorm(.975)*Var.ATE.Log.Ratio)
  CI.ATE.Log.Ratio.L = CI.ATE.Log.Ratio[1]
  CI.ATE.Log.Ratio.U = CI.ATE.Log.Ratio[2]
  
  Across.Strata.Results = data.frame(Tau.ATE.Log.Ratio, Var.ATE.Log.Ratio,
                                     CI.ATE.Log.Ratio.L, CI.ATE.Log.Ratio.U)
  
  # 8. Make plot of 95% CI
  p <- ggplot(Within.Strata.Results, aes_string(x=covariate, y="Mean.Ln_ratio")) +
    geom_point( size = 3, color = 'purple') +
    geom_errorbar(aes(ymax = CI.U.Ratio, ymin = CI.L.Ratio), color = 'purple', alpha = .8) +
    coord_flip() +
    xlab("")
  
  p = p + geom_hline(yintercept= 0, color= 'black', linetype = 'dashed', alpha = .8) + 
    ggtitle(paste("95% CI for Log Ratio for ", crop, " and ", covariate))  + theme(text=element_text(size=10))
  
  return (list(Within.Strata.Results,Across.Strata.Results, p))
}

# Main automation function. Runs basic analysis for each crop type, stratified
# on given covariate (column name).
get_basic_causal_results_all_crops <- function(data, covariate, n_draw, seed) {
  crop_types = sort(unique(data$Crop))
  
  # For each crop type, stratify by the given covariate and save results
  for (i in 1:length(crop_types)) {
    crop = crop_types[i]
    print(crop)
    results = get_basic_causal_results(data, crop,covariate)
    
    if (!dir.exists('./results')) {
      dir.create('./results')
    }
    
    write.csv(results[1], paste('./results/Within.Strata.Results',crop,covariate,'.csv',sep='_'))
    write.csv(results[2], paste('./results/Across.Strata.Results',crop,covariate,'.csv',sep='_'))
    ggsave(paste('./results/CI_',crop,'_',covariate,'.png')) # save the most revent ggplot (p)
  }
}

