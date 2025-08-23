setwd('...')

source("funtions.R")
devtools::source_url("https://raw.githubusercontent.com/Jeremy690/False-Discovery-Rate-via-Data-Splitting/main/code/functions/linear/DS.R")
devtools::source_url("https://raw.githubusercontent.com/Jeremy690/False-Discovery-Rate-via-Data-Splitting/main/code/functions/linear/analysis.R")
devtools::source_url("https://raw.githubusercontent.com/Jeremy690/False-Discovery-Rate-via-Data-Splitting/main/code/functions/linear/fdp_power.R")

library(glmnet)
library(CVXR)
library(copula)
library(stabledist)
library(glasso)
library(corpcor)
library(scalreg)
library(dbscan)

N = 500
p = 1000
r_vector = seq(0, 0.2, 0.05)
lasso_avgs = sub_avgs = lasso_avgs_len = sub_avgs_len = sub_test = lasso_test = matrix(0, length(r_vector), 2)
IFV = list()

for(u in 1:length(r_vector)){
  k = r_vector[u]
  aaa = list()
  lasso_avgs1 = sub_avgs1 = lasso_avgs_len1 = sub_avgs_len1 = sub_test1 = lasso_test1 = matrix(0, N, 2)
  
  for(m in 1:N){
    set.seed(m)
    t1 = Sys.time()
    s <- 12
    n <- 300
    alpha <- 1
    sigma <- 1
    rr = 0.3
    scenario = 2
    
    if(scenario == 1){
      distr = 'uniform'
      gdata = generate_target_data(m, rr, n, p, s, sigma, k, distr)
    } else if(scenario == 2){
      gdata = generate_target_data2(m, rr, n, p, s, sigma, k, distr)
    } else if(scenario == 3){
      distr = 'cauchy'
      gdata = generate_target_data0(m, rr, n, p, s, sigma, k, distr)
    }
    
    X0 = Phi = gdata$X
    y0 = gdata$y
    indice = sample(n, floor(n))
    X_train = X0[indice,]
    y_train = y0[indice]
    X_test = X0[-indice,]
    y_test = y0[-indice]
    signal_true = gdata$signal_true
    nonzero_indices = gdata$non_zero_indices
    err = abs(gdata$corruption)
    S = which(abs(signal_true) > 0)
    Sc = which(abs(signal_true) == 0)
    alpha = 0.05
    
    result_9 = sample_fun(X0, y0)
    subsample_lasso00 = result_9$subsample_lasso00
    val_set = result_9$val_set
    B = 5
    beta = subsample_lasso00
    X01 = X0[-val_set,]
    y01 = y0[-val_set]
    X_val = X0[val_set,]
    y_val = y0[val_set]
    M = matrix(0, B, p)
    qtable = c()
    sr = seq(0.7, 0.85, 0.05)
    
    for(b in 1:B){
      if(b <= 4){
        set.seed(1)
        sr_0 = sr[b]
        subsample_size = floor(n * sr[b])
      } else {
        set.seed(b)
        pos = which(qtable[,2] == qtable[,2][which.min(qtable[,2])])
        sr_0 = qtable[,1][max(pos)]
        subsample_size = floor(n * sr_0)
      }
      
      gradient_l1_norms <- gradient_norm(X01, y01, beta)
      gradient_l1_norms[gradient_l1_norms == 0] <- 1e-10
      what <- exp_normalize(gradient_l1_norms, 1e-04)
      what = what + min(what[!what == 0])
      n0 = length(y01)
      set2 <- sample(n0, subsample_size, replace = FALSE, prob = what)
      X_set5 = X01[set2,]
      y_set5 = y01[set2]
      
      cv.lasso <- cv.glmnet(X_set5, y_set5, alpha = 1)
      fit_w <- glmnet(X_set5, y_set5, alpha = 1, lambda = cv.lasso$lambda.min)
      beta_upate <- as.numeric(coef(fit_w)[-1])
      
      what <- (1 / gradient_l1_norms) / max(1 / gradient_l1_norms)
      set2 <- sample(n0, subsample_size, replace = FALSE, prob = what)
      X_set5 = X01[set2,]
      y_set5 = y01[set2]
      cv.lasso <- cv.glmnet(X_set5, y_set5, alpha = 1)
      fit_w <- glmnet(X_set5, y_set5, alpha = 1, lambda = cv.lasso$lambda.min)
      beta_upate1 <- as.numeric(coef(fit_w)[-1])
      
      esset = cbind(beta, beta_upate, beta_upate1)
      w_init = rep(0, dim(esset)[2])
      
      result <- optim(par = w_init, fn = optfun, method = "BFGS", 
                      y_test = y_val, X_test = X_val, B = esset)
      
      w_optimal = exp(result$par)
      w_optimal = w_optimal / sum(w_optimal)
      beta = esset %*% w_optimal
      
      if(b <= 4){
        qtable = rbind(qtable, c(sr[b], mean((y_val - X_val %*% beta_upate1)^2)))
      }
      
      M[b,] = beta
    }
    
    signal_index = which(!signal_true == 0)
    set2 <- sample(n0, floor(n0 * 0.7), replace = FALSE, prob = what)
    X_set5 = X01[set2,]
    y_set5 = y01[set2]
    num_split = 2
    q = 0.1
    result = DS(X_set5, y_set5, num_split, q)
    sub_test1[m,] = c(result$MDS_power, result$MDS_fdp)
    
    subsample_lasso = beta
    fit <- scalreg(X_set5, y_set5)
    noise_level <- fit$hsigma
    cov0 <- cov.shrink(X_set5)
    rho <- 0.1
    glasso_fit <- glasso(cov0, rho = rho)
    Theta_glasso <- glasso_fit$wi
    beta0 = beta_upate1
    db_lasso1 = beta0 + (1/length(y_set5)) * Theta_glasso %*% (t(X_set5)) %*% (y_set5 - X_set5 %*% beta0)
    
    CI = CI_fun(db_lasso1, Theta_glasso, noise_level, cov0, length(y_set5))
    L2 = CI$L
    U2 = CI$U
    se0 = CI$se
    
    beta0 = beta_upate
    db_lasso2 = beta0 + (1/length(y_set5)) * Theta_glasso %*% (t(X_set5)) %*% (y_set5 - X_set5 %*% beta0)
    
    CI = CI_fun(db_lasso2, Theta_glasso, noise_level, cov0, length(y_set5))
    L1 = CI$L
    U1 = CI$U
    se1 = CI$se
    
    esset = cbind(db_lasso1, db_lasso2)
    w_init = rep(0, dim(esset)[2])
    
    result <- optim(par = w_init, fn = optfun, method = "BFGS", 
                    y_test = y_val, X_test = X_val, B = esset)
    
    w_optimal = exp(result$par)
    w_optimal = w_optimal / sum(w_optimal)
    db_lasso = esset %*% w_optimal
    se = sqrt(w_optimal[1]^2 * se0^2 + w_optimal[2]^2 * se1^2 + 2 * w_optimal[1] * w_optimal[2] * se0 * se1)
    CI = CI_fun2(se, db_lasso)
    L = CI$L
    U = CI$U
    
    avgs = sum((L[S] <= signal_true[S]) & (U[S] >= signal_true[S])) / length(S)
    avgsc = sum((L[Sc] <= signal_true[Sc]) & (U[Sc] >= signal_true[Sc])) / length(Sc)
    sub_avgs_len1[m,] = c(mean(U[S] - L[S]), mean(U[Sc] - L[Sc]))
    sub_avgs1[m,] = c(avgs, avgsc)
    
    num_split = 2
    q = 0.1
    result = DS(X0, y0, num_split, q)
    lasso_test1[m,] = c(result$MDS_power, result$MDS_fdp)
    
    cv.lasso <- cv.glmnet(X0, y0, alpha = 1)
    fit_w <- glmnet(X0, y0, alpha = 1, lambda = cv.lasso$lambda.min)
    lasso <- as.numeric(coef(fit_w)[-1])
    
    rho <- 0.1
    glasso_fit <- glasso(cov0, rho = rho)
    Theta_glasso <- glasso_fit$wi
    db_lasso = lasso + (1/length(y0)) * Theta_glasso %*% (t(X0)) %*% (y0 - X0 %*% lasso)
    
    CI = CI_fun(db_lasso, Theta_glasso, noise_level, cov0, length(y0))
    L = CI$L
    U = CI$U
    avgs = sum((L[S] <= signal_true[S]) & (U[S] >= signal_true[S])) / length(S)
    avgsc = sum((L[Sc] <= signal_true[Sc]) & (U[Sc] >= signal_true[Sc])) / length(Sc)
    lasso_avgs1[m,] = c(avgs, avgsc)
    lasso_avgs_len1[m,] = c(mean(U[S] - L[S]), mean(U[Sc] - L[Sc]))
    
    t2 = Sys.time()
    print(c(u, m, (t2 - t1)))
  }
  
  sub_test[u,] = apply(sub_test1, 2, mean)
  lasso_test[u,] = apply(lasso_test1, 2, mean)
  lasso_avgs[u,] = apply(lasso_avgs1, 2, mean)
  lasso_avgs_len[u,] = apply(lasso_avgs_len1, 2, mean)
  sub_avgs[u,] = apply(sub_avgs1, 2, mean)
  sub_avgs_len[u,] = apply(sub_avgs_len1, 2, mean)
  
  ab = list(sub_test = sub_test, lasso_test = lasso_test, lasso_avgs = lasso_avgs,
            sub_avgs = sub_avgs, lasso_avgs_len = lasso_avgs_len, sub_avgs_len = sub_avgs_len)
}