msefun <- function(x, y) {
  sum((x - y)^2)
}
lossfun = function(x,y,b){
  mean((y - x%*%b)^2)
}

gradient_norm <- function(X, y, beta) {
  
  y_pred <- X %*% beta  
  residuals <- y - y_pred  
  
  residuals <- as.vector(residuals) 
  
  grad <- -2 * X * residuals  
  
  grad_ss <- rowSums(grad^2)  
  return(grad_ss)
}
exp_normalize <- function(x, lambda = 1) {
  x_shifted <- x - min(x)
  
  exp_x <- exp(-lambda * x)
  
  (exp_x - min(exp_x)) / (max(exp_x) - min(exp_x))
  exp_x
}

###############################
#FBOD
###############################
Func.FBOD <- function(data, iter, k.nn) {
  d <- dim(data)[2]
  res <- matrix(NA, nrow = nrow(data), ncol = iter)
  for (i in 1:iter) {
    l <- sample(x = round(d/2, digits = 0):(d-1), size = 1)
    ind <- sample(x = 1:d, size = l, replace = FALSE)
    data.use <- data[, ind, drop = FALSE]
    score <- dbscan::lof(data.use, minPts = k.nn + 1) 
    res[, i] <- score
  }
  res[is.nan(res)] <- NA
  res[is.infinite(res)] <- max(res[!is.infinite(res)], na.rm = TRUE)
  res.final <- rowMeans(res, na.rm = TRUE)
  return(  res.final )
}

###############################
#pilot estimation
###############################
sample_fun = function(X0,y0){
  val_s = c()
  
  i = 1
  while(i <= 5){
    
    set.seed(i)
    X_train = X0
    y_train = y0
    n0 = length(y_train)
    set3= sample(n0,floor(n0/2))
    X_set3 = X_train[set3,]
    y_set3 = y_train[set3]
    
    cv.lasso <- cv.glmnet( X_set3, y_set3, alpha=1) 
    fit_w <- glmnet(X_set3, y_set3  , alpha = 1,lambda = cv.lasso$lambda.min)
    subsample_lasso0 <- as.numeric(coef(fit_w)[-1])
    ser(signal_true,  subsample_lasso0)
    gradient_l1_norms <- gradient_norm(X0, y0, subsample_lasso0)
    dim(X_train)
    set3= sample(n0,floor(n0/2), replace = FALSE,1/gradient_l1_norms)
    X_set3 = X_train[set3,]
    y_set3 = y_train[set3]
    
    cv.lasso <- cv.glmnet( X_set3, y_set3, alpha=1) 
    fit_w <- glmnet(X_set3, y_set3  , alpha = 1,lambda = cv.lasso$lambda.min)
    subsample_lasso_0 <- as.numeric(coef(fit_w)[-1])
    gradient_l1_norms <- gradient_norm(X0, y0, subsample_lasso_0)
    set4 =    order(gradient_l1_norms)[floor(n0*0.2):floor(n0*0.8)]
    val_set =   sample(set4,30)
    
    X_val = X_train[val_set,]
    y_val = y_train[val_set]
    set5= setdiff(set4,val_set)
    X_set5 = X_train[set5,]
    y_set5 = y_train[set5]
    cv.lasso <- cv.glmnet( X_set5, y_set5, alpha=1)
    fit_w <- glmnet(X_set5, y_set5  , alpha = 1,lambda = cv.lasso$lambda.min)
    subsample_lasso00 <- as.numeric(coef(fit_w)[-1])
    val = mean((y_val-X_val%*%subsample_lasso00)^2)
    val_s = c(val_s,val)
    i = i +1 
  }
  seed = which.min(val_s)
  set.seed(seed)
  X_train = X0
  y_train = y0
  n0 = length(y_train)
  set3= sample(n0,floor(n0/2))
  X_set3 = X_train[set3,]
  y_set3 = y_train[set3]
  
  cv.lasso <- cv.glmnet( X_set3, y_set3, alpha=1) 
  fit_w <- glmnet(X_set3, y_set3  , alpha = 1,lambda = cv.lasso$lambda.min)
  subsample_lasso0 <- as.numeric(coef(fit_w)[-1])
  ser(signal_true,  subsample_lasso0)
  gradient_l1_norms <- gradient_norm(X0, y0, subsample_lasso0)
  dim(X_train)
  set3= sample(n0,floor(n0/2), replace = FALSE,1/gradient_l1_norms)
  X_set3 = X_train[set3,]
  y_set3 = y_train[set3]
  
  cv.lasso <- cv.glmnet( X_set3, y_set3, alpha=1) 
  fit_w <- glmnet(X_set3, y_set3  , alpha = 1,lambda = cv.lasso$lambda.min)
  subsample_lasso_0 <- as.numeric(coef(fit_w)[-1])
  gradient_l1_norms <- gradient_norm(X0, y0, subsample_lasso_0)
  set4 =    order(gradient_l1_norms)[floor(n0*0.2):floor(n0*0.8)]
  val_set =   sample(set4,30)
  X_val = X_train[val_set,]
  y_val = y_train[val_set]
  set5= setdiff(set4,val_set)
  X_set5 = X_train[set5,]
  y_set5 = y_train[set5]
  cv.lasso <- cv.glmnet( X_set5, y_set5, alpha=1)
  fit_w <- glmnet(X_set5, y_set5  , alpha = 1,lambda = cv.lasso$lambda.min)
  subsample_lasso00 <- as.numeric(coef(fit_w)[-1])
  abcd = list(
    subsample_lasso00 = subsample_lasso00,
    val_set = val_set
  )
  return(abcd)
}

###############################
#Signal-to-error Ratio (SER)
###############################

ser <- function(signal_true, signal_estimated) {
  error <- sum((signal_true - signal_estimated)^2)
  power <- sum(signal_true^2)
  if (power != 0 && error != 0) {
    return(10 * log10(power / error))
  } else {
    return(Inf)
  }
}

###############################
#generate target dataset
###############################
generate_target_data0 <- function(i,rr,n,p,s,sigma_eps,k,distr){
  set.seed(i)
  sigma <- diag(p) 
  mu <- rep(0, p)
  for (i in 1:(p-1)) {
    sigma[i, i+1] <- rr 
    sigma[i+1, i] <- rr 
  }
  Phi <- mvrnorm(n , mu = mu, Sigma = sigma) 
  signal_true <- rep(0, p)
  nonzero_indices0 <- sample(p, s)
  signal_true[nonzero_indices0] <- sign(rnorm(s))
  measurements_clean <- Phi %*% signal_true
  zero_indices <- sample(n, size = floor((1-k/2)*n), replace = FALSE)
  non_zero_indices <- setdiff(seq_len(n), zero_indices)
  nc = length(non_zero_indices)
  if(distr == 'cauchy'){
    er = rcauchy(nc * p, location = 0.5, scale = 0.5)
  }else if(distr == 'uniform'){
    er = runif(nc * p, min = 0, max = 4)
  }
  Phi[non_zero_indices,] =  matrix(er, nrow = nc, ncol = p)
  zero_indices <- sample(n, size = floor((1-k/2)*n), replace = FALSE)
  non_zero_indices1 <- setdiff(seq_len(n), zero_indices)
  noise <- rnorm(n,0,sigma_eps)
  if(distr == 'cauchy'){
    ery = abs(rstable(length(non_zero_indices), alpha = 1.5, beta = 0, gamma = 3, delta = 0))
  }else if(distr == 'uniform'){
    ery = runif(length(non_zero_indices), min = 3, max = 10)
  }
  
  vector <- rep(0, n)
  vector[non_zero_indices1] <-   ery 
  err = rep(0,n)
  err[c(non_zero_indices1,non_zero_indices)] = 1
  return(list(X = Phi, 
              y =  measurements_clean + noise + vector,
              Sigma = sigma,
              signal_true = signal_true,
              non_zero_indices = nonzero_indices0,
              corruption =  err))
}
generate_target_data <- function(i,rr,n,p,s,sigma_eps,k,distr){
  set.seed(i)
  sigma <- diag(p) 
  mu <- rep(0, p)
  for (i in 1:(p-1)) {
    sigma[i, i+1] <- rr 
    sigma[i+1, i] <- rr 
  }
  Phi <- mvrnorm(n , mu = mu, Sigma = sigma) 
  signal_true <- rep(0, p)
  nonzero_indices0 <- sample(p, s)
  signal_true[nonzero_indices0] <- sign(rnorm(s))
  measurements_clean <- Phi %*% signal_true
  vector <- rep(NA, n)
  zero_indices <- sample(n, size = floor((1-k)*n), replace = FALSE)
  vector[zero_indices] <- 0
  non_zero_indices <- setdiff(seq_len(n), zero_indices)
  nc = length(non_zero_indices)
  if(distr == 'cauchy'){
    er = rcauchy(nc * p, location = 0.5, scale = 0.5)
  }else if(distr == 'uniform'){
    er = runif(nc * p, min = 0, max = 4)
  }
  Phi[non_zero_indices,] =  matrix(er, nrow = nc, ncol = p)
  noise <- rnorm(n,0,sigma_eps)
  if(distr == 'cauchy'){
    ery = abs(rstable(length(non_zero_indices), alpha = 1.5, beta = 0, gamma = 3, delta = 0))
  }else if(distr == 'uniform'){
    ery = runif(length(non_zero_indices), min = 3, max = 10)
  }
  vector[non_zero_indices] <-   ery 
  return(list(X = Phi, 
              y =  measurements_clean + noise + vector,
              Sigma = sigma,
              signal_true = signal_true,
              non_zero_indices = nonzero_indices0,
              corruption =  vector))
}
generate_target_data2 <- function(i,rr,n,p,s,sigma_eps,k,distr){
  set.seed(i)
  sigma <- diag(p) 
  mu <- rep(0, p)
  for (i in 1:(p-1)) {
    sigma[i, i+1] <- rr 
    sigma[i+1, i] <- rr 
  }
  Phi <- mvrnorm(n , mu = mu, Sigma = sigma) 
  signal_true <- rep(0, p)
  
  nonzero_indices0 <- sample(p, s)
  signal_true[nonzero_indices0] <- sign(rnorm(s))
  measurements_clean <- Phi %*% signal_true
  measurements_ob = measurements_clean
  err_indices = NULL
  if(k > 0){
    err_indices <- sample(n, size = floor((k)*n), replace = FALSE)
    Phi_err <- mvrnorm(floor((k)*n) , mu =  rep(0, p), Sigma = sigma) 
    signal_err = runif(p, -1,1)
    y_err <- Phi_err %*% signal_err
    measurements_ob[err_indices] = y_err
    Phi[err_indices,] = Phi_err
  }
  noise <- rnorm(n,0,sigma_eps)
  
  return(list(X = Phi, 
              y =  measurements_ob  + noise,
              Sigma = sigma,
              signal_true = signal_true,
              non_zero_indices = nonzero_indices0,
              corruption =  as.integer(1:n %in% err_indices)))
}

###############################
#generate source data
###############################

create_auxiliary_data <- function(j,signal_true,sigma_eps,nonzero_indices,rr,e,n,p) {
  set.seed(j)  
  signal_true1 <- c(signal_true, rep(0,p-length(signal_true)))
  a = sample(e,1)
  signal_true1[sample(nonzero_indices, a)] <- 0
  signal_true1[sample(setdiff(1:n,nonzero_indices), a)] <- sign(rnorm(a))
  tau = rbinom(1, size = 1, prob = 1/L)
  signal_true1[sample(which(signal_true1==0),e)] <- tau*sign(rnorm(e))
  mu <- rep(0, p)
  
  sigma <- diag(p) 
  for (i in 1:(p-1)) {
    sigma[i, i+1] <- rr 
    sigma[i+1, i] <- rr 
  }
  Sigma <- matrix(rr, nrow = p, ncol = p)
  diag(Sigma) <- 1  
  
  gauss_cop <- normalCopula(P2p(Sigma), dim = p, dispstr = "un")
  Phi_au  <- qnorm(rCopula(n, gauss_cop) )
  measurements_clean <- Phi_au %*% signal_true1
  noise <- rnorm(n,0,sigma_eps)
  return(list(X = Phi_au, 
              y =  measurements_clean + noise,
              signal = signal_true1,
              error = a))
}

###############################
#subsampling_lasso
###############################
subsampling_lasso_fun <- function(X, y,  X_test, y_test, lambda, 
                                  subsample_lasso0 ,
                                  subsample_ratio = 0.5, max_iter = 1000, tol = 1e-6) {
  n <- nrow(X)       
  p <- ncol(X)       
  subsample_size <- floor(n * subsample_ratio)  
  
  X_center <- colMeans(X)
  X_scale <- apply(X, 2, sd)
  X_scale[X_scale == 0] <- 1  
  X_std <- scale(X, center = TRUE, scale = TRUE)
  y_center <- mean(y)
  y_std <- y - y_center
  
  beta <- subsample_lasso0  
  beta_old <-  beta  
  iter <- 0          
  momentum <- 0.9    
  
  lambda_adj <- lambda
  
  compute_loss <- function(beta, X_test, y_test, lambda) {
    residuals <- y_test - X_test %*% beta
    mse <- mean(residuals^2)
    return(mse )
  }
  
  loss_old <- compute_loss(beta, X_test, y_test, lambda)
  lossvec = rep(0,max_iter)
  
  while (iter < max_iter) {
    set.seed(iter)
    gradient_l1_norms <- apply(-X_std, 1, function(x) sum(abs(x^2))) * abs(y_std - X_std %*% beta)
    gradient_l1_norms[gradient_l1_norms == 0] <- 1e-10
    what <- (1 / gradient_l1_norms) / max(1 / gradient_l1_norms)
    set2 <- setdiff(sample(n, subsample_size, replace = FALSE, prob = what), which(folds == 1))
    subsample_idx <- set2
    X_sub <- X_std[subsample_idx, , drop = FALSE]  
    y_sub <- y_std[subsample_idx]                  
    
    residuals <- y_sub - X_sub %*% beta
    
    beta_new <- beta  
    for (j in 1:p) {
      beta_j_old <- beta[j]  
      rho <- sum(X_sub[, j] * residuals) / subsample_size + beta[j]
      
      beta_new[j] <- sign(rho) * max(0, abs(rho) - lambda_adj)
      
      beta_new[j] <- momentum * beta[j] + (1 - momentum) * beta_new[j]
      
      residuals <- residuals + X_sub[, j] * (beta_j_old - beta_new[j])
    }
    
    lossvec[iter] = loss_old
    iter0  = max(iter-100,1)
    vec0  = lossvec[-iter0]
    if (min(vec0) > lossvec[iter0] * 0.9) {  
      cat("Loss increased at iteration", iter + 1, ", reverting to previous beta\n")
      beta_new <- subsample_lasso0
    }
    
    beta <- beta_new
    
    beta_old <- beta
    iter <- iter + 1
  }
  
  if (iter == max_iter) {
    cat("Reached maximum iterations without convergence\n")
  }
  
  return(beta)
}

###############################
#q-aggregation
###############################

aggre_fun <- function(y_train, X_train, B){
  theta.hat <- exp(-colSums((y_train - X_train %*% B)^2) / 2)
  theta.hat <- theta.hat / sum(theta.hat)
  theta.old <- theta.hat
  
  if (sum(is.na(theta.old)) > 0) {
    w_init <- rep(0, ncol(B))
    
    optfun <- function(w, y_test, X_test, B) {
      sum((y_test - X_test %*% (B %*% w))^2)
    }
    
    result <- optim(par = w_init, fn = optfun, method = "BFGS", 
                    y_test = y_train, X_test = X_train, B = B)
    
    w_optimal <- exp(result$par)
    w_optimal <- w_optimal / sum(w_optimal)
    beta <- as.numeric(B %*% w_optimal)
  } else {
    beta <- as.numeric(B %*% theta.hat)
    total.step <- 10
    
    for (ss in 1:total.step) {
      theta.hat <- exp(-colSums((y_train - X_train %*% B)^2) / 2 + 
                         colSums((as.vector(X_train %*% beta) - X_train %*% B)^2) / 8)
      theta.hat <- theta.hat / sum(theta.hat)
      
      beta <- as.numeric(B %*% theta.hat * 1/4 + 3/4 * beta)
      
      theta.old <- theta.hat
    }
  }
  
  return(list(theta=theta.hat,beta=beta))
}

###############################
# 
###############################
optfun <- function(w, y_test, X_test, B) {
  l = length(w)
  n = nrow(X_test)
  w_transformed = exp(w)
  w_transformed = w_transformed / sum(w_transformed)
  length(w_transformed)
  obj =  mean((y_test - as.vector(X_test%*%B %*% w_transformed))^2) + 
    mean((y_test - X_test%*%B )^2%*%w_transformed) + 2*log(l)*sum(abs(w_transformed))/n
  return(obj )
}

###############################
# SDL
###############################
SDL = function(Xs0,ys,X0,y0,r2,r0,w_hat_A){
  L0 = length(w_hat_A)
  hhat =  h0 = rep(0,L0)
  for(j in nu){
    w_hat <- w_hat_A[[j]]
    if(r2>r0){
      X1 = rbind(Xs0[[j]],X0[,posvec[[j]]])
      Y1 = c(ys[[j]],y0)
      p = dim(X1)[2]
      n = dim(X1)[1]
      q <- Variable(p)
      r <- Variable(n)
      tau_x = sqrt((log(p)/n))
      tau_r = sqrt((log(n)/n))
      objective <- Minimize((1/2) * sum((Y1  - X1 %*%q - r)^2) + tau_r* sum(abs(r)))
      problem <- Problem( objective)
      result <- solve(problem)
      x0 = result$getValue(q)
      hhat[j]  <-  2*sum(abs(x0 - w_hat ))
    }else{hhat[j]  <-  sum(abs(rlasso[posvec[[j]]] - w_hat ))}
    h0[j] = sum(abs(signal_true - signal[[j]][1:400]))
  }
  hhat[!c(1:L)%in%nu]=200
  cbind(hhat,h0)
  return(hhat)
}

###############################
# inference
###############################
CI_fun <- function(db_lasso,Theta_glasso,noise_level,cov0,n){
  var_beta <- noise_level * diag(Theta_glasso %*% cov0 %*% t(Theta_glasso)) / n
  var_beta <- pmax(var_beta, 1e-10)  
  se_beta <- sqrt(var_beta)
  p = length(db_lasso)
  z_stats <- db_lasso / se_beta
  p_values <- 2 * (1 - pnorm(abs(z_stats)))
  p_order <- order(p_values)  
  p_sorted <- p_values[p_order]
  alpha_holm <- alpha / (p - seq_len(p) + 1)  
  reject <- numeric(p)
  for (i in seq_len(p)) {
    if (p_sorted[i] > alpha_holm[i]) break
    reject[p_order[i]] <- 1
  }
  alpha_adj <- rep(alpha, p)
  alpha_adj[p_order] <- alpha_holm
  z_holm <- qnorm(1 - alpha_adj / 2)
  
  L <- db_lasso - z_holm * se_beta
  U <- db_lasso + z_holm * se_beta
  list(L = L,
       U = U,
       pval = p_values,
       se = se_beta,
       reject = reject)
}

CI_fun2 <- function(se_beta,db_lasso){
  z_stats <- db_lasso / se_beta
  p_values <- 2 * (1 - pnorm(abs(z_stats)))
  p_order <- order(p_values)  
  p_sorted <- p_values[p_order]
  p = length(db_lasso)
  alpha_holm <- alpha / (p - seq_len(p) + 1)  
  reject <- numeric(p)
  for (i in seq_len(p)) {
    if (p_sorted[i] > alpha_holm[i]) break
    reject[p_order[i]] <- 1
  }
  alpha_adj <- rep(alpha, p)
  alpha_adj[p_order] <- alpha_holm
  z_holm <- qnorm(1 - alpha_adj / 2)
  p_adj <- numeric(p)
  for (i in seq_len(p)) {
    p_adj[p_order[i]] <- max(p_sorted[i] * (p - i + 1), p_adj[p_order[max(1, i-1)]])
  }
  L <- db_lasso - z_holm * se_beta
  U <- db_lasso + z_holm * se_beta
  list(L = L,
       U = U,
       pval = p_adj,
       reject = reject)
}
CI_fun_BH <- function(db_lasso, Theta_glasso, noise_level, cov0,n,alpha = 0.05) {
  var_beta <- noise_level * diag(Theta_glasso %*% cov0 %*% t(Theta_glasso)) / n
  var_beta <- pmax(var_beta, 1e-10)  
  se_beta <- sqrt(var_beta)
  
  p <- length(db_lasso)
  
  z_stats <- db_lasso / se_beta
  p_values <- 2 * (1 - pnorm(abs(z_stats)))
  
  p_sorted <- sort(p_values)  
  indices <- order(p_values)   
  
  reject_bh <- numeric(p)  
  critical_values <- (seq_len(p) / p) * alpha  
  
  significant <- p_sorted <= critical_values
  if (any(significant)) {
    k_max <- max(which(significant))
    reject_bh[indices[1:k_max]] <- 1
  }
  
  if (any(significant)) {
    bh_threshold <- p_sorted[k_max]  
    alpha_adj <- rep(bh_threshold, p)  
  } else {
    alpha_adj <- rep(alpha/p, p)  
  }
  
  z_bh <- qnorm(1 - alpha_adj / 2)
  
  L <- db_lasso - z_bh * se_beta
  U <- db_lasso + z_bh * se_beta
  
  list(L = L,
       U = U,
       pval = p_values,
       se = se_beta,
       reject = reject_bh,
       fdr_threshold = if (any(significant)) p_sorted[k_max] else alpha/p,
       significant_count = sum(reject_bh))

}
