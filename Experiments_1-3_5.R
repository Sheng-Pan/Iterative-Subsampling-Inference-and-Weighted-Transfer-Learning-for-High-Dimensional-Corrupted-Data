# Simulation Study for Experiments 1-3 and 5
# -------------------------------------------------------------------

# Set working directory with fallback option
setwd('...')
# Load required functions
source("MIP.R")
source("functions.R")


# Load required packages
library(glmnet)    # For Lasso regression
library(CVXR)      # For convex optimization
library(copula)    # For copula models
library(stabledist) # For stable distributions
library(glasso)    # For graphical lasso
library(corpcor)   # For correlation estimation
library(scalreg)   # For scaled lasso
library(dbscan)    # For density-based clustering

# Simulation parameters
N <- 30            # Number of Monte Carlo repetitions
p <- 400           # Dimension of the parameter vector
r_vector <- seq(0, 0.2, 0.05)  # Corruption levels

# Initialize result matrices
initialize_matrix <- function() {
  matrix(0, nrow = length(r_vector), ncol = N)
}

ser_lasso <- initialize_matrix()
ser_subsample_lasso_transfer <- initialize_matrix()
ser_subsample_lasso <- initialize_matrix()
ser_MI_lasso <- initialize_matrix()
ser_subsample_lasso_transfer2 <- initialize_matrix()
mse_subsample_lasso_transfer2 <- initialize_matrix()
mse_lasso <- initialize_matrix()
mse_subsample_lasso_transfer <- initialize_matrix()
mse_subsample_lasso <- initialize_matrix()
mse_MI_lasso <- initialize_matrix()
ser_oracle_subsample_lasso_transfer <- initialize_matrix()
ser_transfer_lasso <- initialize_matrix()
ser_transfer_lasso2 <- initialize_matrix()
ser_LGL_lasso <- initialize_matrix()
mse_LGL_lasso <- initialize_matrix()
mse_oracle_subsample_lasso_transfer <- initialize_matrix()
mse_transfer_lasso <- initialize_matrix()
mse_transfer_lasso2 <- initialize_matrix()
fdr_MI <- initialize_matrix()
fdr <- initialize_matrix()
ser_gra_lasso <- initialize_matrix()
mse_gra_lasso <- initialize_matrix()

# Array initialization for signal estimates
signal_array <- array(0, dim = c(p, length(r_vector), N))
singal_lasso <- signal_array
signal_MI <- signal_array
signal_sl_lasso <- signal_array
signal_RSL <- signal_array
signal_ORSL <- signal_array
signal_LGL_lasso <- signal_array

# Additional result matrices
lasso_avgs <- matrix(0, length(r_vector), 2)
sub_avgs <- matrix(0, length(r_vector), 2)
lasso_avgs_len <- matrix(0, length(r_vector), 2)
sub_avgs_len <- matrix(0, length(r_vector), 2)
sub_test <- matrix(0, length(r_vector), 2)
lasso_test <- matrix(0, length(r_vector), 2)

IFV <- list()

# Main simulation loop
for (u in 1:length(r_vector)) {
  k <- r_vector[u]
  aaa <- list()
  
  # Temporary storage for Monte Carlo results
  lasso_avgs1 <- matrix(0, N, 2)
  sub_avgs1 <- matrix(0, N, 2)
  lasso_avgs_len1 <- matrix(0, N, 2)
  sub_avgs_len1 <- matrix(0, N, 2)
  sub_test1 <- matrix(0, N, 2)
  lasso_test1 <- matrix(0, N, 2)
  
  for (m in 1:N) {
    set.seed(m)
    t1 <- Sys.time()
    
    # Data generation parameters
    s <- 12           # Sparsity level
    n <- 300          # Sample size
    sigma <- 1        # Natural noise level
    rr <- 0.3
    scenario <- 2
    
    # Generate data based on scenario
    if (scenario == 1) {
      distr <- 'uniform'  # Experiment 2
      gdata <- generate_target_data(m, rr, n, p, s, sigma, k, distr)
    } else if (scenario == 2) {
      gdata <- generate_target_data2(m, rr, n, p, s, sigma, k, distr)  # Experiment 3
    } else if (scenario == 3) {
      distr <- 'cauchy'
      gdata <- generate_target_data0(m, rr, n, p, s, sigma, k, distr)  # Experiment 1
    }
    
    # Extract data components
    X0 <- gdata$X
    y0 <- gdata$y
    signal_true <- gdata$signal_true
    S <- which(abs(signal_true) > 0)
    Sc <- which(abs(signal_true) == 0)
    alpha <- 0.05
    
    # ###############################
    # #FBOD method
    # ###############################
    data = cbind(y0,X0)
    scores <- Func.FBOD(data = data, iter = 10, k.nn = 5)
    X_sub1 = X0[which(scores <2),]
    y_sub1 = y0[which(scores <2)]
    cv.lasso <- cv.glmnet( X_sub1, y_sub1, alpha=1)
    fit_w <- glmnet(X_sub1, y_sub1 , alpha = 1,lambda = cv.lasso$lambda.min)
    FBOD  <- as.numeric(coef(fit_w)[-1])
    # ###############################
    # #MI method
    # ###############################
    # 
    k_sub = 0.3
    n_subset = k_sub*n+1
    subset_vol = n/n_subset
    q = 1
    alpha = 0.05
    infl =  tryCatch({
      MIP(X0, y0,n,p,q,n_subset,subset_vol,ep=0.01,alpha)
    }, error = function(e) {
      cat("Error occurred:", e$message, "\n")
      return(NULL)
    })
    infl
    aaa[[m]] = infl
    sub1 = setdiff(c(1:n),infl$inf_setfinal)
    X_set3 = X0[sub1,]
    y_set3 = y0[sub1]
    # 
    cv.lasso <- cv.glmnet( X_set3, y_set3, alpha=1)
    fit_w <- glmnet(X_set3, y_set3  , alpha = 1,lambda = cv.lasso$lambda.min)
    HIM <- as.numeric(coef(fit_w)[-1])
    ser(signal_true,  HIM)
    ser_MI_lasso[u,m] =  ser(signal_true,  HIM)
    mse_MI_lasso[u,m] =  msefun(signal_true,  HIM)
    signal_MI[,u,m] = HIM
    
    # ###############################
    # #ISHI method
    # ###############################
    # 
    result_9 = sample_fun(X0,y0)
    subsample_lasso00 = result_9$subsample_lasso00
    val_set = result_9$val_set
    B = 20
    beta = subsample_lasso00
    X01 = X0[-val_set,]
    y01 = y0[-val_set]
    X_val = X0[val_set,]
    y_val = y0[val_set]
    length(y_train)
    M = matrix(0,B,p)
    S = which(abs(signal_true)>0)
    Sc = which(abs(signal_true)==0)
    z_alpha = qnorm(1-alpha/2)
    # 
    qtable = c()
    sr = seq(0.7,0.85,0.05)
    for(b in 1:B){
      
      if(b <= 4){
        set.seed(1)
        sr_0 = sr[b]
        subsample_size = floor(n*sr[b])
      }else{  
        set.seed(b)
        pos = which(qtable[,2] == qtable[,2][which.min(qtable[,2])])
        sr_0 = qtable[,1][max(pos)]
        subsample_size =  floor(n*sr_0)
      }
      gradient_l1_norms <- gradient_norm(X01, y01, beta)
      # 
      gradient_l1_norms[gradient_l1_norms == 0] <- 1e-10
      # 
      what <- exp_normalize( gradient_l1_norms,1e-04)
      what = what + min(what[!what==0])
      # 
      n0 = length(y01)
      set2 <- sample(n0, subsample_size, replace = FALSE, prob = what)
      ######
      X_set5 = X01[set2,]
      y_set5 = y01[set2]
      
      cv.lasso <- cv.glmnet( X_set5, y_set5, alpha=1)
      fit_w <- glmnet(X_set5, y_set5  , alpha = 1,lambda = cv.lasso$lambda.min)
      beta_upate  <- as.numeric(coef(fit_w)[-1])
      ######
      what <- (1 / gradient_l1_norms) / max(1 / gradient_l1_norms)
      
      set2 <- sample(n0, subsample_size, replace = FALSE, prob = what)
      X_set5 = X01[set2,]
      y_set5 = y01[set2]
      cv.lasso <- cv.glmnet( X_set5, y_set5, alpha=1)
      fit_w <- glmnet(X_set5, y_set5  , alpha = 1,lambda = cv.lasso$lambda.min)
      beta_upate1  <- as.numeric(coef(fit_w)[-1])
      
      ##
      esset = cbind(beta,beta_upate,beta_upate1)
      w_init = rep(0,dim(esset)[2])
      lossfun = function(x,y,b){
        mean((y - x%*%b)^2)
      }
      # 
      result <- optim(par = w_init, fn = optfun, method = "BFGS", 
                      y_test = y_val, X_test = X_val, B =  esset)
      
      # 
      w_optimal = exp(result$par)
      w_optimal = w_optimal / sum(w_optimal)
      beta = esset %*% w_optimal
      if(b <= 4){
        qtable = rbind(qtable,c(sr[b],mean((y_val-X_val%*%beta_upate1)^2)))
      }
      print(c(sr_0,ser(signal_true,  beta),mean((y_val-X_val%*%beta)^2)))
      M[b,] =  beta
    }
    ##########################
    subsample_lasso = beta
    set2 <- sample(n0, floor(n0*0.7), replace = FALSE, prob = what)
    X_set5 = X01[set2,]
    y_set5 = y01[set2]
    #noise_level estimation
    fit <- scalreg(X_set5, y_set5)
    noise_level <- fit$hsigma
    #
    cov0 <-  cov.shrink(X_set5)
    rho <- 0.1  
    glasso_fit <- glasso(cov0, rho = rho)
    Theta_glasso <- glasso_fit$wi
    #
    beta0 = beta_upate1
    db_lasso1 = beta0 + (1/length(y_set5))*Theta_glasso%*%(t(X_set5))%*%(y_set5 -X_set5%*%beta0 )
    #
    CI = CI_fun(db_lasso1,Theta_glasso,noise_level,cov0,length(y_set5))
    se0 = CI$se
    
    beta0 = beta_upate
    db_lasso2 = beta0 + (1/length(y_set5))*Theta_glasso%*%(t(X_set5))%*%(y_set5 -X_set5%*%beta0 )
    
    CI = CI_fun (db_lasso2,Theta_glasso,noise_level,cov0,length(y_set5))
    se1 = CI$se
    ##
    esset = cbind(db_lasso1,db_lasso2)
    w_init = rep(0,dim(esset)[2])
    lossfun = function(x,y,b){
      mean((y - x%*%b)^2)
    }
    # 
    result <- optim(par = w_init, fn = optfun, method = "BFGS", 
                    y_test = y_val, X_test = X_val, B =  esset)
    
    w_optimal = exp(result$par)
    w_optimal = w_optimal / sum(w_optimal)
    db_lasso = esset %*% w_optimal
    ser(signal_true,  db_lasso)
    se = sqrt(w_optimal[1]^2*se0^2 + w_optimal[2]^2*se1^2 + 2*w_optimal[1]*w_optimal[2]*se0*se1)
    CI = CI_fun2(se,db_lasso)
    L = CI$L
    U = CI$U
    
    avgs = sum((L[S]<=signal_true[S])&(U[S]>=signal_true[S]))/length(S)
    avgsc = sum((L[Sc]<=signal_true[Sc])&(U[Sc]>=signal_true[Sc]))/length(Sc)
    U[S] - L[S]
    reject = CI$reject
    power = sum(reject[S])/length(S)
    FDR = length(intersect(which(reject==1),Sc))/max(1,sum(reject))
    print(c(avgs,avgsc,power,FDR,mean(U[S] - L[S]),mean(U[Sc] - L[Sc])))
    sub_test1[m,]=c(power,FDR)
    sub_avgs1[m,]=c(avgs,avgsc)
    sub_avgs_len1[m,]=c(mean(U[S] - L[S]),mean(U[Sc] - L[Sc]))
    ser_subsample_lasso[u,m] = ser(signal_true,  subsample_lasso)
    mse_subsample_lasso[u,m] =  msefun(signal_true,  subsample_lasso)
    ###############################
    #lasso
    ###############################
    cv.lasso <- cv.glmnet(X0, y0, alpha=1)  # alpha=1 for Lasso
    fit_w <- glmnet(X0, y0, alpha = 1,lambda =cv.lasso$lambda.min)
    lasso <- as.numeric(coef(fit_w)[-1])
    ser(signal_true,  lasso)
    ser_lasso[u,m] = ser(signal_true,  lasso)
    mse_lasso[u,m] = msefun(signal_true,  lasso)
    singal_lasso[,u,m] = lasso
    ##
    # noise_level = sum((y0-X0%*%lasso)^2)/ (length(y0) - sum(lasso!= 0))
    # fit <- scaled_lasso(X0, y0, cv.lasso$lambda.min)
    # noise_level <- fit$sigma^2
    cov0 <- cov.shrink(X0)
    fit <- scalreg(X0, y0)
    noise_level <- fit$hsigma
    rho <- 0.1  # 正则化参数
    glasso_fit <- glasso(cov0, rho = rho)
    
    # 
    Theta_glasso <- glasso_fit$wi
    db_lasso = lasso+ (1/length(y0))*Theta_glasso%*%(t(X0))%*%(y0 -X0%*%lasso)
    CI = CI_fun(db_lasso,Theta_glasso,noise_level,cov0,length(y0))
    L = CI$L
    U = CI$U
    lasso[S]
    L[S]
    U[S]
    U[S] - L[S]
    avgs = sum((L[S]<=signal_true[S])&(U[S]>=signal_true[S]))/length(S)
    avgsc = sum((L[Sc]<=signal_true[Sc])&(U[Sc]>=signal_true[Sc]))/length(Sc)
    reject = CI$reject
    power = sum(reject[S])/length(S)
    FDR = length(intersect(which(reject==1),Sc))/sum(reject)
    ser(signal_true,  db_lasso )
    lasso_test1[m,]=c(power,FDR)
    lasso_avgs1[m,]=c(avgs,avgsc)
    lasso_avgs_len1[m,]=c(mean(U[S] - L[S]),mean(U[Sc] - L[Sc]))
    ser_lasso[u,m] = ser(signal_true,  lasso)
    mse_lasso[u,m] = msefun(signal_true,  lasso)
    ###############################
    #weighted transfer learning
    ###############################
    ###generate source data
    ns = 300
    e=8
    L = 20
    sigma = 0.1
    nonzero_indices = which(!signal_true == 0 )
    auxiliary_data <- lapply(seq_len(L), function(i) create_auxiliary_data(i+m,signal_true,sigma,nonzero_indices,rr,e,ns,p))
    # 
    get_X_safe <- function(data_list, component_name) {
      if (!is.null(data_list) && is.list(data_list) && component_name %in% names(data_list)) {
        return(data_list[[component_name]])
      } else {
        stop(paste("Component", component_name, "not found in data list."))
      }
    }
    
    #
    Xs <- lapply(auxiliary_data, get_X_safe, component_name = "X")
    ys <- lapply(auxiliary_data, get_X_safe, component_name = "y")
    signal <- lapply(auxiliary_data, get_X_safe, component_name = "signal")
    error <- unlist(lapply(auxiliary_data, get_X_safe, component_name = "error"))
    h0 = rep(0,L)
    for(j in 1:L){
      h0[j] = sum(abs(signal_true - signal[[j]]))
    }
    # lasso on source
    w_hat_A <- list()
    for(j in 1:L){
      cv.lasso <- cv.glmnet(Xs[[j]], ys[[j]], alpha=1)  # alpha=1 for Lasso
      fit_w <- glmnet(Xs[[j]], ys[[j]], alpha = 1, lambda = cv.lasso$lambda.min)
      w <- as.numeric(coef(fit_w)[-1]) # 
      w_hat_A[[j]] <- w
    }
    # source data selection
    hhat <- h0<-rep(0,L)
    signal0 = ww = delta = matrix(0,p,L)
    
    for(j in 1:L){
      X2 = Xs[[j]]
      y2 = ys[[j]] - X2%*%subsample_lasso
      cv.lasso <- cv.glmnet(X2, y2 , alpha=1)
      fit_w <- glmnet(X2, y2  , alpha = 1,lambda = cv.lasso$lambda.min)
      temp <- as.numeric(coef(fit_w)[-1])
      delta[,j] <- temp
      hhat[j] =  sum(abs(temp))
      h0[j] = sum(abs(signal_true - signal[[j]]))
      signal0[,j] = signal[[j]]
      ww[,j] = w_hat_A[[j]]
    }
    #
    cbind(hhat,h0)
    Ahat = intersect(order(hhat)[1:10],which(hhat<15))
    test_data = which.min(hhat)
    h0[Ahat]
    # # #source data aggregation
    gradient_l1_norms <-  gradient_norm(X0, y0, subsample_lasso )
    index1 = order(gradient_l1_norms)[1:floor(n*0.8)]
    X_train = X0[index1,]
    y_train = y0[index1]
    
    if(length(Ahat)==0){
      subsample_tran_lasso = rep(0,p)
    }else{
      if(length(Ahat)==1){
        agg_lasso = ww[,Ahat]
      }else{
        B = ww[,Ahat]
        result = aggre_fun(y_train, X_train, B)
        agg_lasso = result$beta
        ser(signal_true,agg_lasso)
      }
      ##weighted iteration
      subsample_size = floor(n*0.8)
      beta = subsample_lasso
      B = 10
      for(j in 1:B){
        set.seed(j)
        gradient_l1_norms <-  gradient_norm(X0, y0, beta)
        gradient_l1_norms[gradient_l1_norms == 0] <- 1e-10
        tau = median(gradient_l1_norms ) + 5*median(abs(gradient_l1_norms-median(gradient_l1_norms )))
        g0 = gradient_l1_norms#pmin(gradient_l1_norms, tau)
        what <- (1 /  g0) / max(1 /  g0)
        set2 = order(what)[1:subsample_size]
        ss0 = sample(set2,30)
        val_set = ss0
        X_val = X0[val_set,]
        y_val = y0[val_set]
        X_train = X0[-val_set,]
        y_train = y0[-val_set]
        weight =  what[-val_set]
        length(weight)
        dim(X_train)
        y1 = y_train - X_train%*%agg_lasso
        cv.lasso <- cv.glmnet(X_train, y1 , alpha=1, weights = weight)
        fit_w <- glmnet(X_train, y1 , alpha = 1,lambda = cv.lasso$lambda.min, weights = weight)
        delta_hat = as.numeric(coef(fit_w)[-1])
        betahat <- delta_hat + agg_lasso
        ser(signal_true,  betahat)
        beta_upate  <- betahat
        ########################
        c = -log(0.9)/max(g0)
        what <- exp_normalize(  g0,1e-4)
        weight =  what[-val_set]
        y1 = y_train - X_train%*%agg_lasso
        cv.lasso <- cv.glmnet(X_train, y1 , alpha=1, weights = weight)
        fit_w <- glmnet(X_train, y1 , alpha = 1,lambda = cv.lasso$lambda.min, weights = weight)
        lasso0 = as.numeric(coef(fit_w)[-1])+ agg_lasso
        ser(signal_true, lasso0)
        ########################
        mean((y_val - X_val%*%lasso0)^2)
        mean((y_val - X_val%*%betahat)^2)
        if(j==1){
          beta = beta_upate
        }else{
          esset = cbind(beta,lasso0)
          w_init = c(1,1)
          lossfun = function(x,y,b){
            mean((y - x%*%b)^2)
          }
          result <- optim(par = w_init, fn = optfun, method = "L-BFGS-B",
                          y_test = y_val, X_test = X_val, B =  esset)
          
          # 
          w_optimal = exp(result$par)
          w_optimal = w_optimal / sum(w_optimal)
          beta = esset %*% w_optimal
        }
        
        print(ser(signal_true,  beta))
      }
      subsample_tran_lasso = beta
      
    }
    
    
    ser_subsample_lasso_transfer[u,m] = ser(signal_true,  subsample_tran_lasso)
    mse_subsample_lasso_transfer[u,m] = msefun(signal_true,  subsample_tran_lasso)
    signal_RSL[,u,m] = subsample_tran_lasso
    
    # ###############################
    # #oracle weighted transfer learning
    # ###############################
    # 
    hhat <- h0<-rep(0,L)
    signal0 = ww = matrix(0,p,L)
    for(j in 1:L){
      h0[j] = sum(abs(signal_true - signal[[j]]))
      signal0[,j] = signal[[j]]
      ww[,j] = w_hat_A[[j]]
    }
    #
    cbind(hhat,h0)
    Ah = which(h0<9)
    test_data = which.min(h0)
    h0[Ahat]
    ###
    gradient_l1_norms <- gradient_norm(X0, y0, subsample_lasso)
    index1 = order(gradient_l1_norms)[1:floor(n*0.8)]
    X_train = X0[index1,]
    y_train = y0[index1]
    B = ww[,Ah]
    result = aggre_fun(y_train,X_train,B)
    agg_lasso=result$beta
    ###
    
    beta = subsample_lasso
    B = 10
    for(j in 1:B){
      set.seed(j)
      gradient_l1_norms <- gradient_norm(X0, y0, beta)
      gradient_l1_norms[gradient_l1_norms == 0] <- 1e-10
      tau = median(gradient_l1_norms ) + 5*median(abs(gradient_l1_norms-median(gradient_l1_norms )))
      g0 = gradient_l1_norms#pmin(gradient_l1_norms, tau)
      what <- (1 /  g0) / max(1 /  g0)
      ss0 = sample(set2,30)
      sort(ss0)
      val_set = ss0
      X_val = X0[val_set,]
      y_val = y0[val_set]
      X_train = X0[-val_set,]
      y_train = y0[-val_set]
      weight =  what[-val_set]
      length(weight)
      dim(X_train)
      y1 = y_train - X_train%*%agg_lasso
      cv.lasso <- cv.glmnet(X_train, y1 , alpha=1, weights = weight)
      fit_w <- glmnet(X_train, y1 , alpha = 1,lambda = cv.lasso$lambda.min, weights = weight)
      delta_hat = as.numeric(coef(fit_w)[-1])
      betahat <- delta_hat + agg_lasso
      ser(signal_true,  betahat)
      beta_upate  <- betahat
      ########################
      c = -log(0.9)/max(g0)
      what <- exp_normalize(  g0,1e-4)
      weight =  what[-val_set]
      y1 = y_train - X_train%*%agg_lasso
      cv.lasso <- cv.glmnet(X_train, y1 , alpha=1, weights = weight)
      fit_w <- glmnet(X_train, y1 , alpha = 1,lambda = cv.lasso$lambda.min, weights = weight)
      lasso0 = as.numeric(coef(fit_w)[-1])+ agg_lasso
      ser(signal_true, lasso0)
      ########################
      mean((y_val - X_val%*%lasso0)^2)
      mean((y_val - X_val%*%betahat)^2)
      if(j==1){
        beta = beta_upate
      }else{
        esset = cbind(beta,lasso0)
        w_init = c(1,1)
        lossfun = function(x,y,b){
          mean((y - x%*%b)^2)
        }
        result <- optim(par = w_init, fn = optfun, method = "L-BFGS-B",
                        y_test = y_val, X_test = X_val, B =  esset)
        
        # 
        w_optimal = exp(result$par)
        w_optimal = w_optimal / sum(w_optimal)
        beta = esset %*% w_optimal
      }
      
     # print(ser(signal_true,  beta))
    }
    oracle_subsample_tran_lasso = beta
    ser(signal_true,  oracle_subsample_tran_lasso)
    ser_oracle_subsample_lasso_transfer[u,m] = ser(signal_true,  oracle_subsample_tran_lasso)
    mse_oracle_subsample_lasso_transfer[u,m] = msefun(signal_true,  oracle_subsample_tran_lasso)
    signal_ORSL [,u,m] = oracle_subsample_tran_lasso
    print(data.frame(k = k,
                     m=m,
                     ser_lasso = ser_lasso[u,m],
                     ser_subsample_lasso = ser_subsample_lasso[u,m],
                     ser_MI_lasso = ser_MI_lasso[u,m],
                     ser_subsample_tran_lasso  = ser_subsample_lasso_transfer[u,m]))
      # Record computation time
      t2 <- Sys.time()
      print(t2 - t1)
  }
  # Aggregate results across Monte Carlo repetitions
  sub_test[u, ] <- apply(sub_test1[1:30, ], 2, mean)
  lasso_test[u, ] <- apply(lasso_test1, 2, mean)
  lasso_avgs[u, ] <- apply(lasso_avgs1, 2, mean)
  lasso_avgs_len[u, ] <- apply(lasso_avgs_len1, 2, mean)
  sub_avgs[u, ] <- apply(sub_avgs1, 2, mean)
  sub_avgs_len[u, ] <- apply(sub_avgs_len1, 2, mean)
  
  # Create summary results
  a <- data.frame(
    k = r_vector,
    ser_lasso = apply(ser_lasso, 1, mean),
    ser_subsample_lasso = apply(ser_subsample_lasso, 1, mean),
    ser_MI_lasso = apply(ser_MI_lasso, 1, mean),
    ser_LGL_lasso = apply(ser_LGL_lasso, 1, mean),
    ser_subsample_tran_lasso = apply(ser_subsample_lasso_transfer, 1, mean),
    ser_oracle_subsample_tran_lasso = apply(ser_oracle_subsample_lasso_transfer, 1, mean)
  )
  
  b <- data.frame(
    k = r_vector,
    mse_lasso = apply(mse_lasso, 1, mean),
    mse_subsample_lasso = apply(mse_subsample_lasso, 1, mean),
    mse_MI_lasso = apply(mse_MI_lasso, 1, mean),
    mse_LGL_lasso = apply(mse_LGL_lasso, 1, mean),
    mse_subsample_tran_lasso = apply(mse_subsample_lasso_transfer, 1, mean),
    mse_oracle_subsample_tran_lasso = apply(mse_oracle_subsample_lasso_transfer, 1, mean)
  )
  
  
  # Save final results
  results <- list(
    ser_results = a,
    mse_results = b,
    sub_test = sub_test,
    lasso_test = lasso_test,
    lasso_avgs = lasso_avgs,
    sub_avgs = sub_avgs,
    lasso_avgs_len = lasso_avgs_len,
    sub_avgs_len = sub_avgs_len
  )
  
  #saveRDS(results, file = "simulation_results.rds")
}
    
  