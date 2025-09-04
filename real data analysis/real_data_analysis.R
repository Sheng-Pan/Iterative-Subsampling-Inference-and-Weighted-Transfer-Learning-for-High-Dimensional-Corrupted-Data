library(data.table)  # 用于高效处理大数据集的库
library(dplyr)       # 数据操作和转换的库
library(tidyverse)   # 数据处理和可视化的综合库
library(glmnet)
library(scalreg)
library(corpcor)
library(glasso)
source("MIP.R")
source("functions.R")

# run data_preprocessing.R
###############################
#traget and source data
###############################
setwd('D:/gbm/preprocessing_data')

target_data = readRDS('target_data.rds')
source1_data = readRDS('source1_data.rds')
source2_data = readRDS('source2_data.rds')

x_s1_d = as.matrix(source1_data[, -which(names(source1_data) == "Y")])
y_s1 = source1_data$Y

x_s2_d = as.matrix(source2_data[, -which(names(source2_data) == "Y")])
y_s2 = source2_data$Y

X_target = as.matrix(target_data[, -which(names(target_data) == "Y")])
Y_target = target_data$Y
dim(X01)
length(y0)
##########################
#Robust subsampling Lasso
###############################
n0 = length(Y_target)
set.seed(123)
set3= sample(n0,floor(n0*0.8))
X_set3 = X_target[set3,]
y_set3 = Y_target[set3]

# 
cv.lasso <- cv.glmnet( X_set3, y_set3, alpha=1) 
fit_w <- glmnet(X_set3, y_set3  , alpha = 1,lambda = cv.lasso$lambda.min)
subsample_lasso0 <- as.numeric(coef(fit_w)[-1])
gradient_l1_norms <- apply(-X_target, 1, function(x) sum(abs(x^2))) * abs(
  Y_target- X_target %*% subsample_lasso0)
dim(X_target)
set4 =    order(gradient_l1_norms)[floor(n0*0.1):floor(n0*0.9)]#uniform
X_set4 = X_target[set4,]
y_set4 = Y_target[set4]
cv.lasso <- cv.glmnet( X_set4, y_set4, alpha=1) 
fit_w <- glmnet(X_set4, y_set4  , alpha = 1,lambda = cv.lasso$lambda.min)
subsample_lasso0 <- as.numeric(coef(fit_w)[-1])
gradient_l1_norms <- apply(-X_set3, 1, function(x) sum(abs(x^2))) * abs(y_set3- X_set3 %*% subsample_lasso0)
# 
val_set =   sample(set4,20)
X_val = X_target[val_set,]
y_val = Y_target[val_set]
set5= setdiff(set4,val_set)
X_set5 = X_target[set5,]
y_set5 = Y_target[set5]
cv.lasso <- cv.glmnet( X_set5, y_set5, alpha=1)
fit_w <- glmnet(X_set5, y_set5  , alpha = 1,lambda = cv.lasso$lambda.min)
subsample_lasso00 <- as.numeric(coef(fit_w)[-1])
subsample_size = 50
B = 20
beta = subsample_lasso00
X02 = X_target[-val_set,]
y02 = Y_target[-val_set]
length(y02)
dim(X_target)
sr = seq(0.75,0.95,0.05)
qtable = c()
for(j in 1:B){
  n0 = length(y02)
  if(j <= length(sr)){
    set.seed(1)
    sr_0 = sr[j]
    subsample_size = floor(n0*sr[j])
  }else{  
    set.seed(j)
    pos = which(qtable[,2] == qtable[,2][which.min(qtable[,2])])
    sr_0 = qtable[,1][max(pos)]
    subsample_size =  floor(n0*sr_0)
  }
  gradient_l1_norms <- apply(-X02, 1, function(x) sum(abs(x^2))) * abs(y02- X02%*%  beta )
  gradient_l1_norms[gradient_l1_norms == 0] <- 1e-10
  #gradient_l1_norms = gradient_l1_norms^2
  what <- exp_normalize( gradient_l1_norms,1e-4)
  # 假设 folds 全为 0（不排除样本）
  n0 = length(y02)
  set2 <- sample(n0, subsample_size, replace = FALSE, prob = what)
  X_set5 = X02[set2,]
  y_set5 = y02[set2]
  cv.lasso <- cv.glmnet( X_set5, y_set5, alpha=1)
  fit_w <- glmnet(X_set5, y_set5  , alpha = 1,lambda = cv.lasso$lambda.min)
  beta_upate  <- as.numeric(coef(fit_w)[-1])
  ##########################
  ##########################
  what <- (1 / gradient_l1_norms) / max(1 / gradient_l1_norms)
  
  set2 <- sample(n0, subsample_size, replace = FALSE, prob = what)
  X_set5 = X02[set2,]
  y_set5 = y02[set2]
  cv.lasso <- cv.glmnet( X_set5, y_set5, alpha=1)
  fit_w <- glmnet(X_set5, y_set5  , alpha = 1,lambda = cv.lasso$lambda.min)
  beta_upate1  <- as.numeric(coef(fit_w)[-1])
  ##
  esset = cbind(beta,beta_upate,beta_upate1)
  w_init = rep(0,dim(esset)[2])
  lossfun = function(x,y,b){
    mean((y - x%*%b)^2)
  }
  # 进行优化
  result <- optim(par = w_init, fn = optfun, method = "BFGS", 
                  y_test = y_val, X_test = X_val, B =  esset)
  
  # 获得满足约束的最优解
  w_optimal = exp(result$par)
  w_optimal = w_optimal / sum(w_optimal)
  beta = esset %*% w_optimal
  if(j <= length(sr)){
    qtable = rbind(qtable,c(sr[j],mean((y_val-X_val%*%beta_upate1)^2)))
  }
  print(qtable)
}
subsample_lasso = beta
##inference
#set2 <- sample(n0, floor(n0*0.8), replace = FALSE, prob = what)
X_set5 = X_target[set2,]
y_set5 = y0[set2]
fit <- scalreg(X_set5, y_set5)
noise_level <- fit$hsigma
# noise_level = min(noise_level1,noise_level2)
cov0 <-  cov.shrink(X_set5)
# cov0 <- cov.shrink(X0)
rho <- 0.1  #
glasso_fit <- glasso(cov0, rho = rho)
# 
Theta_glasso <- glasso_fit$wi
beta0 = beta_upate1
db_lasso1 = beta0 + (1/length(y_set5))*Theta_glasso%*%(t(X_set5))%*%(y_set5 -X_set5%*%beta0 )

CI = CI_fun(db_lasso1,Theta_glasso,noise_level,cov0,length(y_set5))
L2 = CI$L
U2 = CI$U
se0 = CI$se

beta0 = beta_upate
db_lasso2 = beta0 + (1/length(y_set5))*Theta_glasso%*%(t(X_set5))%*%(y_set5 -X_set5%*%beta0 )

CI = CI_fun (db_lasso2,Theta_glasso,noise_level,cov0,length(y_set5))
L1 = CI$L
U1 = CI$U
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

# 
w_optimal = exp(result$par)
w_optimal = w_optimal / sum(w_optimal)
db_lasso = esset %*% w_optimal
se = sqrt(w_optimal[1]^2*se0^2 + w_optimal[2]^2*se1^2 + 2*w_optimal[1]*w_optimal[2]*se0*se1)
CI = CI_fun2(se,db_lasso)
L = CI$L
U = CI$U
reject = CI$reject
colnames(X_target)[which(CI$reject==0)]
cbind(colnames(X_target),subsample_lasso,db_lasso)
colnames(X_target)[which(CI$reject==1)]
####forest plot
genes <- colnames(X_target)
logFC <- db_lasso
lower <- CI$L
upper <- CI$U
p_value <- CI$pval

data <- data.frame(
  Gene = genes,
  Effect_Size = logFC,
  CI_Lower = lower,
  CI_Upper = upper,
  p_adj =  CI$pval
)
ggplot(data, aes(x = Gene, y = Effect_Size, ymin = CI_Lower, ymax = CI_Upper, color = p_adj < 0.05)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +  # 添加零线
  coord_flip() +  # 翻转坐标轴，使基因名称更易读
  scale_color_manual(values = c("grey", "red"), labels = c("p.adj ≥ 0.05", "p.adj < 0.05")) +  # 颜色标注
  labs(x = "Gene", y = "Effect Size (log2 Fold Change)", title = "Forest Plot of Genes") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6, angle = 0, hjust = 1),  # 调整纵坐标字体大小和角度
    legend.position = "right",
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(title = "Significance"))
###############################
#source data selction
###############################
cv.lasso <- cv.glmnet(x_s1_d,y_s1, alpha=1)  # alpha=1 for Lasso
alpha1_fit<-glmnet(x_s1_d,y_s1,alpha=1,lambda = cv.lasso$lambda.min , standardize = FALSE)
lasso1 = as.numeric(coef(alpha1_fit)[-1])
x1_gene = colnames(x_s1_d)[which( abs(lasso1)>0)]

cv.lasso <- cv.glmnet(x_s2_d,y_s2, alpha=1)  # alpha=1 for Lasso
alpha1_fit<-glmnet(x_s2_d,y_s2,alpha=1,lambda = cv.lasso$lambda.min , standardize = FALSE)
lasso2 = as.numeric(coef(alpha1_fit)[-1])
x2_gene = colnames(x_s2_d)[which( abs(lasso2)>0)]


x0_gene = colnames(X_target)[which( abs(subsample_lasso)>0)]

intersect(x0_gene,x1_gene)
intersect(x0_gene,x2_gene)
length(intersect(x0_gene,x2_gene))
# ###############################
# #aggregate
# ###############################
Ahat=c(1,2)
ww = cbind(lasso1,lasso2)
B = ww[,Ahat]
theta.hat<- exp(-colSums((y_train-X_train%*%B)^2)/2)
theta.hat=theta.hat/sum(theta.hat)
theta.old=theta.hat

beta<-as.numeric(B%*%theta.hat)
beta.ew<-beta
total.step = 10
# theta.old=theta.hat
for(ss in 1:total.step){
  theta.hat<- exp(-colSums((y_train-X_train%*%B)^2)/2+colSums((as.vector(X_train%*%beta)-X_train%*%B)^2)/8)
  theta.hat<-theta.hat/sum(theta.hat)
  beta<- as.numeric(B%*%theta.hat*1/4+3/4*beta)
  if(sum(abs(theta.hat-theta.old))<10^(-3)){break}
  theta.old=theta.hat
}
agg_lasso = beta
# ###############################
# #weigted transfer learning
# ###############################
beta = subsample_lasso
B = 20
for(j in 1:B){
  set.seed(j)
  gradient_l1_norms <- apply(-X_target, 1, function(x) sum(abs(x^2))) * abs(Y_target- X_target %*%  beta )
  gradient_l1_norms[gradient_l1_norms == 0] <- 1e-10
  tau = median(gradient_l1_norms ) + 1.34*median(abs(gradient_l1_norms-median(gradient_l1_norms )))
  g0 = pmin(gradient_l1_norms, tau)
  what <- (1 /  g0) / max(1 /  g0)
  set2 = order(what)[1:subsample_size]
  
  ss0 = sample(length(set2),15)
  X_val = X_target[ss0,]
  y_val = Y_target[ss0]
  X_target2 = X_target[-ss0,]
  Y_target2 = Y_target[-ss0]
  weight =  what[-ss0]
  length(weight)
  dim(X_target2)
  y1 = Y_target2 - X_target2%*%agg_lasso
  cv.lasso <- cv.glmnet(X_target2, y1 , alpha=1, weights = weight)
  fit_w <- glmnet(X_target2, y1 , alpha = 1,lambda = cv.lasso$lambda.min, weights = weight)
  delta_hat = as.numeric(coef(fit_w)[-1])
  betahat <- delta_hat + agg_lasso
  beta_upate  <- betahat
  
  esset = cbind(beta,beta_upate)
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
}
RTL = beta
x0_gene = colnames(X_target)[which( abs(RTL)>1e-04)]
length(x0_gene)
sum(abs(RTL)>0)
intersect(x0_gene,x1_gene)
intersect(x0_gene,x2_gene)
var(y0)
# ###############################
# #CV error analysis
# ###############################
set.seed(123)
folds <- sample(rep(1:10, length.out = nrow(X01))) # 创建k个折叠的索引
mse_lasso = mse_SL = mse_WTL = mse_ISHI=mse_FBOD=rep(0,10)
for(k in 1:10){
  set.seed(k)
  train_set = which(folds!=k)
  X_train = X_target[train_set,]
  y_train = Y_target[train_set]
  X_test1 = X_target[-train_set,]
  y_test1 = Y_target[-train_set]
  #
  cv.lasso <- cv.glmnet(X_train,y_train, alpha=1)  # alpha=1 for Lasso
  alpha1_fit<-glmnet(X_train,y_train,alpha=1,lambda = cv.lasso$lambda.min , standardize = FALSE)
  lasso = as.numeric(coef(alpha1_fit)[-1])
  ###############################
  #FBOD
  ###############################
  data = cbind(y_train,X_train)
  scores <- Func.FBOD(data = data, iter = 10, k.nn = 5)
  X_sub1 = X_train[which(scores <2),]
  y_sub1 = y_train[which(scores <2)]
  if(length(which(scores >=2))>0){
    cv.lasso <- cv.glmnet( X_sub1, y_sub1, alpha=1)
    fit_w <- glmnet(X_sub1, y_sub1 , alpha = 1,lambda = cv.lasso$lambda.min)
    FBOD  <- as.numeric(coef(fit_w)[-1])
    GL2 = FBOD
    print(c(length(which(scores >=2)),  sum(abs(GL2-lasso ))))
  }else{
    GL2=lasso
  }
  
  
  ###############################
  #Robust subsampling Lasso
  ###############################
  #pilot estimator
  n0 = length(y_train)
  set3= sample(n0,floor(n0*0.8))
  X_set3 = X_train[set3,]
  y_set3 = y_train[set3]
  
  cv.lasso <- cv.glmnet( X_set3, y_set3, alpha=1) 
  fit_w <- glmnet(X_set3, y_set3  , alpha = 1,lambda = cv.lasso$lambda.min)
  subsample_lasso0 <- as.numeric(coef(fit_w)[-1])
  gradient_l1_norms <- apply(-X_train, 1, function(x) sum(abs(x^2))) * abs(
    y_train- X_train %*% subsample_lasso0)
  dim(X_train)
  set4 =    order(gradient_l1_norms)[floor(n0*0.1):floor(n0*0.9)]#uniform
  X_set4 = X_train[set4,]
  y_set4 = y_train[set4]
  cv.lasso <- cv.glmnet( X_set4, y_set4, alpha=1) 
  fit_w <- glmnet(X_set4, y_set4  , alpha = 1,lambda = cv.lasso$lambda.min)
  subsample_lasso0 <- as.numeric(coef(fit_w)[-1])
  gradient_l1_norms <- apply(-X_set3, 1, function(x) sum(abs(x^2))) * abs(y_set3- X_set3 %*% subsample_lasso0)
  val_set =   sample(set4,15)
  X_val = X_train[val_set,]
  y_val = y_train[val_set]
  set5= setdiff(set4,val_set)
  X_set5 = X_train[set5,]
  y_set5 = y_train[set5]
  cv.lasso <- cv.glmnet( X_set5, y_set5, alpha=1)
  fit_w <- glmnet(X_set5, y_set5  , alpha = 1,lambda = cv.lasso$lambda.min)
  subsample_lasso00 <- as.numeric(coef(fit_w)[-1])
  subsample_size = 50
  B = 20
  beta = subsample_lasso00
  X02 = X_train[-val_set,]
  y02 = y_train[-val_set]
  length(y02)
  #iterative subsampling
  sr = seq(0.9,1,0.05)
  qtable = c()
  for(j in 1:B){
    n0 = length(y02)
    if(j <= length(sr)){
      set.seed(1)
      sr_0 = sr[j]
      subsample_size = floor(n0*sr[j])
    }else{  
      set.seed(j)
      pos = which(qtable[,2] == qtable[,2][which.min(qtable[,2])])
      sr_0 = qtable[,1][max(pos)]
      subsample_size =  floor(n0*sr_0)
    }
    gradient_l1_norms <- apply(-X02, 1, function(x) sum(abs(x^2))) * abs(y02- X02%*%  beta )
    gradient_l1_norms[gradient_l1_norms == 0] <- 1e-10
    #gradient_l1_norms = gradient_l1_norms^2
    what <- (1 / gradient_l1_norms) / max(1 / gradient_l1_norms)
    # 假设 folds 全为 0（不排除样本）
    
    
    set2 <- sample(n0, subsample_size, replace = FALSE, prob = what)
    X_set5 = X02[set2,]
    y_set5 = y02[set2]
    cv.lasso <- cv.glmnet( X_set5, y_set5, alpha=1)
    fit_w <- glmnet(X_set5, y_set5  , alpha = 1,lambda = cv.lasso$lambda.min)
    beta_upate  <- as.numeric(coef(fit_w)[-1])
    ######
    ######
    what <- exp_normalize( gradient_l1_norms,1e-4)
    # 假设 folds 全为 0（不排除样本）
    
    n0 = length(y02)
    set2 <- sample(n0, subsample_size, replace = FALSE, prob = what)
    X_set5 = X02[set2,]
    y_set5 = y02[set2]
    cv.lasso <- cv.glmnet( X_set5, y_set5, alpha=1)
    fit_w <- glmnet(X_set5, y_set5  , alpha = 1,lambda = cv.lasso$lambda.min)
    beta_upate1  <- as.numeric(coef(fit_w)[-1])
    ######
    ######
    esset = cbind(beta,beta_upate,beta_upate1)
    w_init = rep(0,dim(esset)[2])
    lossfun = function(x,y,b){
      mean((y - x%*%b)^2)
    }
    # 进行优化
    result <- optim(par = w_init, fn = optfun, method = "BFGS", 
                    y_test = y_val, X_test = X_val, B =  esset)
    
    # 获得满足约束的最优解
    w_optimal = exp(result$par)
    w_optimal = w_optimal / sum(w_optimal)
    beta = esset %*% w_optimal
    if(j <= length(sr)){
      qtable = rbind(qtable,c(sr[j],mean((y_val-X_val%*%beta_upate1)^2)))
    }
    #print(qtable)
  }
  subsample_lasso = beta
  # ###############################
  # #weighted transfer learning
  # ###############################
  # aggregation
  Ahat=c(1,2)
  ww = cbind(lasso1,lasso2)
  B = ww[,Ahat]
  lossfun = function(x,y,b){
    mean((y - x%*%b)^2)
  }
  result <- optim(par = c(0,0), fn = optfun, method = "L-BFGS-B",
                  y_test = y_train, X_test = X_train, B =  B)
  # 获得满足约束的最优解
  w_optimal = exp(result$par)
  w_optimal = w_optimal / sum(w_optimal)
  agg_lasso = B %*% w_optimal
  # weighted transfer learning
  beta = subsample_lasso
  B = 10
  for(j in 1:B){
    set.seed(j)
    gradient_l1_norms <- apply(-X_train, 1, function(x) sum(abs(x^2))) * abs(y_train- X_train %*% beta)
    gradient_l1_norms[gradient_l1_norms == 0] <- 1e-10
    g0 = gradient_l1_norms#pmin(gradient_l1_norms, tau)
    what <- (1 / g0) / max(1 / g0)
    n = length(y_train)
    set2 <- sample(n, floor(n*0.9), replace = FALSE, prob = what)
    ss0 = sample(length(set2),30)
    X_val = X_train[ss0,]
    y_val = y_train[ss0]
    X_train2 = X_train[-ss0,]
    y_train2 = y_train[-ss0]
    ########################
    what <- exp_normalize(  g0,1e-4)
    weight =  what[-ss0]
    y1 = y_train2 - X_train2%*%agg_lasso
    cv.lasso <- cv.glmnet(X_train2, y1 , alpha=1, weights = weight)
    fit_w <- glmnet(X_train2, y1 , alpha = 1,lambda = cv.lasso$lambda.min, weights = weight)
    lasso0 = as.numeric(coef(fit_w)[-1])+ agg_lasso
    
    ########################
    if(j==1){
      beta = lasso0
    }else{
      esset = cbind(beta,lasso0)
      w_init = c(2,1)
      lossfun = function(x,y,b){
        mean((y - x%*%b)^2)
      }
      # 进行优化
      result <- optim(par = w_init, fn = optfun, method = "BFGS", 
                      y_test = y_val, X_test = X_val, B =  esset)
      
      # 
      w_optimal = exp(result$par)
      w_optimal = w_optimal / sum(w_optimal)
      beta = esset %*% w_optimal
    }
    
  }
  ###############################
  #MI
  ###############################
  k_sub = 0.3
  n_train = length(y_train)
  n_subset = k_sub*n_train+1
  subset_vol = n_train/n_subset
  q = 1
  alpha = 0.05
  infl =  tryCatch({
    
    MIP(X_train, y_train,n_train,p,q,n_subset,subset_vol,ep=0.01,alpha)
  }, error = function(e) {
    
    cat("Error occurred:", e$message, "\n")
    return(NULL)
  })
  if(length(infl$inf_setfinal)>0){
    sub1 = setdiff(c(1:n_train),infl$inf_setfinal)
    X_set3 = X_train[sub1,]
    y_set3 = y_train[sub1]
    
    # 
    cv.lasso <- cv.glmnet( X_set3, y_set3, alpha=1)
    fit_w <- glmnet(X_set3, y_set3  , alpha = 1,lambda = cv.lasso$lambda.min)
    HIM <- as.numeric(coef(fit_w)[-1])
  }else{
    HIM = lasso
  }
  mse_MI[k] =  mean((y_test1 - (X_test1 %*% HIM))^2)
  # cbind(beta,beta0)
  print(c(mean((y_test1 - (X_test1 %*% subsample_lasso))^2),
          mean((y_test1 - (X_test1 %*% beta))^2),mean((y_test1 - (X_test1 %*% GL2))^2)))
  WTL = beta
  mse_lasso[k] =  mean((y_test1 - (X_test1 %*% lasso))^2)
  #  mse_SL[j] = mean((y_test1 - (X_test1 %*% SL))^2)
  mse_ISHI[k] = mean((y_test1 - (X_test1 %*% subsample_lasso))^2)
  mse_FBOD[k] = mean((y_test1 - (X_test1 %*% GL2))^2)
  mse_WTL[k] = mean((y_test1 - (X_test1 %*% WTL))^2)
}
mean(mse_lasso)

mean(mse_ISHI)
mean(mse_FBOD)
mean(mse_WTL)
mean(mse_MI)

