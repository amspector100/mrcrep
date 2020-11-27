library(MASS)
library(glmnet)
library(knockoff)
library(foreach)
library(doParallel)

set.seed(123456)  # Set seed for reproducibility
n <- 300  # Number of observations
p <- 1000  # Number of predictors included in model
q <- 60  # Number of true predictors


# Create data
nrep = 100
w_mat = matrix(0, p, nrep)
rho_vec = seq(0,.8,.2)
pow_kf = fdr_kf = matrix(0,nrep, length(rho_vec))
pow_gm = fdr_gm= matrix(0,nrep, length(rho_vec))
pow_bh = fdr_bh = matrix(0, nrep, length(rho_vec)) 
amplitude=20


# Iterate through nreps
for(a in 1:length(rho_vec)){
for(r in 1: nrep){
    print("At rep")
    print(r)
    times <- 1:p
    sigma <- 1
    rho = rho_vec[a]
    print("With rho = ")
    print(rho)
    V= matrix(rho, p,p)
    diag(V)= 1
    x = mvrnorm(n, mu = rep(0,p), Sigma=V)
    beta_q = rnorm(q, 0, amplitude)/sqrt(n)
    true_index= sample(p,q)
    beta = numeric(p)
    beta[true_index]=beta_q
    mu = x %*%  beta
    y = mu + rnorm(n)




##########################################################
## Model-X knockoffs, using MVR loss for S-matrix.
##########################################################
    true_index = which(beta!=0)
    my_lasso_stat = function(...) stat.glmnet_coefdiff(..., nlambda=100,cores=10)
    diag_s = (1-rho)*diag(p) # This is the MVR construction asymptotically
    knockoffs = function(X) create.gaussian(X, mu=rep(0,p), Sigma=V, diag_s=diag_s)
    result = knockoff.filter(x, y, knockoffs=knockoffs, statistic = my_lasso_stat)
    if(length(result$selected) == 0){
            pow_kf[r,a]=0
            fdr_kf[r,a]=0
    }else{
    pow_kf[r,a] = length(intersect(result$selected,true_index))/length(true_index)
    fdr_kf[r,a] = 1-(length(intersect(result$selected,true_index)))/length(result$selected)
    }


cat(r,"knockoff ", pow_kf[r,a], fdr_kf[r,a],'\n')


}#end for rep
}#end for amplitude



