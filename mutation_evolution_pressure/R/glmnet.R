library(glmnet)
library(Matrix)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
xfname <- args[1] # sparse matrix
yfname <- args[2]
outfname <- args[3]
mode <- as.integer(args[4]) # if 1, calculate p-value by weight; if 2, calculate p-value by proportion
batch = 1000
if (length(args) ==5) {
  batch = as.integer(args[5])
}
batch_size = 10

n_cores = 7

#local test directory

# setwd("~/Desktop/UCSD Research/2018_quarter3/test_biglasso/understand_lasso/fuse29")
# 
# xfname <- "fuse29_conn_sp.txt" # sparse matrix
# yfname <- "pancan_mut.y-loge20.txt"
# outfname <- "fuse29.pancan_mut_loge20.multiscale"
# n_cores = 2
# mode =2
# batch = 1
# batch_size =2

X <- as.matrix(read.table(xfname, header=F))
X_sp = sparseMatrix(X[,1], X[,2], index1 = F)
realy <- as.matrix(read.table(yfname, header=F))

fit <- glmnet(X_sp, realy, 
              lambda.min = 0.0001, nlambda = 500,
              standardize=F, lower.limit=0)
coef = as.matrix(fit$beta)
coef = round(coef, digits=6)
colnames(coef) = round(fit$lambda, digits = 9)
write.table(coef, file=paste(outfname, ".coef", sep=""), sep="\t")

coef1 = coef[,colSums(coef) >0]
norm_beta = scale(coef1, center=FALSE, scale=colSums(coef1) + fit$a0[colSums(coef) >0])
beta_max = apply(norm_beta, 1, max)
rank = apply(-coef, 2, rank)
rand_max = apply(rank, 1, max)

fitrandom1 <- function(k) {
  y = sample(realy)
  fitr <- glmnet(X_sp, y, lambda = fit$lambda, 
                 standardize=F, lower.limit=0)
  rcoef <- as.matrix(round(fitr$beta, digits=6))
  return (coef <= rcoef)
  # return(as.matrix(fitr$beta))
}

fitrandom2 <- function(k) {
  y = sample(realy)
  fitr <- glmnet(X_sp, y, lambda = fit$lambda, 
                 standardize=F, lower.limits=0)
  rcoef <- as.matrix(fitr$beta)
  rcoef = round(rcoef, digits=6)
  rcoef1 = rcoef[, colSums(rcoef) >0]
  # write coef file
  # write.table(round(rcoef, digits=6), file = paste(outfname, '.coef-random', k, sep=""), sep='\t')
  
  rnorm_beta = scale(rcoef1, center=FALSE, scale=colSums(rcoef1) + fit$a0[colSums(rcoef) >0])
  rbeta_max = apply(rnorm_beta, 1, max)
  return(rbeta_max)
}

fitrandom3 <- function(k) {
  y = sample(realy)
  fitr <- glmnet(X_sp, y, lambda = fit$lambda,
                 standardize=F, lower.limits=0)
  rcoef <- as.matrix(fitr$beta)
  rcoef = round(rcoef, digits=6)
  rcoef1 = rcoef[, colSums(rcoef) >0]
  # write coef file
  # write.table(round(rcoef, digits=6), file = paste(outfname, '.coef-random', k, sep=""), sep='\t')

  rnorm_beta = scale(rcoef1, center=FALSE, scale=colSums(rcoef1) + fit$a0[colSums(rcoef) >0])
  rbeta_sum = apply(rnorm_beta, 1, sum)
  return(rbeta_sum)
}

# fitrandom3 <- function(k) {
#   y = sample(realy)
#   fitr <- glmnet(X_sp, y, lambda = fit$lambda,
#                  standardize = F, lower.limits = 0)
#   rcoef <- as.matrix(fitr$beta)
#   srank <- apply(-rcoef, 2, rank, ties.method='last')
#   srank_max <- apply(srank, 1, min)
#   return(srank_max)
# }

add <- function(x) Reduce("+", x)

# nrand=10000
if (mode == 1) {
  for (i in 1:batch) {
    randcompare_all = mclapply(1:batch_size, fitrandom1, mc.cores=n_cores, mc.cleanup=TRUE)
    if (i == 1) {
      randcompare_sum = add(randcompare_all) 
    }
    else {
      randcompare_sum = randcompare_sum + add(randcompare_all)
    }
  }
  pvals = randcompare_sum / (batch_size * batch)
  # pvals_min = apply(pvals, 1, min) # is taking best p-value fair?
  # outdf = as.data.frame(as.matrix(cbind(coef, pvals_min)))
  outdf = as.data.frame(pvals)
  
  write.table(round(outdf, digits=6), file=paste(outfname, ".weight-p.txt", sep=""), sep="\t", col.names = F)
}
if (mode == 2) {
  for (i in 1:batch) {
    beta_max_all <- mclapply(1:(batch_size), fitrandom2, mc.cores=n_cores, mc.cleanup=TRUE)
    if (i == 1) {
      beta_max_combined = cbind(beta_max, as.data.frame(beta_max_all))
    }
    else {
      beta_max_combined = cbind(beta_max_combined, as.data.frame(beta_max_all))
    }
  }
  write.table(round(beta_max_combined, digits=6), file=paste(outfname, ".impact-w-rand.txt", sep=""), sep="\t", col.names = F)
}

if (mode == 3) {
  for (i in 1:batch) {
    beta_max_all <- mclapply(1:(batch_size), fitrandom3, mc.cores=n_cores, mc.cleanup=TRUE)
    if (i == 1) {
      beta_max_combined = cbind(beta_max, as.data.frame(beta_max_all))
    }
    else {
      beta_max_combined = cbind(beta_max_combined, as.data.frame(beta_max_all))
    }
  }
  write.table(round(beta_max_combined, digits=6), file=paste(outfname, ".impact2-w-rand.txt", sep=""), sep="\t", col.names = F)
}
# if (mode == 3) {
#   for (i in 1:batch) {
#     rank_max_all = mclapply(1:(batch_size), fitrandom3, mc.cores = n_cores, mc.cleanup = TRUE)
#   }
#   if (i == 1) {
#     rank_max_combined = cbind(rank_max, as.data.frame(rank_max_all))
#   }
#   else {
#     rank_max_combined = cbind(rank_max_combined, as.data.frame(rank_max_all))
#   }
#   write.table(round(rank_max_combined), file=paste(outfname, ".rank-w-rand.txt", sep=""), sep="\t", col.names = F)
# }