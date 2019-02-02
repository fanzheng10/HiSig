library(glmnet)
library(Matrix)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
xfname <- args[1]
xfname_cv <- args[2]
realyfname <- args[3]
outfname <- args[4]
cv <-as.integer(args[5])

# xfname <- 'fuse29_conn_sp.txt'
# xfname_cv <- 'fuse29_conn_gt5_sp.txt'
# realyfname <- 'pancan_mut.y-loge20.txt'
# outfname <- 'fuse29.pancan_mut.glmnet-cv5-5'
# cv <-5


X = as.matrix(read.table(xfname, header=F))
X_sp = sparseMatrix(X[,1], X[,2], index1 = F)
X_cv = as.matrix(read.table(xfname_cv, header=F))
X_cv_sp = sparseMatrix(X_cv[,1], X_cv[,2], index1 = F)

realy <- as.vector(read.table(realyfname, header = F))$V1
fit <- glmnet(X_sp, realy, 
              lambda.min = 0.0001, nlambda = 500,
              standardize=F, lower.limit=0)
lambda_arr <- fit$lambda

fitcv <- cv.glmnet(X_cv_sp, realy,
                   lambda = lambda_arr,
                   lambda.min = 0.0001, nlambda = 500,
                   standardize=F, lower.limit=0,
                   nfolds=cv)
coef = round(fit$beta, digits=6)
# ind = which(lambda_arr == fitcv$lambda.min)
ind = which.min(abs(lambda_arr - fitcv$lambda.min/10))
colnames(coef) = lambda_arr
col_best = as.matrix(coef[, ind])
colnames(col_best) =  round(fitcv$lambda.min, 9)
# write.table(col_best, file=paste(outfname, ".best.lambda.txt", sep=""), sep="\t", col.names = T)

# include p-value calculation (use the pre-determined lambda)
n_cores=7
fitrandom <- function(k) {
  y = sample(realy)
  fitr <- glmnet(X_sp, y, lambda = c(fitcv$lambda.min/10), standardize=F, lower.limit=0)
  return(as.matrix(fitr$beta))
}
nrand=10000
randbeta_all =  mclapply(1:nrand, fitrandom, mc.cores=n_cores, mc.cleanup=TRUE)
beta_df = cbind(col_best, as.data.frame(randbeta_all))
write.table(round(beta_df, digits=6), file=paste(outfname, ".best.lambda.wrand.txt", sep=""), sep="\t", col.names = F)

save(fit, fitcv, file = paste(outfname, ".RData", sep=""))




