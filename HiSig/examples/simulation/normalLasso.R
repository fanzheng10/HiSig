library(glmnet)

setwd("~/Desktop/UCSD_Research/2021/hisig_paper/N3_m5_p0.3_1")

realy <- as.matrix(read.table('genescore.tsv', header=F))

# the mode with identity

X <- as.matrix(read.table('sim_conn.txt', header=F))
X_sp = sparseMatrix(X[,1], X[,2], index1 = F)

terms = read.table('terms.txt')
cvfit <- cv.glmnet(x = X_sp, y = realy, alpha=1, lower.limit=0, nfolds=5)
fit <- cvfit$glmnet.fit

var <- which(fit$beta[,which(cvfit$lambda == cvfit$lambda.min)]>0)
var <- as.vector(var[var-length(realy)>0])
terms$coef = 0
terms[var-length(realy),]$coef =  fit$beta[,which(cvfit$lambda == cvfit$lambda.min)][var]
terms$coef = format(terms$coef, digits=5)
write.table(terms, 'normal_lasso_w_identity.txt', quote=F, sep='\t')

# the mode w/o identity
X[,2] = X[,2] - max(X[,1]) -1
X = X[(length(realy) +1):dim(X)[1], ]

X_sp = sparseMatrix(X[,1], X[,2], index1 = F)

terms = read.table('terms.txt')
cvfit <- cv.glmnet(x = X_sp, y = realy, alpha=1, lower.limit=0, nfolds=5)
fit <- cvfit$glmnet.fit

var <- which(fit$beta[,which(cvfit$lambda == cvfit$lambda.min)]>0)
terms$coef = 0
terms[var,]$coef =  fit$beta[,which(cvfit$lambda == cvfit$lambda.min)][var]
terms$coef = format(terms$coef, digits=5)
write.table(terms, 'normal_lasso_no_identity.txt', quote=F, sep='\t')
