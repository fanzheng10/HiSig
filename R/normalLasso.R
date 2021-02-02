library(glmnet)

X <- as.matrix(read.table('sim_conn.txt', header=F))
X_sp = sparseMatrix(X[,1], X[,2], index1 = F)
realy <- as.matrix(read.table('genescore.tsv', header=F))

terms = read.table('terms.txt')
cvfit <- cv.glmnet(x = X_sp, y = realy, alpha=1, lower.limit=0, nfolds=5)
fit <- cvfit$glmnet.fit

var <- which(fit$beta[,which(cvfit$lambda == cvfit$lambda.min)]>0)
var <- as.vector(var[var-length(realy)>0]-length(realy))

sel_terms = terms[var,]

writeLines(as.character(sel_terms), 'normal_lasso.txt')
