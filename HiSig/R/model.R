#' Load a Matrix
#'
#' @export
load_data <- function(xfname, yfname, index1=T) {

  X <- as.matrix(read.table(xfname, header=F))
  if (index1==F) {
    X_sp = sparseMatrix(X[,1], X[,2], index1 = F)
  }
  else {
    X_sp = sparseMatrix(X[,1], X[,2], index1 = T)
  }

  realy <- as.matrix(read.table(yfname, header=F))
  data <-list("design"=X_sp, "response"=realy)
  return(data)
}

hisig_fit <- function(data,
                      dummy=1,
                    lambda.min=0.0001,
                    nlambda=100,
                    pos.only=F,
                    lambda = NULL,
                    random=F) {
  if (pos.only == F) {
    lower.limits = -Inf
  }
  else {
    lower.limits = 0
  }
  if (random) {
    stopifnot(is.null(lambda) == F)
    data$response <- sample(data$response)
    fit <- glmnet(data$design, data$response, lower.limits = lower.limits,
                  lambda = lambda,
                  standardize=F, family='gaussian')
  }
  else {
    fit <- glmnet(data$design, data$response, lower.limits = lower.limits,
                  lambda.min.ratio=lambda.min, nlambda=nlambda,
                  standardize=F, family='gaussian')
  }
  coef = as.matrix(fit$beta)
  coef = round(coef, digits=6)
  colnames(coef) = round(fit$lambda, digits = 9)

  coef = abs(coef)
  coef1 = coef[, abs(colSums(coef)) >0]
  norm_beta = scale(coef1, center=FALSE, scale=colSums(coef1) + abs(fit$a0[colSums(coef) >0]))
  beta_max = apply(norm_beta, 1, max)

  if (random) {
    return(beta_max)
  }
  else {
    out = list("beta.max"=beta_max, "coef"= coef1, "coef.norm" = norm_beta, "lambda" = fit$lambda)
    return(out)
  }
}


hisig_fit_rand <- function(data, lambda, batch=10, batch_size=10, n_cores=detectCores()-1) {
  for (i in 1:batch) {
    beta_max_all <- mclapply(1:batch_size, hisig_fit, data = data, lambda=lambda, random=T, mc.cores=n_cores, mc.cleanup=TRUE)
    if (i == 1) {
      beta_max_combined = as.data.frame(beta_max_all)
    }
    else {
      beta_max_combined = cbind(beta_max_combined, as.data.frame(beta_max_all))
    }
  }
  names(beta_max_combined) = 1:(batch*batch_size)
  return(beta_max_combined)
}
