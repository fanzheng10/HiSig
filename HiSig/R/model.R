#' Load data for HiSig
#'
#' This function loads the design matrix and a response vector, the components
#' of a HiSig model. The rows and columns of the binary design matrix represent
#' genes and gene sets, respectively. If the element at (i, j) is 1, then the
#' i-th gene is a member of the j-th gene set
#' @param xfname A two column numerical file representing a sparse matrix. The
#'   1st/2nd dimension encodes the rows/columns of the binary design matrix.
#' @param yfname A numerical file containing the observed signal (response
#'   variable) for genes.
#' @param genes A vector of genes
#' @param terms A vector of gene sets
#' @param gene.as.term if the design matrix sees individual gene as a gene set, set it to TRUE to correctly compute statistics. Default is FALSE.
#' @return A list containing the `design`, `response`, `genes`, `terms` fields.
#' @export
load_data <- function(xfname, yfname, genes, terms, index1=T, gene.as.term=F) { #TODO: making genes and terms accept file names too

  X <- as.matrix(read.table(xfname, header=F))
  if (gene.as.term) {
    dim2 = length(genes) + length(terms)
  }
  else {
    dim2 = length(terms)
  }
  if (ncol(X)==2) {
    X_sp = sparseMatrix(i=X[,1], j=X[,2], x=rep(1, dim(X)[1]), index1 = index1, dims=c(length(genes), dim2))
  }
  else {
    X_sp = sparseMatrix(i=X[,1], j=X[,2], x=X[,3], index1 = index1, dims=c(length(genes), dim2))
  }


  realy <- as.matrix(read.table(yfname, header=F))
  data <-list("design"=X_sp, "response"=realy, "genes"=genes, "terms"=terms)
  return(data)
}


#' Fit a linear model with glmnet.
#'
#' This function takes a "data" object (containing `design` and `response`), and
#' fit an entire lasso path. It generates a coefficient matrix of each gene set
#' at multiple sampled strengths of regularization (`lambda`).
#' @param data A named list containing `design` and `response`.
#' @param dummy A dummy variable for parallelzation. Default is 1.
#' @param random If TRUE, the fit is for the null model. Instead of `lambda.min`
#'   and `nlambda`, need to provide a single `lambda` value.
#'   Default is FALSE.
#' @param lambda.min A value for `lambda.min.ratio` variable of the `glmnet`
#'   function.
#' @param nlambda A value for `nlambda` variable of the `glmnet` function.
#' @param pos.only If TRUE, force all coefficients to be non-negative. Default
#'   is FALSE.
#' @param lambda Provide a vector of lambda when `random` is TRUE.
#' @param simple.output Control the type of output. Default is FALSE.
#' @return When `simple.output` is TRUE, only return the `beta_max` vector (gene set
#'   impact on responses); otherwise, return more details of the
#'   fit.
#' @export
hisig_fit <- function(data,
                      dummy = 1,
                      lambda.min = 0.0001,
                      nlambda = 100,
                      pos.only = F,
                      lambda = NULL,
                      random = F,
                      simple.output=F) {
  if (pos.only == F) {
    lower.limits = -Inf
  }
  else {
    lower.limits = 0
  }
  if (ncol(data$response) == 1) {
    response <- data$response
  }
  else {
    response <- as.matrix(data$response[,dummy])
  }
  if (random) { # TODO: this is not a very good behavior. as random and simple.output are coupled in the current use case
    stopifnot(is.null(lambda) == F)
    response <- sample(response)
    fit <- glmnet(data$design, response, lower.limits = lower.limits,
                  lambda = lambda,
                  standardize=F, family='gaussian')
  }
  else {
    fit <- glmnet(data$design, response, lower.limits = lower.limits,
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

  if (simple.output) {
    return(beta_max)
  }
  else {
    out = list("beta.max"=beta_max, "coef"= coef1, "coef.norm" = norm_beta, "lambda" = fit$lambda)
    return(out)
  }
}

#' Create the HiSig null model with glmnet.
#'
#' @param data A named list containing `design` and `response`.
#' @param lambda An vector of lambda determined by the main fit. Argument passed to `hisig_fit`.
#' @param batch Parameter of parallelization. The total number of permutation is `batch*batch_size`.
#' @param batch_size Parameter of parallelization.
#' @param pos.only whether restrict to non-negative coefficient, Argument passed to `hisig_fit`.
#' @param n_cores Number of cores to use for parallelization.
#' @return A matrix to quantify the impact of gene sets under the null model.
#' @export
hisig_fit_rand <- function(data, lambda, batch=10, batch_size=10, n_cores=detectCores()-1,
                           pos.only=F) {
  for (i in 1:batch) {
    beta_max_all <- mclapply(1:batch_size, hisig_fit, data = data,
                             lambda=lambda, pos.only = pos.only, random=T,simple.output=T,
                             mc.cores=n_cores, mc.cleanup=TRUE)
    if (i == 1) {
      beta_max_combined = as.data.frame(beta_max_all)
    }
    else {
      beta_max_combined = cbind(beta_max_combined, as.data.frame(beta_max_all)) #TODO: this can be changed to apply
    }
  }
  names(beta_max_combined) = 1:(batch*batch_size)
  return(beta_max_combined)
}

#' Create the HiSig model (multi-sample mode) with glmnet. Can be used to process both real and permuted data.
#'
#' @param data A named list containing `design` and `response`.
#' @param lambda.min See `hisig_fit`.
#' @param lambda See `hisig_fit`.
#' @param pos.only See `hisig_fit`.
#' @param batch In this function, if `random=F`, the number of batches is determined by the sample size. Otheriwse, same as `hisig_fit_rand`.
#' @param batch_size See `hisig_fit_rand`.
#' @param n_cores See `hisig_fit_rand`.
#' @param random If true, shuffle the input response vector
#' @return A matrix to quantify the impact of gene sets.
# @export this function is not ready yet
hisig_fit_ms <- function(data,
                         lambda.min = 0.0001,
                         nlambda = 100,
                         pos.only = F,
                         lambda = NULL,
                         random=F,
                         # shuffle.row=F, # default is to shuffle column (i.e. shuffle the sample label)
                         batch=10, # only effective for generating null model
                         batch_size=10, n_cores=detectCores()-1) {
  nsample = ncol(data$response)
  data_s <- data
  data_rand <- data
  if (random==F) {
    batch = ceiling(nsample/batch_size)
  }
  for (i in 1:batch) {

    if (random) {
      data_rand$response = t(apply(data$response, 1, sample))
      data_s$response <- data_rand$response[,sample(ncol(data_rand$response), batch_size)]
    }
    else {
      data_s$response = data$response[,((i-1)*batch_size+1):min(i*batch_size, nsample)]
    }
    beta_max_all <- mclapply(1:batch_size, hisig_fit, data = data_s,
                             lambda.min=lambda.min, nlambda=nlambda, pos.only = pos.only, simple.output=T,
                             mc.cores=n_cores, mc.cleanup=TRUE)
    if (i == 1) {
      beta_max_combined = as.data.frame(beta_max_all)
    }
    else {
      beta_max_combined = cbind(beta_max_combined, as.data.frame(beta_max_all)) #TODO: this can be changed to apply
    }
  }
  names(beta_max_combined) = 1:(batch*batch_size)
  return(beta_max_combined)
}
