% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{hisig_fit_ms}
\alias{hisig_fit_ms}
\title{Create the HiSig model (multi-sample mode) with glmnet. Can be used to process both real and permuted data.}
\usage{
hisig_fit_ms(
  data,
  lambda.min = 0.0001,
  nlambda = 100,
  pos.only = F,
  lambda = NULL,
  random = F,
  batch = 10,
  batch_size = 10,
  n_cores = detectCores() - 1
)
}
\arguments{
\item{data}{A named list containing \code{design} and \code{response}.}

\item{lambda.min}{See \code{hisig_fit}.}

\item{pos.only}{See \code{hisig_fit}.}

\item{lambda}{See \code{hisig_fit}.}

\item{random}{If true, shuffle the input response vector}

\item{batch}{In this function, if \code{random=F}, the number of batches is determined by the sample size. Otheriwse, same as \code{hisig_fit_rand}.}

\item{batch_size}{See \code{hisig_fit_rand}.}

\item{n_cores}{See \code{hisig_fit_rand}.}
}
\value{
A matrix to quantify the impact of gene sets.
}
\description{
Create the HiSig model (multi-sample mode) with glmnet. Can be used to process both real and permuted data.
}
