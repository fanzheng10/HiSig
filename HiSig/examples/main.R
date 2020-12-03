#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
library(HiSig)

parser <- ArgumentParser()
parser$add_argument("--x", required=T, help="the design matrix")
parser$add_argument("--y", required=T, help="the response vector")
parser$add_argument("--g", required=T, help='the gene list')
parser$add_argument("--t", required=T, help="the system list")
parser$add_argument("--o", required=T, help="the output name")
parser$add_argument('--permute', default=1000, type="integer", help='number of permutation')
parser$add_argument('--nlambda', default=100, type="integer", help='number of lambda sampled along the lasso path')
args <- parser$parse_args()

genes = readLines(args$g)
terms = readLines(args$t)
data <- load_data(args$x, args$y, genes, terms, index1 = F)

hisig_out <- hisig_fit(data, nlambda=args$nlambda)
beta_max_rand <- hisig_fit_rand(data, hisig_out$lambda, batch=args$permute/10)

impact <- cbind(hisig_out$beta.max, beta_max_rand)

df <- parse_hisig(data, impact, term.names = terms, gene.names = genes)
output(df, paste0(args$o,'.txt'))
