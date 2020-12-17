#!/usr/bin/env Rscript
devtools::install_github("fanzheng10/HiSig", subdir="HiSig")
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
data <- load_data(args$x, args$y, genes, terms, index1 = T)

beta_sample <- hisig_fit_ms(data, nlambda = args$nlambda)
beta_null <- hisig_fit_ms(data, random = T, batch = args$permute/10, nlambda = args$nlambda)

result_ms = parse_hisig_ms(beta_sample, beta_null)

row.names(result_ms[[1]]) = terms
row.names(result_ms[[2]]) = terms

result_ms[[1]] = round(result_ms[[1]],4)
result_ms[[2]] = round(result_ms[[2]],4)

write.table(result_ms[[1]], file=paste0(args$o, '_nes.txt'), quote=F, sep='\t', col.names = F)
write.table(result_ms[[2]], file=paste0(args$o, '_q.txt'), quote=F, sep='\t', col.names = F )
