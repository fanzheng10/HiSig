# a quick implementation of GSEA preranked test
library(Matrix)
library(dplyr)

f_g = 'genes.txt'
f_t = 'terms.txt'
f_gs = 'genescore.tsv'
f_conn = 'sim_conn.txt'
eps_std = 0.001
perm = 1000

out = 'gsea_results.txt'

# probably need multi-core too

g_list = read.table(f_g, header=F)
t_list = read.table(f_t, header=F)
g_score = read.table(f_gs, header=F)
eps = rnorm(nrow(g_list), sd=eps_std)
eps = eps - min(eps)

C <- as.matrix(read.table(f_conn, header=F))
C_sp = sparseMatrix(C[,1], C[,2], index1 = F)
C_sp = C_sp[, (nrow(g_list) + 1):dim(C_sp)[2]]

g_score = g_score + eps

g_score$gene = g_list$V1
g_score_sort <- g_score[order(g_score$V1, decreasing = T),]

cal.ks <- function(i, S) {
  r = g_score_sort
  rownames(r) = r$gene
  NR = sum(r[S,]$V1) # here assume input are non-negative
  r_edge = r[1:i,]
  r_hit = r_edge[r_edge$gene %in% S, ]
  n_hit = dim(r_hit)[1]
  n_miss = i - n_hit
  rj = sum(r_hit$V1)
  Phit = rj/NR
  N = dim(r)[1]
  NH = length(S)
  Pmiss = n_miss/(N-NH)
  return(Phit - Pmiss)
}


calES <- function(x) { # x is a set of gene names
  ES = mapply(cal.ks, 1:dim(g_list)[1], rep(list(x), dim(g_list)[1]))
  # plot(1:length(ES), ES, type='l')
  return(max(ES))
}

calESnull <- function(n) {
  x = sample(g_score$gene, n)
  ES = mapply(cal.ks, 1:dim(g_list)[1], rep(list(x), dim(g_list)[1]))
  return(max(ES))
}


tsize_uniq = unique(Matrix::colSums(C_sp)) # create null distribution for these term sizes

# add parallel here
null_ES = lapply(1:length(tsize_uniq), function(i) sapply(rep(tsize_uniq[i], perm), calESnull))

cal.p <- function(i) {
  s = g_list[C_sp[,i],]
  ES = calES(s)
  size = length(s)
  size_ind = which(tsize_uniq == size)
  print(size)
  print(size_ind)
  ESnull <- null_ES[[size_ind]]
  return(mean(ESnull > ES))
}

# gset <- g_list[C_sp[,1],] # this is test
# cal.p(gset)

all.p = sapply(1:dim(C_sp)[2], cal.p)
t_list$p = all.p

write.table(t_list, file=out, quote=F, sep="\t")
