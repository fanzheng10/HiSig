library(fabricatr)
library(optparse)

option_list = list( # don't use R -f, but Rscript
  make_option(c("-n", "--nlayer"),
              type="integer", default=4),
  make_option(c("-m", "--size_per_layer"),
              type="integer", default=5),
  make_option(c("-p", "--prob_pos"),
                type="double", default=0.1),
  make_option(c("--log"), action='store_true'), # only applies when simulating poission distribution
  make_option(c('--gaussian'), action='store_true') # if true, will simulate Gaussian distribution (with positive and negative)
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

nlayers = opt$nlayer +1
M = rep(opt$size_per_layer, nlayers)
layer_names = paste0('L', 1:nlayers)
p = opt$prob_pos

#local test
# nlayers = 3
# M = rep(2, nlayers)
# layer_names = paste0('L', 1:nlayers)
# p = 0.2

if (opt$gaussian) {
  hier_data <-
    fabricate(
      L = add_level(
        N = M[[1]],
        select = sample(c(draw_binary(prob=p, N=N/2), -1*draw_binary(prob=p, N=N/2)))
      )
    )
} else {
  hier_data <-
    fabricate(
      L = add_level(
        N = M[[1]],
        select = draw_binary(prob=p, N=N)
      )
    )
}

names(hier_data)[names(hier_data) == 'L'] <- 'L1'
names(hier_data)[names(hier_data) == 'select'] <- 'L1.select'
# names(hier_data)[names(hier_data) == 'weight'] <- 'L1.weight'

for (x in 2:nlayers) {

  if (opt$gaussian) {
    hier_data <-
      fabricate(
        L = add_level(
          N = M[[1]],
          select = sample(c(draw_binary(prob=p, N=N/2), -1*draw_binary(prob=p, N=N/2)))
        )
      )
  } else {
    hier_data <-
      fabricate(
        L = add_level(
          N = M[[1]],
          select = draw_binary(prob=p, N=N)
        )
      )
  }
  names(hier_data)[names(hier_data) == 'L'] <- paste0('L', x)
  names(hier_data)[names(hier_data) == 'select'] <- paste0('L', x,'.select')
  # names(hier_data)[names(hier_data) == 'weight'] <-paste0('L', x, '.weight')

}

# add weights to selected
hier_data$weight = 0
N = nrow(hier_data)
for (x in 1:nlayers) {
  cols = paste0('L', x, '.select')
  if (opt$gaussian == F) {
    weight = rpois(N, 1)* hier_data[[cols]] # this is ok, Poisson process is memoryless
  }
  else {
    # flip half of the selected ones
    weight = rnorm(N, mean=2, sd=1) * hier_data[[cols]] # need to set L.select pos and neg
  }
  hier_data$weight <- hier_data$weight + weight
}

if (opt$gaussian==F && opt$log ==T) {
  hier_data$weight <- log1p(hier_data$weight)
}


write.conn<-function(nlayers) {
  # attach(hier_data)
  # cluster names
  clust_names = c()
  for (i in 1:(nlayers-1)) {
    coln = paste0('L', i)
    clust_names <- c(clust_names, paste0('L', i, '_', unique(hier_data[[coln]])))
  }

  # write conn
  outf = 'sim_conn.txt'
  X<-paste0(1:nrow(hier_data)-1, "\t", 1:nrow(hier_data)-1)
  for (i in 1:(nlayers-1)) {
    coln = paste0('L', i)
    x = match( paste0('L', i, '_', hier_data[[coln]]), clust_names) + nrow(hier_data) -1
    x = paste0(1:nrow(hier_data)-1, "\t", x)
    X <-c(X, x)
  }
  filecon <- file(outf)
  writeLines(X, filecon)
  close(filecon)

  # write terms
  outf = 'terms.txt'
  filecon <- file(outf)
  writeLines(clust_names, filecon)
  close(filecon)

  # write truth
  sel_clusts = c()
  for (i in 1:(nlayers-1)) {
    coln = paste0('L', i)
    cols = paste0('L', i, '.select')
    sel = unique(hier_data[hier_data[[cols]] == 1,][[coln]])
    if (length(sel) > 0) {
      sel = paste0(coln, '_', sel)
      sel_clusts <- c(sel_clusts, sel)
    }
  }
  outf = 'selected_terms.txt'
  filecon <- file(outf)
  writeLines(sel_clusts, filecon)
  close(filecon)

}

write.score <- function() {
  X <- hier_data$weight
  outf = 'genescore.tsv'
  filecon <- file(outf)
  writeLines(format(X, digits=3), filecon)
  close(filecon)
}

write.gene <- function() {
  X <-paste0('gene', 1:nrow(hier_data))
  outf = 'genes.txt'
  filecon <- file(outf)
  writeLines(X, filecon)
  close(filecon)
}

write.conn(nlayers)
write.score()
write.gene()

