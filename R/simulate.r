library(fabricatr)

nlayers = 4
M = c(5,5,5,5)
layer_names = paste0('L', 1:4)
p = 0.4

hier_data <-
  fabricate(
    L = add_level(
      N = M[[1]],
      select = draw_binary(prob=p, N=N)
      # weight = rpois(N, 1) * select,
      # sum_weight = weight
    )
  )
names(hier_data)[names(hier_data) == 'L'] <- 'L1'
names(hier_data)[names(hier_data) == 'select'] <- 'L1.select'
# names(hier_data)[names(hier_data) == 'weight'] <- 'L1.weight'
hier_data$sum_weight = 0

for (x in 2:nlayers) {
  hier_data <-
    fabricate(
      data = hier_data,
      L = add_level(
        N = M[[x]],
        select = draw_binary(prob=p, N=N) # make new selection
      )
    )
  names(hier_data)[names(hier_data) == 'L'] <- paste0('L', x)
  names(hier_data)[names(hier_data) == 'select'] <- paste0('L', x,'.select')
  # names(hier_data)[names(hier_data) == 'weight'] <-paste0('L', x, '.weight')
  
}

# add weights to selected 
hier_data$weight = 0
N = nrow(hier_data)
for (x in 1:nlayers) {
  cols = paste0('L', x, '.select')
  weight = rpois(N, 1)* hier_data[[cols]]
  hier_data$weight <- hier_data$weight + weight
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
  X <- paste0(hier_data$sum_weight)
  outf = 'genescore.tsv'
  filecon <- file(outf)
  writeLines(X, filecon)
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

