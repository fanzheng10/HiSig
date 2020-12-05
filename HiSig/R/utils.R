#' convert the regulon from sparse matrix format to the Aracne/Viper format
#' @param data The object created by the 'load_data' function
#' @param design The output file name
#' @return NA
#' @export
sparse2aracne <- function(data, out) {
  X_sp <- which(as.matrix(data$design != 0), arr.ind = T)
  for (i in 1:length(data$terms)) {
    target_ids = X_sp[X_sp[,2] == i,1]
    targets = data$genes[target_ids]

    spm_ids = which(X_sp[,2] ==i) # ids in the sparse matrix
    mor = data$design@x[spm_ids]

    line = paste0(data$terms[i], paste0('\t', targets, '\t', mor, collapse=''))
    if (i ==1) {
      write(line, file=out)
    }
    else {
      write(line, file=out, append=T)
    }
  }
}
