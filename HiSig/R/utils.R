#' convert the regulon from sparse matrix format to the Aracne/Viper format
#' @param data The object created by the 'load_data' function
#' @param design The output file name
#' @return regulon A regulon object in viper package
#' @export
sparse2aracne <- function(data) {
  X_sp <- which(as.matrix(data$design != 0), arr.ind = T)
  regulon = list()
  for (i in 1:length(data$terms)) {
    target_ids = X_sp[X_sp[,2] == i,1]
    targets = data$genes[target_ids]

    spm_ids = which(X_sp[,2] ==i) # ids in the sparse matrix
    mor = data$design@x[spm_ids]

    regul_element = list()
    regul_element$tfmode = mor
    names(regul_element$tfmode) = targets
    regul_element$likelihood = rep(1, length(mor))
    # line = paste0(data$terms[i], paste0('\t', targets, '\t', mor, collapse=''))
    # if (i ==1) {
    #   write(line, file=out)
    # }
    # else {
    #   write(line, file=out, append=T)
    # }
    regulon[[i]] <- regul_element
  }
  names(regulon) <- data$terms
  return(regulon)
}
