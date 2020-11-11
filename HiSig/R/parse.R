#' Load a Matrix
#'
#' @export

getterminfo <- function(i, design, array) {
  info = array[which(design[, i])]
  return(info)
}


parse_hisig <- function(data, impact, term.names, gene.names, signal2=NA, gene.as.term=F, showedge=20) {

  data$response = round(data$response, 2)
  ngenes = length(gene.names)
  nterms = length(term.names)

  stopifnot(ngenes==length(data$response))
  term.sizes = colSums(data$design)

  # calculate p.values

  impact.main = impact[,1]
  impact.rand = impact[,2:dim(impact)[2]]
  pval <- rowSums(impact.rand > impact.main)/dim(impact)[1]
  qval = p.adjust(pval, method='BH')

  df = cbind(impact.main, pval, qval)
  df = as.data.frame(df)
  names(df) = c("Selection.pressure", "p", "q")

  if (gene.as.term) {
    df = df[(ngenes+1):dim(impact)[2], ]
  }
  df$System.names = term.names
  df$System.sizes = term.sizes
  nsys = dim(df)[1]

  # list the genes in the systems
  allranks = rank(-abs(data$response), ties.method = 'max')

  if (gene.as.term) {
    genesets = lapply(1:nsys, getterminfo, data$design[,(ngenes+1):dim(impact)[1]], gene.names)
    generanks = lapply(1:nsys, getterminfo, data$design[,(ngenes+1):dim(impact)[1]], allranks)
    signal1 = lapply(1:nsys, getterminfo, data$design[,(ngenes+1):dim(impact)[1]], data$response)
  }
  else {
    genesets = lapply(1:nsys, getterminfo, data$design, gene.names)
    generanks = lapply(1:nsys, getterminfo, data$design, allranks)
    signal1 = lapply(1:nsys, getterminfo, data$design, data$response)
  }
  gene_sort_ids = lapply(signal1, function(x) order(abs(x), decreasing=T)) # all sorted by signal1

  genesets_sorted = lapply(1:nsys, function(x) genesets[[x]][gene_sort_ids[[x]]] )
  generanks_sorted = lapply(1:nsys, function(x) generanks[[x]][gene_sort_ids[[x]]] )
  signal1_sorted = lapply(1:nsys, function(x) signal1[[x]][gene_sort_ids[[x]]] )

  genesets_sorted = as.character(lapply(genesets_sorted, function(x) paste(x[1:showedge], collapse = '|')))
  generanks_sorted = as.character(lapply(generanks_sorted, function(x) paste(x[1:showedge], collapse = '|')))
  signal1_sorted = as.character(lapply(signal1_sorted, function(x) paste(x[1:showedge], collapse='|')))

  df$Genes = genesets_sorted
  df$Gene.ranks = generanks_sorted
  df$Signal = signal1_sorted

  if (is.na(signal2) ==F) {
    stopifnot(ngenes=length(signal2))
    if (gene.as.term) {
      signal2 = lapply(1:nsys, getterminfo, data = data$design[,(ngenes+1):dim(impact)[1]], array=signal2)
    }
    else {
      signal2 = lapply(1:nsys, getterminfo, data = data$design, array = signal2)

    }
    signal2_sorted = lapply(1:nsys, function(x) signal2[[x]][gene_sort_ids[[x]]] )
    signal2_sorted = as.character(lapply(signal2_sorted, function(x) paste(x[1:showedge], collapse='|')))
    df$Signal.2 = signal2_sorted
  }
  df = as.data.frame(df)
  df %>% filter(Selection.pressure > 0)
  df %>% arrange(q, p, desc(Selection.pressure))

  return(df)
}

output <- function(df, out) {
  colorder= c("System.names", "System.sizes", "Selection.pressure", "p", "q", "Genes", "Gene.ranks", "Signal")
  if ("Signal.2" %in% names(df) ) {
    colorder = c(colorder, "Signal.2")
  }
  write.table(df[,colorder], file=out, quote=F, row.names = F, sep="\t")
}
