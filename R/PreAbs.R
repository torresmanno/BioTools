#' get.pre.abs
#'
#' function to obtain the matrix of present and absent of a list of genes.
#'
#' @param ortho.tables blast tables with the orthologuos hits.
#' @param merge.by string
#'
#' @return a data frame with the present of absent of genes.
#' @export
#'
#' @examples get.pre.abs(ortho.tables = ortho.table, merge.by = "Gene")
get.pre.abs <- function(ortho.tables, merge.by = "Gene"){
  column.names <- names(ortho.tables)
  outlist <- list()
  for (i in seq_along(ortho.tables)){
    outlist[[i]] <- ortho.tables[[i]][,c(1,2)]
    names(outlist[[i]]) = c(merge.by, column.names[i])
  }
  merged.data.frame <- Reduce(function(...) merge(..., by= merge.by, all = T), outlist)

  return(merged.data.frame)
}
