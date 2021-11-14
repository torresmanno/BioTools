#' Title ReadBlastTables
#' Function to convert the tables made by blast in an R objet
#'
#' @param files files path where the table is stored
#'
#' @return a list with the data frames
#' @export
#'
#' @examples ReadBlastTables("FZB42_168_blast.tsv")
#'
#' # To read multple files with a pattern
#'
#' strains <- c("168", "DSM7", "Ames", "ATCC_14579")
#' sapply(strains, function(s){
#'     ReadBlastTables(paste0("Blast/FZB42_", s, "_blast.tsv"))
#' }, simplify = F, USE.NAMES = T)
#'

ReadBlastTables <- function(files){
    t <- try(read.table( files,
                         header = F,
                         fill = T,
                         stringsAsFactors = F),silent = T)
    if(inherits(t, "try-error")){
      return(NULL)
    } else{
      return(t)
    }
}

#' ReadBlastTables.inpath
#'
#' Convert the tables made by blast stored in a directory in an R object
#'
#' @param inpath path where the tables are stored
#' @param pattern patter of the tables files
#'
#' @return a list with the data frames
#' @export
#'
#' @examples ReadBlastTables.inpath(inpath = "raw.bt/", pattern = ".tsv")
ReadBlastTables.inpath <- function(inpath,pattern){
  strains <- gsub(pattern = pattern, "", list.files(inpath, pattern = pattern))
  paths <- list.files(inpath, pattern = pattern ,full.names = T)
  blast.table.list <- sapply(paths, function(s){

    t <- try(read.table(s,
                        header = F,
                        fill = T,
                        stringsAsFactors = F), silent = T)
    if(inherits(t, "try-error")){
      return(NULL)
    } else{
      return(t)
    }
  }, simplify = F, USE.NAMES = T)
  names(blast.table.list) <- strains

  return(blast.table.list)
}


ReadBlastTables.batch <- function(strain, inpath, pattern){
  blast.table.list <- sapply(strain, function(s){
    t <- try(read.table( paste0(inpath, s, pattern),
                         header = F,
                         fill = T,
                         stringsAsFactors = F),silent = T)
    if(inherits(t, "try-error")){
      return(NULL)
    } else{
      return(t)
    }
  }, simplify = F, USE.NAMES = T)

  return(blast.table.list)
}
