#' orthoBlast
#'
#' remove duplicated in query and subject and subset by a coverage and id.cutoff
#'
#' @param raw.blast raw blast tables obtained with ReadBlastTables functino
#' @param coverage numeric coverage cutoff
#' @param id.cutoff numeric, identity percent cutoff
#'
#' @return return a list with the data frames with the tables
#' @export
#'
#' @examples orthoBlast(raw.blast = rawblast, coverage = 50, id.cutoff = 70)
#'
#' # if the tables are in a list
#'
#' sapply(raw.blast.list, orthoBlast, coverage = 70, id.cutoff = 70, simplify = F, USE.NAMES = T)

orthoBlast <- function(raw.blast, coverage, id.cutoff){
  df.out <- raw.blast %>%
    filter(.[[13]] >= coverage) %>%
    filter(.[[3]] >= id.cutoff) %>%
    arrange(.[[1]], desc(.[[12]])) %>%
    group_by(.[[1]]) %>% slice(1) %>%
    ungroup()
  return(df.out)
}

#' clean.by.sub
#'
#' clean the blast tables by subject and coverage and identity percent cutoff
#'
#' @param df data frame to be subset
#' @param cov coverage cutoff to apply
#' @param id percent identity cutoff to apply
#'
#' @return return the data frame subseted
#' @export
#'
#' @examples clean.by.sub(df = blast.table, cov = 50, id = 70)
clean.by.sub <- function(df, cov, id){
  df <- df[df[[13]] >= cov, ]
  df <- df[df[[3]] >=id,]
  df <- df[order(df[[2]], df[[12]], decreasing = T), ]
  return(df[!duplicated(df[[2]]),])}

#' orthology.pairs
#'
#' remove duplicated in query and subject and subset by a coverage and id.cutoff but return the a table with the pair of orthologs
#'
#' @param df data freame of the blast table
#' @param id numeric, percent identity cutoff
#' @param cov numeric, coverage cutoff
#'
#' @return data frame of two columns with the orthologous pairs
#' @export
#'
#' @examples orthology.pairs(blast.table, id = 70, cov = 50)
orthology.pairs <- function(df, id, cov){
  df.out <- df %>%
    filter(.[[13]] >= cov) %>%
    filter(.[[3]] >= id) %>%
    arrange(.[[1]], desc(.[[12]])) %>%
    distinct(.[[1]], .keep_all = T)%>%
    select(c(1,2))

  return(df.out)
}

#' clean.by.query
#'
#' remove the duplicated query hits and apply a cutoff in identity percent and coverage
#'
#' @param df data frame with the blast table
#' @param id numeric, identity percent cutoff
#' @param cov numeric, coverage cutoff
#'
#' @return a data freame with the blast table cleaned
#' @export
#'
#' @examples clean.by.query(df = blast.table, id = 70, cov = 50)
clean.by.query <- function(df, id, cov){
  df <- df[df[[13]] >= cov, ]
  if(is.null(df)){
    stop()
  }
  df <- df[df[[3]] >=id,]
  df <- df[order(df[[1]], df[[12]], decreasing = T), ]
  return(df[!duplicated(df[[1]]),])}
