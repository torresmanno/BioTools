#' download_genomes
#'
#' @param df Assemblie table of NCBI with header
#' @param type "genome", "protein", "cds" or "feature" to download
#' @param outpath path where the files will be stored
#'
#' @importFrom R.utils gunzip
#' @return
#' @export
#'
#' @examples
#'

download_genomes <- function(df, type = "genome", outpath){

  dir.create(outpath, showWarnings = F, recursive = T)
  url.base <- df$ftp_path
  url.file <- strsplit(as.character(url.base), "/")
  url.file <- sapply(url.file, '[[', 10)
  if(type=="genome"){
    Genome.file <- paste0(outpath,"/", df$strain, "_genomic.fna.gz")
    url.genome <- paste0(url.base, "/", url.file, "_genomic.fna.gz")
  }else if(type =="protein"){
    Genome.file <- paste0(outpath, "/", df$strain, "_protein.faa.gz")
    url.genome <- paste0(url.base, "/", url.file, "_protein.faa.gz")
  }else if(type == "rna"){
    Genome.file <-paste0(outpath, "/", df$strain, "_rna_from_genomic.fna.gz")
    url.genome <-paste0(url.base, "/", url.file, "_rna_from_genomic.fna.gz")
  }else if(type == "feature"){
    Genome.file <-paste0(outpath, "/", df$strain, "_feature_table.txt.gz")
    url.genome <-paste0(url.base, "/", url.file, "_feature_table.txt.gz")
  }
  for (i in seq_along(url.file)){
    if(!file.exists(gsub(".gz", "",Genome.file[i]))){
      download.file(url=url.genome[i], destfile = Genome.file[i])
      gunzip(Genome.file[i])
    }
  }
}

download_genomes2 <- function (df, type = "genome", outpath)
{
  dir.create(outpath, showWarnings = F, recursive = T)
  url.base <- df$ftp_path
  url.file <- strsplit(as.character(url.base), "/")
  url.file <- sapply(url.file, "[[", 10)
  if (type == "genome") {
    Genome.file <- paste0(outpath, "/", df$strain, "_genomic.fna.gz")
    url.genome <- paste0(url.base, "/", url.file, "_genomic.fna.gz")
  }
  else if (type == "protein") {
    Genome.file <- paste0(outpath, "/", df$strain, "_protein.faa.gz")
    url.genome <- paste0(url.base, "/", url.file, "_protein.faa.gz")
  }
  else if (type == "cds") {
    Genome.file <- paste0(outpath, "/", df$strain, "_cds_from_genomic.fna.gz")
    url.genome <- paste0(url.base, "/", url.file, "_cds_from_genomic.fna.gz")
  }
  else if (type == "feature") {
    Genome.file <- paste0(outpath, "/", df$strain, "_feature_table.txt.gz")
    url.genome <- paste0(url.base, "/", url.file, "_feature_table.txt.gz")
  }
  for (i in seq_along(url.file)) {
    if(!file.exists(gsub(".gz", "",Genome.file[i]))){
      tryCatch({
        download.file(url=url.genome[i], destfile = Genome.file[i])
        R.utils::gunzip(Genome.file[i])
      },
      error = function(e) try({
        if(grepl("GCF", url.genome[i])){
          download.file(url=gsub("GCF", "GCA", url.genome[i]), destfile = Genome.file[i])
        }else {
          download.file(url=gsub("GCA", "GCF", url.genome[i]), destfile = Genome.file[i])
        }
        R.utils::gunzip(Genome.file[i])
      })
      )

    }
  }
}

get_strain_name <- function(df){
  df$strain <- gsub("strain=", "", df$infraspecific_name)
  df$strain[is.na(df$strain)] <- ""
  no.strain.logic <- (df$strain == ""|df$strain == "n/a" |df$strain == "type strain: N")
  if(sum(no.strain.logic) !=0){
    df$strain[no.strain.logic] <- df$isolate[no.strain.logic]
    df$strain[is.na(df$strain)] <- ""
  }
  no.strain.logic <- df$strain == ""
  if(sum(no.strain.logic) !=0){
    df$strain[no.strain.logic]  <- df[[1]][no.strain.logic]
  }
  remove_sym <-function(s){
    s <- gsub(" ", "_",s)
    s <- gsub("#", "",s)
    s <- gsub("/", "_",s)
    s <- gsub(";.*", "",s)
    s <- gsub('\\)', "",s)
    s <- gsub('\\(', "_",s)
    s <- gsub(':', "_",s)
    s <- gsub('!', "",s)
    s <- gsub('\\.', "_",s)
    s <- gsub('-', "_", s)
    s <- gsub('=', "_", s)
    s <- gsub('__', "_", s)
    s <- gsub('__', "_", s)
    s <- gsub('__', "_", s)

  }
  df$strain <- remove_sym(df$strain)
  return(df)
}


remove.dup <- function(df){
  dup.df <- df[df$strain %in% df$strain[duplicated(df$strain)],]
  # tiene en cuenta si hay duplicados dentro del mismo genero
  dup.l <- lapply(unique(dup.df$Genera), function(g) dup.df[dup.df$Genera == g,])
  dup.df2 <- do.call(rbind, lapply(dup.l, function(d) d[d$strain %in% d$strain[duplicated(d$strain)],]))

  dup <- sapply(unique(dup.df2$strain), function(s) dup.df2[dup.df2$strain == s,],USE.NAMES = T, simplify = F)

  best.candidate <- sapply(dup, function(d){
    b.c <- d[[1]][d$refseq_category != "na"]
    if(length(b.c)==0){
      b.c <- d[[1]][d$assembly_level == "Complete Genome"]
    }
    if(length(b.c) == 0){
      b.c <- d[[1]][d$assembly_level == "Chromosome"]
    }
    if(length(b.c) == 0){
      b.c <- d[[1]][d$assembly_level == "Scaffold"]
    }
    if(length(b.c) != 1){
      d <- d[order(d$seq_rel_date, decreasing = T),]
      b.c <- d[[1]][1]
    }
    return(b.c)
  })

  dup.str <- dup.df2[[1]]
  str.to.remove <- dup.str[!dup.str %in% best.candidate]

  return(df[!df[[1]] %in% str.to.remove, ])
}
