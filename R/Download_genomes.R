#' download_genomes
#'
#' @param df Assemblie table of NCBI with header
#' @param type "genome", "protein", "cds" or "feature" to download
#' @param outpath path where the files will be stored
#'
#' @return
#' @export
#'
#' @examples
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
  }else if(type == "cds"){
    Genome.file <-paste0(outpath, "/", df$strain, "_cds_from_genomic.fna.gz")
    url.genome <-paste0(url.base, "/", url.file, "_cds_from_genomic.fna.gz")
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
