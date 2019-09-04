#' make.DB
#'
#' Make the Blast Data Base with makeblastdb program
#'
#' @param path path where the sequence are stored
#' @param pattern pattern of the sequences files
#' @param outpath path where the database will be stored
#' @param dbtype string "prot" or "nucl"
#'
#' @return Data Bases files
#' @export
#'
#' @examples make.DB(path = "Fasta/", pattern = ".faa", outpath = "DB/", dbtype = "prot")
make.DB <- function(path="FASTA/", pattern=".faa", outpath="DB/", dbtype){

  input.genomes <- list.files(path = path, pattern = pattern)
  numfiles <- length(input.genomes)
  outnames.db <- paste(outpath, unlist(sapply(input.genomes, strsplit, split = pattern)),
                       ".db", sep = "")
  ## create commands function
  cmdCreate.db <- function(infile, outfile){
    paste("makeblastdb -in ", path, infile, " -out ",
          outfile, " -dbtype ", dbtype, sep = "")
  }
  # create output folder
  dir.create(outpath, showWarnings = F, recursive = T)
  ## create actual commands
  cmds.db <- mapply(FUN = cmdCreate.db, infile = input.genomes, outfile = outnames.db)
  ## run commands
  sapply(cmds.db, system)
}


#' blast
#'
#' run blast within R with diferent paramete
#'
#' @param type blast type to be performed. Ex blastn
#' @param query file with the fasta sequences of the query
#' @param subject file with the fasta sequence for the subject
#' @param output output file
#' @param evalue evalue cutoff
#'
#' @return make the blast searches and are stored in a directory
#' @export
#'
#' @examples blast("blastn", "FZB42_cds.fasta", "168_genomic.fasta", "FZB42_168_blast.tsv")
#'
#' # To run for more than one strain
#' strains <- c("168", "DSM7", "Ames", "ATCC_14579")
#'
#' sapply(strains, function(s) {
#'     blast(type = "blastn",
#'           query = "FZB42_cds.fasta",
#'           subject = paste0("Genomes/", s, "_genomic.fna"),
#'           output = paste0("blast/FZB42_", s, "_blast.tsv"))
#'     })

blast <- function(type, query, subject, output, evalue = "1e-5"){
  out.fmt = "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp'"
  cmd.blast <- paste0(type, " -query ", query," -evalue ", evalue, " -subject ", subject, " -outfmt ", out.fmt, " -max_target_seqs 1 -out ", output)
}


#' blast.db
#'
#' @param type blast type to be performed. Ex blastn
#' @param query file with the fasta sequences of the query
#' @param db file of the db subject
#' @param output output file
#' @param evalue evalue cutoff
#' @param threads number of threads
#'
#'
#' @return
#' @export
#'
#' @examples blast.db("blastn", "FZB42_cds.fasta", "168_genomic.db", "FZB42_168_blast.tsv")
#'
#' # To run for more than one strain
#' strains <- c("168", "DSM7", "Ames", "ATCC_14579")
#'
#' sapply(strains, function(s) {
#'     blast.db(type = "blastn",
#'              query = "FZB42_cds.fasta",
#'              db = paste0("Genomes/", s, "_genomic.db"),
#'              output = paste0("blast/FZB42_", s, "_blast.tsv"))
#'     })
blast.db <- function(type, query, db, output, evalue = "1e-5", threads = 1){
  out.fmt = "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp'"
  cmd.blast <- paste0(type, " -query ", query," -evalue ", evalue, " -db ", db, " -outfmt ", out.fmt, " -max_target_seqs 1 -num_threads ", threads, " -out ", output)
}
