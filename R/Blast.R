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

# make.BLAST --------------------------------------------------------------

make.BLAST <- function(path="./DB/", queryfile="markers", BEST=T, blast.type="blastn", evalue="1e-5", outpath="./BLAST/NRPS/TablaI/", pattern = "_representative.db.nsq|_reference.db.nsq|_type.db.nsq") {
  if(BEST==T){
    targets <- " -max_target_seqs 1"
  }else {
    targets <- " "
  }

  ## create output folder
  dir.create(outpath, showWarnings = F, recursive = T)
  # prepare data info
  db.list <- list.files(path, pattern = pattern)  # Look for blastdb
  db.infile <- paste(path, unlist(sapply(db.list, strsplit, split = "*.nsq")), sep = "") # define db to blast
  outnames.blast <- paste(outpath, unlist(sapply(db.list, strsplit, split = "*.db.nsq")), # define name of outfile
                          ".blast.tsv", sep = "")
  ## create commands function
  cmdCreate.blast <- function(db, queryfile, outfile){
    paste(blast.type, " -query ", queryfile, " -evalue ", evalue," -db ",  db,  " -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp'", targets, " -num_threads 4 -out ",
          outfile, sep = "")
  }
  ## create actual commands
  cmds.blast <- mapply(FUN = cmdCreate.blast, db = db.infile, queryfile, outfile = outnames.blast)
  ## run commands
  sapply(cmds.blast, system)
  ##move your data to ./Retrive.Blast Directory
}


blast_bu <- function(query.fasta, strains, db.path, output, blast.type = "blastn", pattern = "_genomic.db", outfmt = 6, evalue = "1e-5"){
  if (outfmt == 6){
    out.fmt = "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp'"
  }else {
    out.fmt = outfmt
  }
  cmd.blast <- function(){
    paste0(blast.type, " -query ", query.fasta," -evalue ", evalue, " -db ",db.path, strains, pattern," -outfmt ", out.fmt, " -num_threads 8 -max_target_seqs 1 -out ", output, strains, ".tsv")
  }
  cmds <- cmd.blast()
  dir.create(output, showWarnings = F, recursive = T)

  null <- sapply(cmds, system)
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
