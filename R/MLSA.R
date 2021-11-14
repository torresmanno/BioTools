#' get.core
#'
#' Funciton to create the matrix with the position of the genes in the core or pangeonme
#'
#' @param ortho.tables Named list of Data frame. Table with the orthologues blast hits like the obtained by blast
#' @param merge.by Character. Select the column name to be merged
#' @param all Logical. If the table should obtain the pangenome or the core.
#'
#' @return return a data frame with the position of the genes in the genomes.
#' @export
#'
#' @examples get.core(ot, merge.by = "Gene", all = F)
get.core <- function(ortho.tables, merge.by = "Gene", all = F){
  column.names <- names(ortho.tables)
  outlist <- list()
  for (i in seq_along(ortho.tables)){
    outlist[[i]] <- ortho.tables[[i]][,c(1,2,9,10)]
    names(outlist[[i]]) = c(
      merge.by,
      paste(column.names[i], "Contig" , sep = "."),
      paste(column.names[i], "Start" , sep = "."),
      paste(column.names[i], "End" , sep = ".")
    )
  }
  merged.data.frame <- Reduce(function(...) merge(..., by= merge.by, all = all), outlist)

  return(merged.data.frame)
}

#' get.gen
#'
#' Function to extract the get.core genes in a given genome.
#'
#'
#' @param Genome.path path where the fasta files are stored
#' @param Genome.pattern pattern to the fasta files
#' @param Genome Character string with the name of the genomes
#' @param HK.table data frame obtained with the get.core function
#' @param path.multifasta path where the fasta files with the genes will be stored
#'
#' @return create fasta files with the genes of the given genome
#' @export
#'
#' @examples get.gen(Genome.path = "Genome.FASTA", Genome.pattern = "_genomic.fna", Genome = "FZB42", HK.table = Core, path.multifasta = "mutif/")

get.gen <- function(Genome.path= "../Genome.FASTA", Genome.pattern = "_genomic.fna", Genome, core.table ,path.multifasta = "./Genes/"){
  library(seqinr)

  merged.blast.result <- core.table

  sec.genomes <- seqinr::read.fasta(file = paste0(Genome.path, "/", Genome, Genome.pattern), seqtype = "DNA")

  getGen <- function(Gene) {
    Gen.num <- match(Gene,merged.blast.result[[1]])
    Gene.coding <- merged.blast.result[Gen.num,paste0(Genome,".Start")] <  merged.blast.result[Gen.num,paste0(Genome,".End")]
    ##identificar nombre contig
    Gene.contig <- as.character(merged.blast.result[Gen.num,paste(Genome,".Contig", sep="")])
    ##generar lista con secuencias
    ###Generate an object with the selected contig
    object.sec  <- sec.genomes[grep(paste(Gene.contig, "_tag", sep=""),paste(names(sec.genomes), "_tag", sep=""), fixed = T )]
    ### Generate an object with the Gene sequence
    if (Gene.coding == T) {
      Gene.sec <- seqinr::getFrag(object.sec, merged.blast.result[Gen.num,paste(Genome,".Start", sep="")], merged.blast.result[Gen.num,paste(Genome,".End", sep="")])
    }    else {
      Gene.sec.rc <- seqinr::getFrag(object.sec, merged.blast.result[Gen.num,paste(Genome,".End", sep="")], merged.blast.result[Gen.num,paste(Genome,".Start", sep="")])
      Gene.sec <- rev(seqinr::comp(seqinr::getSequence(Gene.sec.rc[[1]])))
    }
    ## create folder
    dir.create(path.multifasta, showWarnings = F, recursive = T)
    ## write fasta
    seqinr::write.fasta(sequences = Gene.sec, names = Genome, file.out = paste(path.multifasta, "/", Gene, ".fasta", sep=""), open = "a", nbchar = 60, as.string = FALSE)
  }
  sapply(X = merged.blast.result[[1]], getGen)
}

get.gen.norestrictive <- function(Genome.path= "../Genome.FASTA", Genome.pattern = "_genomic.fna", Genome, HK.table ,path.multifasta = "./HK.multif"){
  library(seqinr)

  merged.blast.result <-  HK.table[,c("Gene", paste0(Genome, ".Contig"),paste0(Genome, ".Start"), paste0(Genome,".End"))]
  merged.blast.result <- merged.blast.result[!is.na(merged.blast.result[,2]),]

  sec.genomes <- read.fasta(file = paste0(Genome.path, "/", Genome, Genome.pattern), seqtype = "DNA")

  getGen <- function(Gene) {
    Gen.num <- match(Gene,merged.blast.result[[1]])
    Gene.coding <- merged.blast.result[Gen.num,paste0(Genome,".Start")] <  merged.blast.result[Gen.num,paste0(Genome,".End")]
    ##identificar nombre contig
    Gene.contig <- as.character(merged.blast.result[Gen.num,paste0(Genome,".Contig")])
    ##generar lista con secuencias
    ###Generate an object with the selected contig
    object.sec  <- sec.genomes[grep(paste(Gene.contig, "_tag", sep=""),paste(names(sec.genomes), "_tag", sep=""), fixed = T )]
    ### Generate an object with the Gene sequence
    if (Gene.coding == T) {
      Gene.sec <- getFrag(object.sec, merged.blast.result[Gen.num,paste(Genome,".Start", sep="")], merged.blast.result[Gen.num,paste(Genome,".End", sep="")])
    }    else {
      Gene.sec.rc <- getFrag(object.sec, merged.blast.result[Gen.num,paste(Genome,".End", sep="")], merged.blast.result[Gen.num,paste(Genome,".Start", sep="")])
      Gene.sec <- rev(comp(getSequence(Gene.sec.rc[[1]])))
    }
    ## create folder
    dir.create(path.multifasta, showWarnings = F)
    ## write fasta
    write.fasta(sequences = Gene.sec, names = Genome, file.out = paste(path.multifasta, "/", Gene, ".fasta", sep=""), open = "a", nbchar = 60, as.string = FALSE)
  }
  sapply(X = merged.blast.result[[1]], getGen)
}

#' Title
#'
#' @param path
#' @param threads
#'
#' @return
#' @export
#'
#' @examples
clustalo <- function(path = "HK.multif", threads = 8){
  path.to.files <- list.files(path, full.names = T)
  file.name <- list.files(path)
  dir.create(paste0(path,"/aln"), showWarnings = F)
  clines <- paste0("~/bin/clustalo -i ", path, "/", file.name," -o ", path,"/aln/", gsub("fasta","aln.fasta", file.name), " -t DNA --threads ", threads)
  sapply(clines, system)
}

mafft <- function(path = "multif", threads = 8){
  path.to.files <- list.files(path, pattern = ".fasta" full.names = T)
  file.name <- list.files(path)
  dir.create(paste0(path,"/aln"), showWarnings = F)
  clines <- paste0("mafft --auto ", path, "/", file.name," > ", path,"/aln/", gsub("fasta","aln.fasta", file.name))
  cl<- parallel::makeCluster(threads)
  parallel::parSapply(cl = cl, clines, system)
  parallel::stopCluster(cl)
}

#' Title
#'
#' @param path2Gblock
#' @param type
#' @param path
#' @param pattern
#' @param out.pattern
#'
#' @return
#' @export
#'
#' @examples
gBlock.R <- function(path2Gblock ="~/bin/Gblocks_0.91b/", type="d", path= ".", pattern= ".aln.fasta", out.pattern = ".gb"){
  input.aln <- list.files(path =  path, pattern = pattern, full.names = T)
  numfiles <- length(input.aln)
  paste(path2Gblock, "Gblocks ", input.aln," -t=", type, " -e=", out.pattern, " >> /dev/null", sep="")
}

#' Title
#'
#' @param path
#' @param pattern
#'
#' @return
#' @export
#'
#' @examples
remove_space <- function(path, pattern) {
  input.aln <- list.files(path = path, pattern = pattern, full.names = T)
  for(i in 1:length(input.aln)){
    fasta <- readLines(input.aln[i])
    fasta <- gsub(" ","", fasta)
    writeLines(fasta,con = input.aln[i])
  }
}

#' Title
#'
#' @param path
#' @param pattern
#' @param n
#'
#' @return
#' @export
#'
#' @examples
Concatenate_amas <- function(amas.file, path, pattern, n) {
  input.aln <- list.files(path = path, pattern = pattern, full.names = T)
  concatenate.core.cmd <- paste(input.aln, collapse = " ")
  concatenate.core.cmd <- paste0("python3.9 ", amas.file, " concat -f fasta -d dna -i ", concatenate.core.cmd, " -u fasta -p ", path, n,"_genes.partition.txt -t ",path, n, "_genes.fasta", sep = " ")
  system(concatenate.core.cmd)
  concatenate.file <- gsub("=", " = ", paste("DNA,", readLines(con= paste0(path, n, "_genes.partition.txt")), sep=""))
  writeLines(text=concatenate.file, con= paste0(path, n, "_genes.partition.txt"))
}
