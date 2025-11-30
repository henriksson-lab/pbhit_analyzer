
library(ggplot2)
library(stringr)

################################################################################
####################### Generate count matrices from raw FASTQ #################
################################################################################


################################################################################
#' Generate count data for barseq FASTQ
#' 
#' @param pooldir Directory representing one pool. See README for more information 
#' @return Nothing
#' 
count_reads_barseq <- function(
    pooldir
) {

    seqbefore <- "TAGTCGCAGTAGGCGG"
    
    fastqdir <- file.path(pooldir, "fastq")
    bcfile <- file.path(pooldir, "bc.csv")
    countfile <- file.path(pooldir,"counts.RDS")
    
    frank_bc <- read.csv("barcode_to_gene_210920_FRANK.csv")
    frank_bc$sgrna <- str_split_fixed(frank_bc$gene,"\\|",2)[,1]
    frank_bc$seq <- str_to_upper(frank_bc$barcode)
    
    #Subset by the BCs expected here
    usedbc <- read.csv(bcfile,sep="\t")
    usedbc <- frank_bc[frank_bc$sgrna %in% usedbc$gene,]
    bclength <- str_length(usedbc$seq[1]) #  TCTTTTCCCAG
    
    #R1 needs reverse complement
    usedbc$seq <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(usedbc$seq)))
    
    #For each of the FASTQ files
    list_bclist <- list()
    for(onef in list.files(fastqdir)){
#      if(str_ends(onef,"R1_001.fastq.gz")){
      if(str_detect(onef,"_R1_") && str_ends(onef, "fastq.gz")){
        print(onef)
        
        #Extract reads that seem relevant. Using zcat speeds up parsing a lot vs doing it all in R
        onep <- pipe(paste("zcat",file.path(fastqdir,onef),"| grep TAGTCGCAGTAGGCGG"))
        li <- readLines(onep)
        close(onep)
        
        bclist <- str_split_fixed(li, seqbefore,2)[,2]
        bclist <- data.frame(bc=str_sub(bclist,1,bclength))
        bclist <- sqldf::sqldf("select count(bc) as cnt, bc from bclist group by bc order by cnt desc")
        bclist$file <- onef
        
        list_bclist[[onef]] <- bclist
      }
    }
    counts <- do.call(rbind,list_bclist)
    counts <- reshape2::acast(counts, bc~file, value.var = "cnt", fill = 0)
    
    #Produce count table
    counts <- counts[order(rowSums(counts), decreasing = TRUE),]
    counts <- counts[rownames(counts) %in% usedbc$seq,]
    colnames(counts) <- str_sub(colnames(counts),1,11)  #The current approach uses the first part of the name. In the future, better to split at _S instead
    
    rownames(usedbc) <- usedbc$seq
    rownames(counts) <- usedbc[rownames(counts),]$sgrna
    
    saveRDS(counts, countfile)
#  }
  
}


################################################################################
#' Generate count matrix for PbHiT
#' 
#' @param pooldir Directory representing one pool. See README for more information 
#' @return Nothing
#' 
count_reads_pbhit <- function(
    pooldir
    #curpool, 
    #allpooldir
) {

  #This sequence is right before the gRNA sequence
  seqbefore <- "CAATATTATT"
  
  #pooldir <- file.path(allpooldir, curpool)
  fastqdir <- file.path(pooldir, "fastq")
  bcfile <- file.path(pooldir, "bc.csv")     
  countfile <- file.path(pooldir,"counts.RDS")
  
  usedbc <- read.csv(bcfile,sep="\t")
  bclength <- str_length(usedbc$seq[1])
  
  #For each of the FASTQ files
  list_bclist <- list()
  for(onef in list.files(fastqdir)){
    print(onef)
    if(str_detect(onef,"_R1_") && str_ends(onef, "fastq.gz")){
    #if(str_ends(onef,"R1_001.fastq.gz")){
      
      #Extract reads that seem relevant. Using zcat speeds up parsing a lot vs doing it all in R
      onep <- pipe(paste("zcat",file.path(fastqdir,onef),"| grep CAATATTATT"))
      li <- readLines(onep)
      close(onep)
      
      #Count barcodes
      bclist <- str_split_fixed(li, seqbefore,2)[,2]
      bclist <- data.frame(bc=str_sub(bclist,1,bclength))
      bclist <- sqldf::sqldf("select count(bc) as cnt, bc from bclist group by bc order by cnt desc")
      bclist$file <- onef
      
      list_bclist[[onef]] <- bclist
    }
  }
  
  #Get counts across all files into one big matrix. Using do.call
  #is much faster than cbind-ing each file
  counts <- do.call(rbind,list_bclist)
  counts <- reshape2::acast(counts, bc~file, value.var = "cnt", fill = 0)
  
  #Subset to only the correct barcodes  
  counts <- counts[order(rowSums(counts), decreasing = TRUE),]
  counts <- counts[rownames(counts) %in% usedbc$seq,]
  
  #Keep name up to _S##_ (i.e. get rid of sample number and beyond)
  colnames(counts) <- str_split_fixed(colnames(counts),"_S",2)[,1]

  rownames(usedbc) <- usedbc$seq
  rownames(counts) <- usedbc[rownames(counts),]$sgrna
  
  saveRDS(counts, countfile)  
}


