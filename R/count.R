
library(ggplot2)
library(stringr)

################################################################################
# Generate count matrix for barseq
# 
count_reads_barseq <- function(curpool) {

  for(curpool in listpools){
    
    print(curpool)
    
    # re.compile("GGCGG"+"(\w{8,16})"+"CTGAC")
    
    seqbefore <- "TAGTCGCAGTAGGCGG"
    
    allpooldir <- "/corgi/otherdataset/ellenbushell/crispr_pools"
    pooldir <- file.path(allpooldir, curpool)
    fastqdir <- file.path(pooldir, "fastq")
    bcfile <- file.path(pooldir, "used_bc.csv")
    countfile <- file.path(pooldir,"counts.RDS")
    
    frank_bc <- read.csv("/corgi/otherdataset/ellenbushell/barcode_to_gene_210920_FRANK.csv")
    frank_bc$sgrna <- str_split_fixed(frank_bc$gene,"\\|",2)[,1]
    frank_bc$seq <- str_to_upper(frank_bc$barcode)
    
    #Subset by the BCs expected here
    usedbc <- read.csv(bcfile,sep="\t")
    usedbc <- frank_bc[frank_bc$sgrna %in% usedbc$gene,]
    bclength <- str_length(usedbc$seq[1]) #  TCTTTTCCCAG
    
    #R1 needs RC
    usedbc$seq <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(usedbc$seq)))
    
    
    list_bclist <- list()
    for(onef in list.files(fastqdir)){
      if(str_ends(onef,"R1_001.fastq.gz")){
        print(onef)
        onep <- pipe(paste("zcat",file.path(fastqdir,onef),"| grep TAGTCGCAGTAGGCGG"))
        li <- readLines(onep)
        close(onep)
        
        #revcomp BC for barseq R1??
        
        bclist <- str_split_fixed(li, seqbefore,2)[,2]
        bclist <- data.frame(bc=str_sub(bclist,1,bclength))
        bclist <- sqldf::sqldf("select count(bc) as cnt, bc from bclist group by bc order by cnt desc")
        bclist$file <- onef
        
        list_bclist[[onef]] <- bclist
      }
    }
    counts <- do.call(rbind,list_bclist)
    counts <- reshape2::acast(counts, bc~file, value.var = "cnt", fill = 0)
    
    #
    counts <- counts[order(rowSums(counts), decreasing = TRUE),]
    counts <- counts[rownames(counts) %in% usedbc$seq,]
    colnames(counts) <- str_sub(colnames(counts),1,11)  #dangerous. split by _S instead?
    
    rownames(usedbc) <- usedbc$seq
    rownames(counts) <- usedbc[rownames(counts),]$sgrna
    
    saveRDS(counts, countfile)
  }
  
}


################################################################################
####################### Generate count matrices -- CRISPR ######################
################################################################################

################################################################################
# Generate count matrix for PbHiT
# 
count_reads_pbhit <- function(curpool, allpooldir="/corgi/otherdataset/ellenbushell/crispr_pools") {
  seqbefore <- "CAATATTATT"
  
  pooldir <- file.path(allpooldir, curpool)
  fastqdir <- file.path(pooldir, "fastq")
  bcfile <- file.path(pooldir, "bc.csv")
  countfile <- file.path(pooldir,"counts.RDS")
  
  usedbc <- read.csv(bcfile,sep="\t")
  bclength <- str_length(usedbc$seq[1])
  
  list_bclist <- list()
  for(onef in list.files(fastqdir)){
    print(onef)
    if(str_ends(onef,"R1_001.fastq.gz")){
      onep <- pipe(paste("zcat",file.path(fastqdir,onef),"| grep CAATATTATT"))
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
  
  #getting out the right barcodes  
  counts <- counts[order(rowSums(counts), decreasing = TRUE),]
  counts <- counts[rownames(counts) %in% usedbc$seq,]
  
  #Keep name up to _S##
  colnames(counts) <- str_split_fixed(colnames(counts),"_S",2)[,1]

  rownames(usedbc) <- usedbc$seq
  rownames(counts) <- usedbc[rownames(counts),]$sgrna
  
  saveRDS(counts, countfile)  
}




