
library(ggplot2)
library(stringr)

#previous counter in /corgi/websites/malariascreenviewer/countcrispr.R 


# Counts/Total
# Counts/Reference
# Normalized growth rate


# Weight by parasitemia
# Shrinkage by average depth

################################################################################
####################### Generate count matrices -- barseq ######################
################################################################################


listpools <- c(
  #  "barseq_minipool2",  #just a cnt file!
  "barseq_priming_Candidatepool1",
  "barseq_priming_Candidatepool2",
  "barseq_priming_barseqpool1",
  "barseq_priming_barseqpool3",
  "barseq_priming_barseqpool4",
  "barseq_slowhires_2023dec",
  "barseq_slowpool_1",
  "barseq_slowpool_2"
)
#listpools <- c("barseq_slowhires_2023dec")
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


################################################################################
####################### Hack for CRISPR minipool2 ##############################
################################################################################

if(FALSE){
  counts <- read.csv("/corgi/otherdataset/ellenbushell/crispr_pools/barseq_minipool2/counts.csv", row.names = "X")
  counts <- counts[rownames(counts) != "_other" & rownames(counts) != "_err",]
  colnames(counts) <- str_split_fixed(colnames(counts),"\\.",2)[,1]
  saveRDS(counts, "/corgi/otherdataset/ellenbushell/crispr_pools/barseq_minipool2/counts.RDS")
}



################################################################################
####################### Generate count matrices -- CRISPR ######################
################################################################################

listpools <- c("cr_2024march_half1","cr_2024march_p1","cr_2024march_p12","cr_2024march_p2")

# "tags" = EB_barseq_slowpool_CRISPR


#listpools <- c("tags")
listpools <- c("2023march_screen_noD4")
listpools <- c("2023march_screen")
#listpools <- c("aug_p192","aug_p24","aug_p96")
listpools <- c("2023aug_p192","2023aug_p24","2023aug_p96","2023jan_tags", "2023march_screen", "2023march_screen_noD4")
listpools <- c("2023dec_ligpool")
for(curpool in listpools){
  
  print(curpool)
  
  seqbefore <- "CAATATTATT"
  
  allpooldir <- "/corgi/otherdataset/ellenbushell/crispr_pools"
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
  #colnames(counts) <- str_sub(colnames(counts),1,11) ##### only applies if
  
  rownames(usedbc) <- usedbc$seq
  rownames(counts) <- usedbc[rownames(counts),]$sgrna
  
  saveRDS(counts, countfile)
}



listpools <- c("cr_2024march_half1","cr_2024march_p1","cr_2024march_p12","cr_2024march_p2")
for(curpool in listpools){
  print(curpool)
  allpooldir <- "/corgi/otherdataset/ellenbushell/crispr_pools"
  pooldir <- file.path(allpooldir, curpool)
  countfile <- file.path(pooldir,"counts.RDS")
  counts <- readRDS(countfile)
  print(rowMeans(counts))
}


