library(shiny)

library(stringr)
library(ggplot2)
library(patchwork)
library(sqldf)
library(reshape2)
library(cowplot)


print("======= reading data to be cached in memory ================ ")

if(FALSE){
  getID <- function(s){
    s <- str_split(s,";")[[1]]
    w <- which(str_starts(s,"ID"))
    if(length(w)>0){
      str_replace_all(s[w[1]],"ID=","")
    } else {
      ""
    }
  }
  getDESC <- function(s){
    s <- str_split(s,";")[[1]]
    w <- which(str_starts(s,"description"))
    if(length(w)>0){
      str_replace_all(s[w[1]],"description=","")
    } else {
      ""
    }
  }
  getNAME <- function(s){
    s <- str_split(s,";")[[1]]
    w <- which(str_starts(s,"Name"))
    if(length(w)>0){
      str_replace_all(s[w[1]],"Name=","")
    } else {
      ""
    }
  }
  gff <- read.table("PlasmoDB-67_PbergheiANKA.gff.gz",sep="\t", quote="")
  gff <- gff[str_detect(gff$V9,"description"),]$V9
  gff <- gff[which(str_detect(gff,"ID"))]
  geneinfo <- data.frame(
    gene=sapply(gff, getID),
    geneDesc=sapply(gff, getDESC),
    geneName=sapply(gff, getNAME)
  )
  geneinfo$geneDesc <- str_replace_all(geneinfo$geneDesc,"%2C","")
  geneinfo <- unique(geneinfo)
  geneinfo  
  saveRDS(geneinfo, "geneinfo.rds")
  
} else {
  #Above is slow enough that we rather not do it every time the docker container
  #is restarted!
  geneinfo <- readRDS("geneinfo.rds")
}







print("======= reading data to be cached in memory ================ ")

all_samplemeta <- readRDS("samplemeta.rds")
all_grstats <- readRDS("grstats.rds")
all_timecourses <- readRDS("timecourses.rds")
all_coverage_stat <- readRDS("coverage_stat.rds")
all_composite <- readRDS("composite.rds")


print("======= renaming pools ================ ")

pools_renamed <- list(
#  cr_2023march_pools="cr_2023march_pools",
  cr_2023march_screen="22x",
  cr_2024march_half1="48x",
  cr_2024march_p1="96x",
#  cr_2024march_p2="p2",
  cr_2024march_p12="192x"
)

renamepool <- function(list_in){
  list_out <- list()
  for(n in names(list_in)){
    print(n)
    if(n %in% names(pools_renamed)){
      list_out[[pools_renamed[[n]]]] <- list_in[[n]]
    } else {
      print(paste("missing",n))
    }
  }
  list_out
}


all_samplemeta <- renamepool(all_samplemeta)
all_grstats <- renamepool(all_grstats)
all_timecourses <- renamepool(all_timecourses)
all_coverage_stat <- renamepool(all_coverage_stat)



print("========== global done ================")
