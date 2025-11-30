source("count.R")
source("grstats.R")

#Original code at https://github.com/henriksson-lab/malaria_crispr2024 ;
#this reposity provides cleaned up code suitable for reuse by other labs

################################################################################
##################### Estimate counts from raw data ############################
################################################################################





################################################################################
############################ Fit growth rates ##################################
################################################################################



################################################################################
########################## Visualize the results ###############################
################################################################################

# see start.R in 


#how to use stuff

#vignette dir?




if(FALSE) {
  
  
  
  if(FALSE){
    allpooldir <- "/corgi/otherdataset/ellenbushell/crispr_pools"
    for(curpool in listpools){
      #...
    }
  }
  
  

  count_reads_pbhit(
    curpool, 
    allpooldir="/corgi/otherdataset/ellenbushell/crispr_pools"
  )
    
    
    
  ### how to put in??
  list_cond_compare <- data.frame(
    cond1=c("P BL6",     "NP BL6",    "P BL6"),
    cond2=c("P RAG1KO",  "NP RAG1KO", "P IFNy")
  )
  
  outdir <- "/corgi/websites/malaria_crispr2024"
  
  
}