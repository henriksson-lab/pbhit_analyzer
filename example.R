source("count.R")
source("grstats.R")

#Original code at https://github.com/henriksson-lab/malaria_crispr2024 ;
#this reposity provides cleaned up code suitable for reuse by other labs

################################################################################
##################### Estimate counts from raw data ############################
################################################################################

#To get some example data, you can download all our FASTQs from https://zenodo.org/records/14003839
allpooldir <- "./example_raw_data"

for(curpool in list.files(allpooldir)) {
  pooldir <- file.path(allpooldir,curpool)
  if(dir.exists(pooldir)){
    count_reads_pbhit(
      pooldir
    )
  }
}


################################################################################
############################ Fit growth rates ##################################
################################################################################

#This is ready-made count data, if you just want to test GR-analysis from counts (unzip the included file).
#Set to folder you 
allpooldir <- "./example_crispr_pools"

#This puts the output in a location for the visualizer to pick up
outdir <- "./viz"

#No comparisons in example data
list_cond_compare <- data.frame(
  cond1=c(),
  cond2=c()
)

#if we had data to compare, it would look like this
if(FALSE){
  list_cond_compare <- data.frame(
    cond1=c("P BL6",     "NP BL6",    "P BL6"),
    cond2=c("P RAG1KO",  "NP RAG1KO", "P IFNy")
  )
}

#List of pools, unless you wish to exclude any
listpools <- list.files(allpooldir)
  
#Perform all comparisons.
#Note that this WILL give a list of errors; this happens when certain data points
#do not allow a growth rate curve to be fitted. When this happens, we exclude
#this data as it is likely too noisy to be of use
make_grstats(
    listpools=listpools, 
    allpooldir=allpooldir, 
    list_cond_compare=list_cond_compare, 
    outdir=outdir
) 


################################################################################
########################## Visualize the results ###############################
################################################################################

# See start.R in viz/
#
# Note that localhost:8080 is the default place for the Shiny server



