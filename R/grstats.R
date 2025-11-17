library(minpack.lm)
library(ggplot2)
library(stringr)
library(umap)


### how to put in??
list_cond_compare <- data.frame(
  cond1=c("P BL6",     "NP BL6",    "P BL6"),
  cond2=c("P RAG1KO",  "NP RAG1KO", "P IFNy")
)

outdir <- "/corgi/websites/malaria_crispr2024"

################################################################################
# Perform statistics on count tables
make_grstats <- function(listpools, allpooldir, list_cond_compare, outdir) {
  
  timecourses <- list()
  all_grstats_per_grna <- list()
  all_grstats <- list()
  list_samplemeta <- list()
  all_coverage_stat <- list()
  for(curpool in listpools){
    
    print(curpool)
    
    pooldir <- file.path(allpooldir, curpool)
    countfile <- file.path(pooldir,"counts.RDS")
    samplemetafile <- file.path(pooldir,"sampleinfo.txt")
    controlmetafile <- file.path(pooldir,"list_control.csv")
    
    #### Read sample metadata
    samplemeta <- read.csv(samplemetafile, sep = "\t")[,1:2]
    colnames(samplemeta) <- c("sampleid","samplename")
    samplemeta$day <- str_sub(str_split_fixed(samplemeta$samplename, "_",5)[,4],2)
    samplemeta$is_input <- str_count(samplemeta$samplename,"input")>0
    samplemeta$day <- as.integer(samplemeta$day)
    samplemeta$mouse_ref <- str_split_fixed(samplemeta$samplename, "_",5)[,5]
    samplemeta$genotype <- str_split_fixed(samplemeta$samplename, "_",5)[,3] ##"wt"
    samplemeta$primed <- str_split_fixed(samplemeta$samplename, "_",5)[,2]
    
    if(sum(samplemeta$is_input)>0){
      samplemeta$day[samplemeta$is_input] <- NA 
      samplemeta$mouse_ref[samplemeta$is_input] <- NA 
      samplemeta$genotype[samplemeta$is_input] <- NA 
      samplemeta$primed[samplemeta$is_input] <- NA 
    }  
    
    ### Read count table
    counts <- readRDS(countfile)
    
    #### Read info about the cloning
    cloningfile <- file.path(pooldir,"cloning.csv")
    if(file.exists(cloningfile)){
      print("Reading cloning.csv")
      allgeneconstructs <- read.csv(cloningfile,sep="\t")  
      allgeneconstructs$gene <- str_split_fixed(allgeneconstructs$grna,"gRNA",2)[,1]
      allgeneconstructs$genecat[allgeneconstructs$genecat==""] <- "Unstudied"
    } else {
      print("No cloning.csv -- constructing equivalent")
      
      allgeneconstructs <- data.frame(
        grna=rownames(counts),
        gene=str_split_fixed(rownames(counts),"00gRNA",2)[,1]  ### is this ok??? hack!
      )
      allgeneconstructs$genecat <- "Unstudied"
      allgeneconstructs$genewiz <- "NA"
      allgeneconstructs$ligationwell <- "NA"
      
      geneinfotable <- read.csv("/corgi/otherdataset/ellenbushell/gene_description.csv",sep="\t")
      allgeneconstructs$genecat[allgeneconstructs$gene %in% geneinfotable$gene[geneinfotable$genedesc=="Dispensable"]] <- "Dispensable"
      allgeneconstructs$genecat[allgeneconstructs$gene %in% geneinfotable$gene[geneinfotable$genedesc=="Slow"]] <- "Slow"
    }
    
    grna_dispensible <- allgeneconstructs$grna[allgeneconstructs$genecat=="Dispensable"]
    genes_dispensible <- allgeneconstructs$gene[allgeneconstructs$genecat=="Dispensable"]
    
    
    ### Extract FIRST input
    if(sum(samplemeta$is_input)>0){
      print("Using input sample as coverage")
      input_sampleid <- samplemeta$sampleid[samplemeta$is_input][1]
      coverage_stat <- data.frame(
        grna=rownames(counts),
        cnt=counts[,input_sampleid]
      )
    } else {
      print("Using average sample as coverage")
      coverage_stat <- data.frame(
        grna=rownames(counts),
        cnt=rowSums(counts)
      )
    }
    coverage_stat <- merge(allgeneconstructs,coverage_stat)
    all_coverage_stat[[curpool]] <- coverage_stat
    
    ### Make pseudocounts
    counts <- counts + 1
    
    #Gather total count
    rownames(samplemeta) <- samplemeta$sampleid
    samplemeta <- samplemeta[colnames(counts),]
    samplemeta$total_count <- colSums(counts)
    
    #Filter bad samples
    count_stats <- colSums(counts)
    bad_libs <- count_stats<0 #was 800 #was 10000, later 2000 (still lost a lot of samples in barseq, email)
    if(sum(bad_libs)>0){
      print("bad libs! here are counts before")
      print(count_stats)
      print("Removing:")
      print(colnames(counts)[bad_libs])
      counts <- counts[,!bad_libs]
      #print("bad libs! here are counts after")
      #count_stats <- colSums(counts)
      #print(count_stats)
    }
    
    ######## Figure out what control to use
    do_control_compensation <- TRUE
    if(file.exists(controlmetafile)){
      print("using list of controls-file")
      list_controls <- read.csv(controlmetafile)[,1]
    } else {
      print("using dispensable as controls")
      list_controls <- unique(rownames(counts)[rownames(counts) %in% grna_dispensible])   #was genes_dispensible, bug!
    }
    print(list_controls)
    
    
    #Normalize each library by depth
    for(i in 1:ncol(counts)){
      counts[,i] <- counts[,i]/sum(counts[,i]) 
    }
    
    
    this_timecourse <- list()
    
    ### Keep these counts for visualization
    #### Merge metadata with counts; remove input
    longcnt_sf <- reshape2::melt(as.matrix(counts))
    colnames(longcnt_sf) <- c("grna","sampleid","y") 
    longcnt_sf <- merge(longcnt_sf, samplemeta[!is.na(samplemeta$day),])
    longcnt_sf$gene <- str_split_fixed(longcnt_sf$grna,"gRNA",2)[,1]
    #longcnt_sf$count_type <- "Count/AllCount"
    longcnt_sf <- merge(longcnt_sf,unique(allgeneconstructs[,c("grna","genecat")]))
    this_timecourse[["Count/AllCount"]] <- longcnt_sf
    
    
    print(paste("Detected genes:", length(unique(longcnt_sf$gene))))
    print(paste("Detected sgRNAs:", length(unique(longcnt_sf$grna))))
    
    if(length(list_controls)==0){
      print("No control genes in list, so just normalizing by total")
    } else {
      #Normalize each library by sum of controls
      print("Have control genes, so normalizing by control sum")
      for(i in 1:ncol(counts)){
        counts[,i] <- counts[,i]/sum(counts[rownames(counts) %in% list_controls,i]) 
      }
      
      ### Keep these counts for visualization
      #### Merge metadata with counts; remove input
      longcnt_control <- reshape2::melt(as.matrix(counts))
      colnames(longcnt_control) <- c("grna","sampleid","y") 
      longcnt_control <- merge(longcnt_control, samplemeta[!is.na(samplemeta$day),])
      longcnt_control$gene <- str_split_fixed(longcnt_control$grna,"gRNA",2)[,1]
      longcnt_control$count_type <- "Count/ControlCount"
      longcnt_control <- merge(longcnt_control,unique(allgeneconstructs[,c("grna","genecat")]))
      this_timecourse[["Count/ControlCount"]] <- longcnt_control
    }
    
    timecourses[[curpool]] <- this_timecourse
    
    #Align samplemeta with counts. Compute UMAP
    rownames(samplemeta) <- samplemeta$sampleid
    samplemeta <- samplemeta[colnames(counts),]
    umap.settings <- umap.defaults
    umap.settings$n_neighbors <- min(umap.settings$n_neighbors, ncol(counts))
    cnt.umap <- umap(t(counts), config=umap.settings)
    samplemeta$umap1 <- cnt.umap$layout[,1]
    samplemeta$umap2 <- cnt.umap$layout[,2]
    if(FALSE){
      ggplot(samplemeta, aes(umap1,umap2, label=paste(sampleid)))+ geom_point(color="gray") + geom_text()
    }
    #table(sort(samplemeta$samplename))

    ####### Merge metadata with counts.
    ####### remove input samples from TC
    longcnt <- reshape2::melt(as.matrix(counts))
    colnames(longcnt) <- c("grna","sampleid","cnt") 
    longcnt <- merge(longcnt, samplemeta[!is.na(samplemeta$day),])
    longcnt$gene <- str_split_fixed(longcnt$grna,"gRNA",2)[,1]
    
    ######## Compute GRs (RGRs because of previous normalization already)
    
    #For each phenotype
    fitted_gr <- list()
    for(thepheno in unique(longcnt$primed)){
      #print(thepheno)
      sub_longcnt2 <- longcnt[longcnt$primed==thepheno,,drop=FALSE]
      #For each genotype
      for(thegeno in unique(longcnt$genotype)){
        #print(thegeno)
        sub_longcnt3 <- sub_longcnt2[sub_longcnt2$genotype==thegeno,,drop=FALSE]
        #For each mouse
        for(themouse in unique(longcnt$mouse_ref)){
          #print(themouse)
          sub_longcnt4 <- sub_longcnt3[sub_longcnt3$mouse_ref==themouse,,drop=FALSE]
          #For each grna
          for(thegrna in unique(longcnt$grna)){
            #print(thegrna)
            sub_longcnt5 <- sub_longcnt4[sub_longcnt4$grna==thegrna,,drop=FALSE]
            curstate <- paste(curpool, themouse, thegrna, thegeno, thepheno)

            if(nrow(sub_longcnt5)>0 & nrow(sub_longcnt2)>0){
              
              
              result <- try({
                
                ############ Figure out weights              
                #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5527095/  8x over 48h
                #https://link.springer.com/article/10.1007/s00436-009-1435-8   bergei , doubling rate is every 6-8h. so 2*2*2*2=16x per day
                malaria_growth_rate <- sqrt(8) ##per day
                malaria_rel_amount <- data.frame(day=1:5)
                malaria_rel_amount$approx_amount <- malaria_growth_rate**malaria_rel_amount$day
                malaria_rel_amount$weight <- 1 #sqrt(malaria_rel_amount$day) ###sqrt(malaria_rel_amount$approx_amount)
                #plot(malaria_rel_amount$weight)
                #weight <- 1/(1/malaria_rel_amount$approx_amount)  #variance is ~1/approx_amount ;
                #https://en.wikipedia.org/wiki/Inverse-variance_weighting
                
                sub_longcnt5$weight <- sub_longcnt5$day
                
                # fit the model 
                x <- sub_longcnt5$day
                y <- sub_longcnt5$cnt
                if(FALSE){
                  start_values
                  plot(x,y)
                }
                
                
                if(FALSE){
                  ############### Exponential growth; R standard nls
                  start_values <- c(a=mean(y), b=0)
                  fit <- nls(y ~ a * exp(b * x),
                             start = start_values,
                             algorithm = "port",
                             control = nls.control(maxiter = 1000))  
                  slope_mean <- as.data.frame(coef(summary(fit)))[2,"Estimate"]
                  slope_sd <- as.data.frame(coef(summary(fit)))[2,"Std. Error"]
                } 
                
                if(TRUE) {
                  ########## logistic growth, but with pre-fixed maximum capacity to 1
                  #https://www.usu.edu/math/powell/ysa-html/node8.html
                  
                  scalefactor <- max(y)
                  y <- y/scalefactor
                  
                  ############### Levenberg-Marquardt nls
                  #The more detailed help here. Levenberg-Marquardt from MINPACK
                  #?minpack.lm::nls.lm
                  start_values <- c(a=mean(y), b=0)
                  fit <- nlsLM(y ~ a * exp(b * x) / (1 + a * exp(b * x) ),
                               start = start_values,
                               algorithm = "port",
                               control = nls.control(maxiter = 1000))  #was 1000
                  slope_mean <- as.data.frame(coef(summary(fit)))[2,"Estimate"]
                  slope_sd <- as.data.frame(coef(summary(fit)))[2,"Std. Error"]
                }
                
                
                if(FALSE){
                  ########## logistic growth, nlsLM; fit maximum capacity  ## previous favourite
                  #https://www.usu.edu/math/powell/ysa-html/node8.html
                  
                  start_values <- c(a=mean(y), b=0, d=max(y))
                  fit <- nlsLM(y ~ a * exp(b * x) / (1 + a/d * exp(b * x) ),
                               start = start_values,
                               algorithm = "port",
                               control = nls.control(maxiter = 1000))  #was 1000
                  slope_mean <- as.data.frame(coef(summary(fit)))[2,"Estimate"]
                  slope_sd <- as.data.frame(coef(summary(fit)))[2,"Std. Error"]
                }
                
                
                if(FALSE) {
                  #Exponential growth, nlsLM
                  #Comparison of methods https://rdrr.io/rforge/nlsr/f/inst/doc/nlsr-nls-nlsLM.pdf
                  
                  start_values <- c(a=mean(y), b=0)
                  fit <- nlsLM(y ~ a * exp(b * x),
                               start = start_values,
                               algorithm = "port",
                               control = nls.control(maxiter = 1000))  #was 1000
                  slope_mean <- as.data.frame(coef(summary(fit)))[2,"Estimate"]
                  slope_sd <- as.data.frame(coef(summary(fit)))[2,"Std. Error"]
                }
                
                if(FALSE){
                  #Exponential growth, nlsr::nlxb
                  #Comparison of methods https://rdrr.io/rforge/nlsr/f/inst/doc/nlsr-nls-nlsLM.pdf
                  start_values <- c(a=mean(y), b=0)
                  fit <- nlsr::nlxb(y ~ a * exp(b * x),
                                    start = start_values
                  ) 
                  slope_mean <-  as.double(coef(fit)[2])
                  thesummary <- summary(fit)
                  slope_sd <- thesummary$Sd[2]  #disagrees a fair bit with above
                }
                
                
                fitted_gr[[paste(themouse, thegrna, thegeno, thepheno)]] <- data.frame(
                  mouse=themouse,
                  grna=thegrna,
                  phenotype=thepheno,
                  genotype=thegeno,
                  gr_sd=slope_sd,
                  gr_mean=slope_mean
                )
                #print(paste("Worked",curstate))
              }, silent = FALSE) ## Ignore those we cannot fit
              
            }
            
          }
        }
        
      }
    }
    fitted_gr <- do.call(rbind, fitted_gr)
    if(is.null(fitted_gr)){
      stop("fitted_gr null, so empty; no models converged")
    }
    
    print(paste("Fitted # grna:",length(unique(fitted_gr$grna))))
    print(paste("Fitted # gene:",length(unique(str_split_fixed(fitted_gr$grna,"gRNA",2)[,1]))))
    

    ##############################################################################
    ################### Calculate FCs to control and between conditions ##########
    ##############################################################################
    
    computeLogpFromVar <- function(fc, sd){
      p <- pnorm(fc,0,sd=sd)
      ind_flip <- which(p>0.5) #handles NA
      p[ind_flip] <- 1 - p[ind_flip]
      logp <- -log10(p)
      logp
    }
    
    ############## Function to get variance for one guide
    computeVar <- function(fitted_gr){
      
      ##### PER GRNA CONSTRUCT: Calculate FC vs control
      grna_var <- sqldf::sqldf("select sum(gr_sd*gr_sd) as totalvar, count(gr_sd) as cnt, avg(gr_mean) as fc, grna from fitted_gr group by grna")
      grna_var$sd <- sqrt(grna_var$totalvar/grna_var$cnt)
      grna_var$logp <- computeLogpFromVar(grna_var$fc,grna_var$sd)
      
      grna_var <- merge(grna_var,allgeneconstructs) #maybe an outer join later? TODO
      grna_var <- grna_var[order(grna_var$genecat, grna_var$fc),]
      all_grstats_per_grna[[curpool]] <- grna_var
      
      ##### PER gene: Calculate FC vs control
      gene_var <- sqldf::sqldf("select sum(sd*sd) as totalvar, count(sd) as cnt, avg(fc) as fc, gene, genecat from grna_var group by gene")
      gene_var$sd <- sqrt(gene_var$totalvar/gene_var$cnt)
      gene_var$logp <- computeLogpFromVar(gene_var$fc,gene_var$sd)
      
      list(
        grna_var=grna_var,
        gene_var=gene_var
      )
    }
    
    ######### All comparisons vs control  
    list_all_grstats <- list()
    list_all_grstats_per_grna <- list()
    #For each phenotype
    for(thepheno in unique(longcnt$primed)){
      #For each genotype
      for(thegeno in unique(longcnt$genotype)){
        print(paste(thepheno, thegeno, "vs control"))
        vscontrol <- computeVar(fitted_gr[fitted_gr$genotype==thegeno & fitted_gr$phenotype==thepheno,])
        list_all_grstats_per_grna[[paste(thepheno, thegeno)]] <- vscontrol$grna_var
        list_all_grstats[[paste(thepheno, thegeno)]] <- vscontrol$gene_var
      }
    }
    
    
    ######### All pairwise comparisons
    tocompare <- names(list_all_grstats)
    for(curcond1 in tocompare){
      for(curcond2 in tocompare){
        if(curcond1!=curcond2){
          var1 <- list_all_grstats[[curcond1]]
          var2 <- list_all_grstats[[curcond2]]
          
          var1 <- data.frame(sd1=var1$sd, fc1=var1$fc, gene=var1$gene, genecat=var1$genecat)
          var2 <- data.frame(sd2=var2$sd, fc2=var2$fc, gene=var2$gene)
          
          varcomp <- merge(var1,var2)
          varcomp$fc <- varcomp$fc1 - varcomp$fc2
          varcomp$sd <- sqrt(varcomp$sd1**2 + varcomp$sd2**2)
          varcomp$logp <- computeLogpFromVar(varcomp$fc,varcomp$sd)
          
          print(paste(curcond1, "-", curcond2))
          list_all_grstats[[paste(curcond1, "-", curcond2)]] <- varcomp
        }
      }
    }
    
    
    ##############################################################################
    ############ Add comparisons between conditions ##############################
    ##############################################################################
    
    list_scatter <- list()
    for(i in 1:nrow(list_cond_compare)){
      if(list_cond_compare$cond1[i] %in% names(list_all_grstats) & list_cond_compare$cond2[i] %in% names(list_all_grstats)){
        compname <- sprintf("%s -- %s",list_cond_compare$cond1[i],list_cond_compare$cond2[i])
        print(compname)
        stat1 <- list_all_grstats[[list_cond_compare$cond1[i]]]
        stat2 <- list_all_grstats[[list_cond_compare$cond2[i]]]
        
        toplot <- merge(
          data.frame(
            genedesc=stat1$genecat,
            gene=stat1$gene,
            sd1=stat1$sd,
            p1=stat1$logp,
            fc1=stat1$fc),
          
          data.frame(
            gene=stat2$gene,
            sd2=stat2$sd,
            p2=stat2$logp,
            fc2=stat2$fc)
        )
        
        #Estimate p-value of difference, assuming a normal distribution
        toplot$diff_fc <- toplot$fc1 - toplot$fc2
        toplot$diff_sd <- sqrt(toplot$sd1**2 + toplot$sd2**2)
        toplot$diff_log_p <- computeLogpFromVar(toplot$diff_fc, toplot$diff_sd)
        
        list_scatter[[compname]] <- toplot
      }
      
    }
    
    all_grstats[[curpool]] <- list(
      volcano=list_all_grstats,
      stats_per_grna=list_all_grstats_per_grna,
      scatterplot=list_scatter
    )
    
    list_samplemeta[[curpool]] <- samplemeta
  }
  
  
  
  
  saveRDS(all_grstats, file=file.path(outdir,"grstats.rds"))
  saveRDS(timecourses, file=file.path(outdir,"timecourses.rds"))
  saveRDS(list_samplemeta, file=file.path(outdir,"samplemeta.rds"))
  saveRDS(all_coverage_stat, file=file.path(outdir,"coverage_stat.rds"))
}



