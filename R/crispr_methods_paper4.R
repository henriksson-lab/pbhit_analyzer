if(FALSE){
  install.packages("plotROC")
}

library(ggtext)
library(ggplot2)
library(grid)
library(egg)
library(plotROC)

################################################################################
############## Read data #######################################################
################################################################################

### Renaming, for display
pools_renamed <- list(
  cr_2023march_pools="cr_2023march_pools",
  cr_2023march_screen="22x",
  cr_2024march_half1="48x",
  cr_2024march_p1="96x",
  cr_2024march_p2="p2",
  cr_2024march_p12="192x"
)

all_samplemeta <- readRDS("/corgi/websites/malaria_crispr2024/samplemeta.rds")
all_grstats <- readRDS("/corgi/websites/malaria_crispr2024/grstats.rds")
all_timecourses <- readRDS("/corgi/websites/malaria_crispr2024/timecourses.rds")
all_coverage_stat <- readRDS("/corgi/websites/malaria_crispr2024/coverage_stat.rds")



################################################################################
############## Fig ????  Barchart, relative abundance, ligation tests ##########
################################################################################

counts <- readRDS("/corgi/otherdataset/ellenbushell/crispr_pools/cr_2023march_pools/counts.RDS")
meta <- read.csv("/corgi/otherdataset/ellenbushell/crispr_pools/cr_2023march_pools/sampleinfo.txt", sep="\t")

expected_gene <- read.csv("/corgi/otherdataset/ellenbushell/crispr_pools/cr_2023march_pools/expected.csv", sep="\t")
expected_gene[is.na(expected_gene)] <- FALSE
expected_gene <- rbind(
  data.frame(gene=expected_gene$gene[expected_gene$in8 =="TRUE"], ligationwell="pool8"),
  data.frame(gene=expected_gene$gene[expected_gene$in12=="TRUE"], ligationwell="pool12"),
  data.frame(gene=expected_gene$gene[expected_gene$in22=="TRUE"], ligationwell="pool22")
)

colnames(counts) <- meta$User.ID
coverage_stat <- melt(counts)
colnames(coverage_stat) <- c("grna","ligationwell","cnt")
coverage_stat <- coverage_stat[coverage_stat$cnt>0,]
coverage_stat$mouse <- str_split_fixed(coverage_stat$ligationwell,"_",5)[,5]
coverage_stat$ligationwell <- str_split_fixed(coverage_stat$ligationwell,"_",5)[,1]
coverage_stat$gene <- str_split_fixed(coverage_stat$grna,"gRNA",2)[,1]
coverage_stat <- merge(expected_gene, coverage_stat)
ligpools <- sqldf::sqldf("select mouse, grna, ligationwell, sum(cnt) as cnt from coverage_stat group by grna, ligationwell, mouse")
ligpools_sum <- sqldf::sqldf("select mouse, ligationwell, sum(cnt) as cnt_tot from coverage_stat group by ligationwell, mouse")
ligpools <- merge(ligpools_sum, ligpools)
ligpools$frac <- ligpools$cnt/ligpools$cnt_tot
ligpools$ligationwell <- factor(ligpools$ligationwell, levels = c("pool8","pool12","pool22"))

ligpools$well_rep <- paste(ligpools$ligationwell,str_replace_all(ligpools$mouse,"m","#"))
ligpools$well_rep <- factor(ligpools$well_rep, levels = rev(c(
  "pool8 #1","pool8 #2","pool8 #3",
  "pool12 #1","pool12 #2","pool12 #3",
  "pool22 #1","pool22 #2","pool22 #3"
  )))

if(FALSE){
  current_pool <- "cr_2023march_pools"
  ggplot(ligpools, aes(well_rep, frac*100, fill=grna)) + geom_bar(stat="identity", position = "stack") + 
    xlab("") + ylab("Fraction (%)")+
    theme_bw() + 
    theme(legend.position = "none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  ggsave(sprintf("/corgi/websites/malaria_crispr2024/plots_crispr/fraction_per_ligpool %s A.pdf",current_pool), width = 7, height = 3)  
}

current_pool <- "cr_2023march_pools"
ggplot(ligpools, aes(well_rep, frac*100, fill=grna)) + geom_bar(stat="identity", position = "stack") + 
  coord_flip()+
  xlab("") + ylab("Fraction (%)")+
  theme_bw() + 
  theme(legend.position = "none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave(sprintf("/corgi/websites/malaria_crispr2024/plots_crispr/fraction_per_ligpool %s B.pdf",current_pool), width = 4, height = 2)  



################################################################################
############## Fig ????  Barchart, relative abundance, pools ###################
################################################################################

for(current_pool in c("cr_2024march_half1","cr_2024march_p1","cr_2024march_p12")){ #"cr_2024march_p2"
  #current_pool <- "cr_2024march_half1"
  samplemeta <- all_samplemeta[[current_pool]]
  coverage_stat <- all_coverage_stat[[current_pool]]

  
  #Per ligation pool
  if(FALSE){
    ligpools <- sqldf::sqldf("select ligationwell, sum(cnt) as cnt from coverage_stat group by ligationwell")
    ligpools <- ligpools[ligpools$ligationwell!="spikein",]
    ligpools$frac <- ligpools$cnt/sum(ligpools$cnt)
    #ligpools$ligationwell <- str_sub(ligpools$ligationwell,2)
    ggplot(ligpools, aes(ligationwell, frac*100)) + geom_bar(stat="identity") + 
      geom_hline(yintercept=100/nrow(ligpools), color="blue")+
      xlab("") + ylab("Fraction (%)")+
      theme_bw() + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  }
  
  
  #Per grna
  ligpools <- sqldf::sqldf("select gene, ligationwell, sum(cnt) as cnt from coverage_stat group by gene, ligationwell")
  ligpools <- ligpools[ligpools$ligationwell!="spikein",]
  ligpools$frac <- ligpools$cnt/sum(ligpools$cnt)
  numpools <- length(unique(coverage_stat$ligationwell))
  ggplot(ligpools, aes(ligationwell, frac*100, fill=gene)) + geom_bar(stat="identity", position = "stack") + 
    geom_hline(yintercept=100/numpools, color="blue")+
    xlab("") + ylab("Fraction (%)")+
    theme_bw() + 
    theme(legend.position = "none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    

  ggsave(sprintf("/corgi/websites/malaria_crispr2024/plots_crispr/fraction_per_ligpool %s.pdf",current_pool), width = numpools*0.3, height = 3)  
}

#current_pool <- "cr_2023march_screen"

################################################################################
####################### Fig xXX. Scatter plot of sgRNA 1 vs 2 ##################
################################################################################

listpools <- c(
  "cr_2024march_half1",
  "cr_2024march_p1",
#  "cr_2024march_p2"
  "cr_2024march_p12"
)
listplot <- list()
for(curpool in listpools){
  print(curpool)
  grstats <- all_grstats[[curpool]]$stats_per_grna$`NP BL6`
  
  g1 <- grstats[str_ends(grstats$grna,"1"),]
  g2 <- grstats[str_ends(grstats$grna,"2"),]
  toplot <- merge(
    data.frame(gene=g1$gene, fc1=g1$fc, genecat=g1$genecat),
    data.frame(gene=g2$gene, fc2=g2$fc)
  )
  toplot$genecat <- factor(toplot$genecat, levels=c("Dispensable","Essential","Slow growers","Unstudied"))
  plot_range <- range(c(toplot$fc1, toplot$fc2))
  onep <- ggplot(toplot) + 
    xlab(paste("RGR sgRNA #1", pools_renamed[[curpool]]))+
    ylab("RGR sgRNA #2")+
    xlim(plot_range)+
    ylim(plot_range)+
    theme_bw()+
    theme(legend.position = "none")+
    geom_smooth(method = "lm", aes(fc1,fc2),color="black")+
    geom_point(aes(fc1,fc2,color=genecat),size=3)  +
    scale_color_manual(values = c("chartreuse4", "red", "dodgerblue", "#999999")) #"Dispensible","Essential","Slow growers","Unstudied"
  listplot[[curpool]] <- onep
}  
totp <- egg::ggarrange(plots=listplot)
totp
ggsave(plot=totp, "/corgi/websites/malaria_crispr2024/plots_crispr/scatterplot_sgrna_fc.pdf", width = 4, height = 10)



################################################################################
############### Fig xxx Comparison of screen RGRs, scatter plots ###############
################################################################################
  

listpools <- c(
#  "cr_2024march_half1",
  "cr_2024march_p1",
#  "cr_2024march_p2",
  "cr_2024march_p12"
)
listplot <- list()
for(curpool in listpools){
  print(curpool)

  grstats1 <- all_grstats[["cr_2024march_half1"]]$volcano$`NP BL6`
  grstats2 <- all_grstats[[curpool]]$volcano$`NP BL6`
  
  toplot <- merge(
    data.frame(
      gene=grstats1$gene,
      fc1=grstats1$fc,
      genecat=grstats1$genecat
    ),  
    data.frame(
      gene=grstats2$gene,
      fc2=grstats2$fc
    )
  )
  toplot$genecat <- factor(toplot$genecat, levels=c("Dispensable","Essential","Slow growers","Unstudied"))
  plot_range <- range(c(toplot$fc1, toplot$fc2))
  onep <- ggplot(toplot, aes(fc1,fc2,color=genecat)) + 
    geom_point() + 
    xlab(paste("RGR",pools_renamed[["cr_2024march_half1"]])) +
    ylab(paste("RGR",pools_renamed[[curpool]]))+
    #xlim(plot_range)+
    #ylim(plot_range)+
    theme_bw()+
    theme(legend.position = "none")+
    geom_smooth(method = "lm", aes(fc1,fc2),color="black")+  #95% confidence interval
    geom_point(aes(fc1,fc2,color=genecat),size=3)  +
    scale_color_manual(values = c("chartreuse4", "red", "dodgerblue", "#999999")) #"Dispensible","Essential","Slow growers","Unstudied"
  listplot[[curpool]] <- onep
  
  
  cor(toplot$fc1,toplot$fc2) ####
}  
totp <- egg::ggarrange(plots=listplot)
totp

ggsave(plot=totp, "/corgi/websites/malaria_crispr2024/plots_crispr/scatterplot_pool_reproducibility.pdf", width = 3, height = 5)




################################################################################
################## Fig xxx Composite analysis of RGR ###########################
################################################################################


listpools <- c(
  "cr_2024march_half1",
  "cr_2024march_p1"
# , "cr_2024march_p2"
#  "cr_2024march_p12"
)

#### Merge all pools
poolstats <- NULL
for(curpool in listpools){
  print(curpool)
  grstats <- all_grstats[[curpool]]$volcano$`NP BL6`
  grstats$pool <- curpool
  grstats$fc <- grstats$fc - mean(grstats$fc[grstats$genecat=="Dispensable"])
  #grstats$fc <- grstats$fc / -mean(grstats$fc[grstats$genecat=="Essential"])  #this will not work for p2
  poolstats <- rbind(poolstats, grstats)
}

avgpool <- sqldf::sqldf("select count(*) as cnt, sum(sd*sd) as totvar, avg(fc) as fc, gene, genecat from poolstats group by gene")
avgpool$sd <- sqrt(avgpool$totvar/avgpool$cnt)

#### Distribution of Essential
fc_ess <- avgpool$fc[avgpool$genecat=="Essential"]
mean_ess <- mean(fc_ess)
sd_ess <- sd(fc_ess)
dist_ess <- data.frame(
  x=seq(from=min(avgpool$fc),to=max(avgpool$fc), by=0.01)
)
dist_ess$p <- dnorm(dist_ess$x,mean=mean_ess, sd=sd_ess)
dist_ess$type<-"Essential"


#### Distribution of Dispensable
fc_disp <- avgpool$fc[avgpool$genecat=="Dispensable"]
mean_disp <- mean(fc_disp)
sd_disp <- sd(fc_disp)
dist_disp <- data.frame(
  x=seq(from=min(avgpool$fc),to=max(avgpool$fc), by=0.01)
)
dist_disp$p <- dnorm(dist_disp$x,mean=mean_disp, sd=sd_disp)
dist_disp$type<-"Dispensable"


#### Produce plot of all pools -- 1/d
if(FALSE){
  avgpool$genecat <- factor(avgpool$genecat, levels=c("Dispensable","Essential","Slow growers","Unstudied"))
  ggplot(avgpool, aes(fc, 1/sd, color=genecat)) + 
    geom_point()+
    geom_line(data=rbind(dist_disp, dist_ess),aes(x,p,color=type))+
    xlab("RGR") +
    theme_bw()+
    theme(legend.position = "none") +
    scale_color_manual(values = c("chartreuse4", "red", "dodgerblue", "#999999")) #"Dispensable","Essential","Slow growers","Unstudied"
  ggsave("/corgi/websites/malaria_crispr2024/plots_crispr/composite_histograms.pdf", width = 10, height = 10)
}

#### Produce plot of all pools -- rank vs fc
if(FALSE){
  avgpool$genecat <- factor(avgpool$genecat, levels=c("Dispensable","Essential","Slow growers","Unstudied"))
  #avgpool$fc/sd
  ggplot(avgpool, aes(rank(fc/sd), fc, color=genecat)) + 
    geom_point()+
    #  geom_line(data=rbind(dist_disp, dist_ess),aes(x,p,color=type))+
    xlab("rank(RGR/sd)") +
    ylab("RGR") +
    theme_bw()+
    theme(legend.position = "none") +
    scale_color_manual(values = c("chartreuse4", "red", "dodgerblue", "#999999")) #"Dispensable","Essential","Slow growers","Unstudied"
  #ggsave("/corgi/websites/malaria_crispr2024/plots_crispr/composite_histograms.pdf", width = 10, height = 10)
}


################### ROC curve for 
tovery <- avgpool[avgpool$genecat %in% c("Dispensable", "Essential"),]
tovery$obs <- tovery$genecat=="Dispensable"
ggplot(tovery, aes(d = obs, m = fc)) + geom_roc(labels=FALSE) +
  xlab("FPR")+ylab("TPR")+
  theme_bw()+
  theme(legend.position = "none") +
  geom_abline(slope=1, intercept=0,linetype=3)

ggsave("/corgi/websites/malaria_crispr2024/plots_crispr/composite_roc.pdf", width = 4, height = 3)




################## 3-way separation of composite phenotype

#### Logistic regression to determine gene type
totrain <- avgpool[avgpool$genecat %in% c("Dispensable","Essential","Slow growers"),]
totrain$is_ess <- totrain$genecat=="Essential"
totrain$is_disp <- totrain$genecat=="Dispensable"
logistic_model_ess <- glm(is_ess ~ fc, data = totrain, family = "binomial")
logistic_model_disp <- glm(is_disp ~ fc, data = totrain, family = "binomial")

dist_logistic <- data.frame(fc=seq(from=min(avgpool$fc),to=max(avgpool$fc), by=0.01))
dist_logistic$p_ess <- predict(logistic_model_ess, dist_logistic, type = "response")
dist_logistic$p_disp <- predict(logistic_model_disp, dist_logistic, type = "response")

cutoff_p <- 0.5
cutoff_fc_ess <- max(dist_logistic$fc[dist_logistic$p_ess>cutoff_p])
cutoff_fc_disp <- min(dist_logistic$fc[dist_logistic$p_disp>cutoff_p])

dist_logistic$scaled_p_ess <- dist_logistic$p_ess*max(1/avgpool$sd)
dist_logistic$scaled_p_disp <- dist_logistic$p_disp*max(1/avgpool$sd)

avgpool$genecat <- factor(avgpool$genecat, levels=c("Dispensable","Essential","Slow growers","Unstudied"))


df_composite <- list(
  model=logistic_model,
  dist_logistic=dist_logistic,
  avgpool=avgpool,
  cutoff_fc_ess=cutoff_fc_ess,
  cutoff_fc_disp=cutoff_fc_disp
)

ggplot(df_composite$avgpool, aes(fc, 1/sd, color=genecat)) + 
  geom_point()+
  geom_line(data=df_composite$dist_logistic,aes(fc,scaled_p_ess),color="red")+
  geom_line(data=df_composite$dist_logistic,aes(fc,scaled_p_disp),color="chartreuse4")+
  geom_vline(xintercept = df_composite$cutoff_fc_ess, color="red")+
  geom_vline(xintercept = df_composite$cutoff_fc_disp, color="chartreuse4")+
  xlab("RGR") +
  theme_bw()+
  theme(legend.position = "none") +
  scale_color_manual(values = c("chartreuse4", "red", "dodgerblue", "#999999")) #"Dispensable","Essential","Slow growers","Unstudied"
ggsave("/corgi/websites/malaria_crispr2024/plots_crispr/composite_histograms_with_cutoff.pdf", width = 10, height = 10)


saveRDS(df_composite, "/corgi/websites/malaria_crispr2024/composite.rds")

df_composite$cutoff_fc_ess

avgpool$gene[avgpool$genecat %in% c("Unstudied") & avgpool$fc < df_composite$cutoff_fc_ess]
avgpool$gene[avgpool$genecat %in% c("Unstudied") & avgpool$fc > df_composite$cutoff_fc_disp]

totrain$is_ess <- totrain$genecat=="Essential"
totrain$is_disp <- totrain$genecat=="Dispensable"


################################################################################
################## Fig xxx ROC curve for each screen ###########################
################################################################################


listpools <- c(
  "cr_2024march_half1",
  "cr_2024march_p1",
  #"cr_2024march_p2", #no essentials!
  "cr_2024march_p12"
)

##### Plot each separately
if(FALSE){
  listplot <- list()
  for(curpool in listpools){
    print(curpool)
    
    grstats <- all_grstats[[curpool]]$volcano$`NP BL6`
    tovery <- grstats[grstats$genecat %in% c("Dispensable", "Essential"),]
    tovery$obs <- (tovery$genecat=="Dispensable")+0
  
    onep <- ggplot(tovery, aes(d = obs, m = fc)) + geom_roc(labels=FALSE) + 
      xlab(paste("FPR",pools_renamed[[curpool]]))+
      ylab("TPR")+
      geom_abline(slope=1, intercept=0,linetype=3)
    
    listplot[[curpool]] <- onep
  }
  totp <- egg::ggarrange(plots=listplot)
  totp
  
  ggsave(plot = totp, "/corgi/websites/malaria_crispr2024/plots_crispr/roc_perscreen.pdf", width = 4, height = 8)
}

####### Plot all in one
alldat <- NULL
allauc <- NULL
for(curpool in listpools){
  grstats <- all_grstats[[curpool]]$volcano$`NP BL6`
  tovery <- grstats[grstats$genecat %in% c("Dispensable", "Essential"),]
  tovery$obs <- (tovery$genecat=="Dispensable")+0
  tovery$screen <- pools_renamed[[curpool]]
  
  #tovery <- tovery[order(tovery$genecat),]
  tovery <- tovery[order(tovery$fc),]
  rocob <- pROC::roc(as.integer(tovery$genecat=="Dispensable"),tovery$fc)
  print(curpool)
  print(rocob$auc)
  allauc <- rbind(allauc, data.frame(screen=pools_renamed[[curpool]], auc=pROC::auc(rocob)))
  #print(pROC::auc(rocob))
  
  alldat <- rbind(alldat,tovery)
}

alldat$screen <- factor(alldat$screen, levels=c("48x","96x","192x"))
ggplot(alldat, aes(d = obs, m = fc, color=screen)) + geom_roc(labels=FALSE) + 
  xlab(paste("FPR"))+
  ylab("TPR")+
  theme_bw()+
  #theme(legend.position = "none") +
  geom_abline(slope=1, intercept=0,linetype=3)
ggsave("/corgi/websites/malaria_crispr2024/plots_crispr/roc_perscreen_inone.pdf", width = 4, height = 3)


#Pool of 48X = black, RGB (0,0,0) Hex: #000000
#  Pool of 96X = dark grey, RGB (77,77,77) Hex: #4D4D4D
#  Pool of 192X = light grey, RGB (179,179,179) Hex: #B3B3B3

allauc$screen <- factor(allauc$screen, levels=c("48x","96x","192x"))
ggplot(allauc, aes(screen, auc)) + 
  theme_bw()+
  theme(legend.position = "none") +
  geom_bar(stat="identity")
ggsave("/corgi/websites/malaria_crispr2024/plots_crispr/roc_perscreen_barplot.pdf", width = 4, height = 3)


  
  
##################################################################################
# Fig xxxx. Rank screen plots ##
##################################################################################

  


#### Produce plot of all pools -- rank vs fc
if(FALSE){
  
  
  listpools <- c(
    "cr_2023march_screen",
    "cr_2024march_half1",
    "cr_2024march_p1",
    #"cr_2024march_p2", #no essentials!
    "cr_2024march_p12"
  )
  
  
  listplot <- list()
  for(curpool in listpools){
    print(curpool)
    
    grstats <- all_grstats[[curpool]]$volcano$`NP BL6`
    
    
    grstats$genecat <- factor(grstats$genecat, levels=c("Dispensable","Essential","Slow growers","Unstudied"))
    onep <- ggplot(grstats, aes(rank(fc/sd), fc, color=genecat)) + 
      geom_point()+
      #  geom_line(data=rbind(dist_disp, dist_ess),aes(x,p,color=type))+
      xlab(paste("rank(RGR/sd)", pools_renamed[[curpool]])) +
      ylab("RGR") +
      theme_bw()+
      theme(legend.position = "none") +
      scale_color_manual(values = c("chartreuse4", "red", "dodgerblue", "#999999")) #"Dispensable","Essential","Slow growers","Unstudied"

    listplot[[curpool]] <- onep
  }
  totp <- egg::ggarrange(plots=listplot)
  totp
  
#  ggsave(plot = totp, "/corgi/websites/malaria_crispr2024/plots_crispr/roc_perscreen.pdf", width = 4, height = 8)
}


#ggsave("/corgi/websites/malaria_crispr2024/plots_crispr/composite_histograms.pdf", width = 10, height = 10)

  


##################################################################################
# Fig xxxx. Normalized abundance over time, per gRNA; all of them, dirty format ##
##################################################################################


#for(curpool in c("cr_2024march_half1","cr_2024march_p1","cr_2024march_p12")){ #,"cr_2024march_p2""
for(curpool in c("cr_2024march_half1")){
  print(curpool)
  grstats <- all_timecourses[[curpool]]$`Count/ControlCount`
  list_all_plot <- list()
  for(curgene in unique(grstats$gene)){
    
    p1 <- ggplot(grstats[grstats$gene==curgene,], aes(x=day, y=y, linetype=grna, color=grna, group=paste(grna,mouse_ref))) + 
      geom_line(color="black")+
      #scale_color_manual(values = c("black","black"))+
      theme_bw()+
      theme(legend.position = "none")+
      #theme(axis.title.x = element_text(colour = geneinfo$xlabcol[geneinfo$gene==curgene]))+
      xlab(curgene)+
      ylab("Count/Disp. gene count")
    #p1
    list_all_plot[[curgene]] <- p1
  }
  ptot <- egg::ggarrange(plots=list_all_plot, ncol=1)
  ggsave(sprintf("/corgi/websites/malaria_crispr2024/plots_crispr/all_lineplot %s.pdf",curpool), width = 3, height = 2*length(list_all_plot), plot = ptot, limitsize=FALSE)
}  


  
  
  
################################################################################
################# Fig xxx. "Volcano" plots #####################################
################################################################################



  
for(current_pool in c("cr_2023march_screen","cr_2024march_half1","cr_2024march_p1","cr_2024march_p12")){ #,"cr_2024march_p2"
  
  highlight_genes <- c("PBANKA_1037800", "PBANKA_1401600", "PBANKA_0817800")   #PBANKA_1322400 is a slowgrower
  
  if(current_pool=="cr_2023march_screen"){
    highlight_genes <- c()
  }
  
  #current_pool <- "cr_2024march_half1"
  print(current_pool)
  thecond <- "NP BL6"
  grstats <- all_grstats[[current_pool]]
  toplot <- grstats$volcano[[thecond]]
    print(dim(toplot))
  
  toplot$y <- 1/toplot$sd
  toplot$pool <- current_pool
  
  yname <- paste("inverse s.d.")
  toplot$genecat <- factor(toplot$genecat, levels=c("Dispensable","Essential","Slow growers","Unstudied"))
  ggplot(toplot, aes(fc, y, label=gene, color=genecat)) + 
    geom_point() + 
    xlab("RGR") + 
    ylab(yname) +
    scale_color_manual(values = c("chartreuse4", "red", "dodgerblue", "#999999"))+ #"Dispensible","Essential","Slow growers","Unstudied"
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    ggrepel::geom_text_repel(data=toplot[toplot$gene %in% highlight_genes,],min.segment.length=0,nudge_x = 0.5, nudge_y = 0.3, show.legend = FALSE)
  ggsave(sprintf("/corgi/websites/malaria_crispr2024/plots_crispr/volcano %s.pdf",current_pool), width = 7, height = 6)
}

#1322400 should be red, essential





################################################################################
# Fig xxxx. Normalized abundance over time, per gRNA; for genes to highlight ###
################################################################################


#curpool <- "cr_2024march_half1"
for(curpool in c("cr_2024march_half1","cr_2024march_p1","cr_2024march_p12")){ #,"cr_2024march_p2""
  print(curpool)
  grstats <- all_timecourses[[curpool]]$`Count/ControlCount`
  list_all_plot <- list()
  for(curgene in intersect(highlight_genes,grstats$gene)){
    
    p1 <- ggplot(grstats[grstats$gene==curgene,], aes(x=day, y=y, linetype=grna, color=grna, group=paste(grna,mouse_ref))) + 
      geom_line(color="black")+
      #scale_color_manual(values = c("black","black"))+
      theme_bw()+
      theme(legend.position = "none")+
      #theme(axis.title.x = element_text(colour = geneinfo$xlabcol[geneinfo$gene==curgene]))+
      xlab(curgene)+
      ylab("Count/Disp. gene count")
    #p1
    list_all_plot[[curgene]] <- p1
  }
  ptot <- egg::ggarrange(plots=list_all_plot, nrow=1)
  ggsave(sprintf("/corgi/websites/malaria_crispr2024/plots_crispr/highlight_lineplot %s.pdf",curpool), width = 6, height = 2, plot = ptot, limitsize=FALSE)
}  



  
  



################################################################################
################# Fig 5. Relative abundance over time, per gRNA ################
################################################################################



#curpool <- "cr_2024march_half1"
listpools <- c("cr_2023march_screen","cr_2024march_half1")
for(curpool in listpools){

  print(curpool)
  
  grstats <- all_timecourses[[curpool]]$`Count/AllCount`
  
  #Figure out where in a grid all plots should go
  geneinfo <- unique(grstats[,c("gene","genecat")])
  geneinfo <- geneinfo[order(geneinfo$genecat),]
  geneinfo$i <- 1:nrow(geneinfo)
  geneinfo <- merge(geneinfo,sqldf::sqldf("select min(i) as mini, genecat from geneinfo group by genecat"))
  geneinfo$grid_y <- geneinfo$i-geneinfo$mini+1
  geneinfo$grid_x <- as.integer(factor(geneinfo$genecat))
  
  #Figure out coloring based on gene category  
  geneinfo$xlabcol <- "gray" #Other etc
  geneinfo$xlabcol[geneinfo$genecat == "Essential"] <- "red"
  geneinfo$xlabcol[geneinfo$genecat == "Dispensable"] <- "chartreuse4"
  geneinfo$xlabcol[geneinfo$genecat == "Slow"] <- "dodgerblue"
  geneinfo$xlabcol[geneinfo$genecat == "Slow growers"] <- "dodgerblue"
  
  pdf(sprintf("/corgi/websites/malaria_crispr2024/plots_crispr/gridpanel_per_gene %s.pdf",curpool), width = 7, height = max(geneinfo$grid_y)*1.1)
  grid.newpage()
  vp <- viewport(layout = grid.layout(max(geneinfo$grid_y), max(geneinfo$grid_x)))
  pushViewport(vp)
  for(curgene in unique(grstats$gene)){
    
    p1 <- ggplot(grstats[grstats$gene==curgene,], aes(x=day, y=y, linetype=grna, color=grna, group=paste(grna,mouse_ref))) + 
      geom_line(color="black")+
      #scale_color_manual(values = c("black","black"))+
      theme_bw()+
      theme(legend.position = "none")+
      theme(axis.title.x = element_text(colour = geneinfo$xlabcol[geneinfo$gene==curgene]))+
      xlab(curgene)+
      ylab("Rel.count")
    
    p1
    
    print(p1, vp = viewport(
      layout.pos.row = geneinfo$grid_y[geneinfo$gene==curgene], 
      layout.pos.col = geneinfo$grid_x[geneinfo$gene==curgene]))
    
  }
  dev.off()

}



  