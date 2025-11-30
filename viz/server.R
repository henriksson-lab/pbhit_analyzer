library(plotly)
library(Cairo)
library(naturalsort)

options(shiny.usecairo=T)



if(TRUE){
  options(shiny.sanitize.errors = FALSE)
  # logging level DEBUG
  logging::basicConfig(level = 10)
  # write logging output to the stderr file
  logging::addHandler(logging::writeToFile, logger = '', file = stderr())
}



if(FALSE){
  #To run this app
  library(shiny)
  runApp(".")
}



server <- function(input, output, session) {

  observeEvent(input$grstats_pool,{
    current_pool <- input$grstats_pool
    grstats <- all_grstats[[current_pool]]
    
    updateSelectizeInput(session, 'grstats_scatter', choices = names(grstats$scatterplot), server = TRUE)
    
    grstats <- all_timecourses[[current_pool]]
    all_gr_type <- names(grstats)
    updateSelectizeInput(session, 'grstats_gene', choices = c("",unique(grstats[[all_gr_type[1]]]$gene)), server = TRUE)
    updateSelectizeInput(session, 'grstats_units', choices = all_gr_type, server = TRUE)

  })
  


  ################################################################################
  ########### Sample metadata ####################################################
  ################################################################################

  output$plotSamplemetaUmap <- renderPlot(height=700, {

    current_pool <- input$samplemeta_pool
    #print(current_pool)
    
    samplemeta <- all_samplemeta[[current_pool]]
    coverage_stat <- all_coverage_stat[[current_pool]]
    
    #print(samplemeta)
    
    samplemeta$day <- sprintf("d%s", samplemeta$day)
    p1 <- ggplot(samplemeta, aes(umap1,umap2,color=mouse_ref))+geom_point()
    p2 <- ggplot(samplemeta, aes(umap1,umap2,color=day))+geom_point()
    p3 <- ggplot(samplemeta, aes(umap1,umap2,color=is_input))+geom_point()
    p4 <- ggplot(samplemeta, aes(umap1,umap2,color=genotype))+geom_point()
    p5 <- ggplot(samplemeta, aes(umap1,umap2,color=primed))+geom_point()
    p6 <- ggplot(samplemeta, aes(umap1,umap2,color=total_count))+geom_point()
    ptot <- p1/p2|p3/p4|p5/p6 #|p7
    
    coverage_stat <- coverage_stat[order(coverage_stat$cnt, decreasing = TRUE),]
    coverage_stat$grna <- factor(coverage_stat$grna, levels=coverage_stat$grna)
    covstatplot <- ggplot(coverage_stat, aes(grna,cnt,color=genecat)) + 
      geom_point() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      scale_color_discrete(name = "PlasmoGEM Phenotype")

    #Natural order on genewiz column    
    coverage_stat <- coverage_stat[naturalsort::naturalorder(coverage_stat$genewiz),]  #R install naturalsort
    coverage_stat$genewiz <- factor(coverage_stat$genewiz, levels=naturalsort::naturalsort(unique(coverage_stat$genewiz)))
    
    wellplot <- ggplot(coverage_stat, aes(genewiz, ligationwell, color=log10(cnt))) + 
      geom_point() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    ptot/(covstatplot|wellplot)

  })

  ################################################################################
  ########### For a single screen - volcano ######################################
  ################################################################################

  get_current_volcano <- function(){
    current_pool <- input$grstats_pool
    print(current_pool)
    grstats <- all_grstats[[current_pool]]
    thecond <- "NP BL6"                   
    
    if(thecond %in% names(grstats$volcano)){
      grstats$volcano[[thecond]]
    } else {
      data.frame()
    }
  }
  
  output$plot_grstats_volcano <- renderPlotly({
    
    
    thecond <- "NP BL6" 
    toplot <- get_current_volcano()
    
    if(input$grstats_y=="-Log10 p, different from control genes"){
      toplot$y <- toplot$logp
      yname <- paste("-log10 pval")
    } else {
      toplot$y <- 1/toplot$sd
      yname <- paste("inverse s.d.")
    }
    
    toplot$genecat <- factor(toplot$genecat, levels=c("Dispensable","Essential","Slow growers","Unstudied"))

    if(nrow(toplot)>0){
      
      #Fix names of genes
      #if(any(str_length(toplot$gene)==str_length("PBANKA_071020"))){
      #  toplot$gene <- paste0(toplot$gene, "0")
      #}
      
      
      if(input$grstats_show_gene_name){
        theplot <- ggplot(toplot, aes(fc, y, label=gene, color=genecat)) + 
          geom_point(color="gray") + 
          geom_text() +
          xlab(paste("RGR")) + 
          ylab(yname) +
          scale_color_discrete(name = "PlasmoGEM Phenotype")
        
      } else {
        
        #For hover, add gene symbol
        toplot <- merge(toplot,geneinfo,all.x=TRUE)
        toplot$gene <- paste0(toplot$gene, "\nName: ", toplot$geneName, "\nDescription: ", toplot$geneDesc)
        
        
        theplot <- ggplot(toplot, aes(fc, y, label=gene, color=genecat)) + 
          geom_point() + 
          xlab(paste("RGR")) + 
          ylab(yname) +
          scale_color_manual(values = c("chartreuse4", "red", "dodgerblue", "#999999"), name = "PlasmoGEM Phenotype") #"Dispensible","Essential","Slow growers","Unstudied"
        #https://sape.inf.usi.ch/quick-reference/ggplot2/colour
      }
    } else {
      theplot <- ggplot() + theme_void()
    }
    theplot  %>% ggplotly(source="plot_grstats_volcano") %>% event_register("plotly_click")
  })
  
  ## Callback for clicking on the top-left plot
  observeEvent(
    eventExpr = event_data("plotly_click", source = "plot_grstats_volcano"),
    handlerExpr = {
      
      toplot <- get_current_volcano()

      event_data <- event_data("plotly_click", source = "plot_grstats_volcano")

      if(input$grstats_y=="-Log10 p, different from control genes"){
        toplot$y <- toplot$logp
      } else {
        toplot$y <- 1/toplot$sd
      }
      
      toplot$dx <- (toplot$fc - event_data$x)
      toplot$dy <- (toplot$y - event_data$y)
      toplot$dist <- toplot$dx**2 + toplot$dy**2
      toplot <- toplot[order(toplot$dist, decreasing = FALSE),]
      
      clicked_gene <- toplot$gene[1] 
      clicked_gene <- str_split_fixed(clicked_gene," ",3)[1] ## remove gene name in hover
      updateSelectInput(session, "grstats_gene", selected = clicked_gene)
    }
  )  
  
  
################################################################################
########### Comparison of two screens - scatter or volcano #####################
################################################################################

  get_current_scatter <- function(){
    current_pool <- input$grstats_pool
    grstats <- all_grstats[[current_pool]]
    thecond <- input$grstats_scatter
    
    if(thecond %in% names(grstats$scatterplot)){
      grstats$scatterplot[[thecond]]
    } else {
      data.frame()
    }
  }
  
  
  
  output$plot_grstats_scatterplot <- renderPlotly({
    
    current_pool <- input$grstats_pool
    grstats <- all_grstats[[current_pool]]
    thecond <- input$grstats_scatter
    
    represent_as <- input$grstats_scatter_type  ##### possibly store xlab and ylab names
    
    
    if(thecond %in% names(grstats$scatterplot)){
      toplot <- grstats$scatterplot[[thecond]]
      thecond2 <- str_split_fixed(thecond," / ",2)  #hopefully works
      cond1 <- thecond2[1]
      cond2 <- thecond2[2]
      
      
      if(represent_as=="RGR scatter plot"){
        
        fc_range <- range(c(toplot$fc1, toplot$fc2))
        theplot <- ggplot(toplot, aes(fc1,fc2, label=gene, color=genedesc)) + geom_point(color="gray") + geom_text()+#size=1) +
          xlab(paste("RGR",cond1)) + ylab(paste("RGR",cond2)) +
          xlim(fc_range[1], fc_range[2]) + ylim(fc_range[1], fc_range[2])
        
      } else {
        
        theplot <- ggplot(toplot, aes(diff_fc, diff_log_p, label=gene, color=genedesc)) + 
          geom_point(color="gray") + 
          geom_text() +
          xlab(paste("RGR",thecond)) + 
          ylab(paste("-log10 pval",thecond))
        
      }
      
    } else {
      print("missing comparison cond")
      theplot <- ggplot() + theme_void()
    }
    theplot %>% ggplotly(source="plot_grstats_scatterplot") %>% event_register("plotly_click")
  })
  
  
  
  ## Callback for clicking on the top-right plot
  observeEvent(
    eventExpr = event_data("plotly_click", source = "plot_grstats_scatterplot"),
    handlerExpr = {
      event_data <- event_data("plotly_click", source = "plot_grstats_scatterplot")
      #print(event_data)
      clicked_gene <- get_current_scatter()$gene[event_data$pointNumber+1]  #plotly seems to do 0-indexing
      updateSelectInput(session, "grstats_gene", selected = clicked_gene)
      }
  )  
  
  

  ################################################################################
  ########### GRstats - timecourse ###############################################
  ################################################################################
  
  
  plotTC <- function(
    grstats,
    grstats_avg_grna, grstats_avg_mouse, grstats_avg_genotype, grstats_avg_treatment,
    grstats_gene, grstats_colorby
  ){
    
    ########### Average together based on user input
    
    if(grstats_avg_grna){
      grstats <- sqldf::sqldf(
        "select day, avg(y) as y, gene, primed, genotype, mouse_ref from grstats group by mouse_ref, gene, day, primed, genotype")
      grstats$grna <- paste(grstats$gene,"*",sep="")
    }
    
    if(grstats_avg_mouse){
      grstats <- sqldf::sqldf(
        "select day, avg(y) as y, gene, grna, primed, genotype from grstats group by gene, grna, day, primed, genotype")
      grstats$mouse_ref <- "m*"
    }
    
    if(grstats_avg_genotype){
      grstats <- sqldf::sqldf(
        "select day, avg(y) as y, gene, grna, primed from grstats group by gene, grna, day, primed")
      grstats$primed <- "g*"
    }
    
    if(grstats_avg_treatment){
      grstats <- sqldf::sqldf(
        "select day, avg(y) as y, gene, grna, genotype from grstats group by gene, grna, day, genotype")
      grstats$primed <- "t*"
    }
    
    grstats$group <- paste(grstats$grna, grstats$mouse_ref, grstats$primed, grstats$genotype)
    
    ######## Only show one gene, optionally
    current_gene <- input$grstats_gene
    if(current_gene!=""){
      grstats <- grstats[grstats$gene==current_gene,,drop=FALSE]
    }
    
    
    ######## Decide coloring strategy
    grstats$colorby <- grstats$mouse_ref
    if(grstats_colorby=="Gene"){
      grstats$colorby <- grstats$gene
    }
    if(grstats_colorby=="Genotype"){
      grstats$colorby <- grstats$genotype
    }
    if(grstats_colorby=="Treatment"){
      grstats$colorby <- grstats$primed
    }
    if(grstats_colorby=="Genotype+Treatment"){
      grstats$colorby <- paste(grstats$genotype,grstats$primed)
    }
    if(grstats_colorby=="Genetic construct"){
      grstats$colorby <- paste(grstats$grna)
    }
    if(grstats_colorby=="Genotype+Treatment+Genetic construct"){
      grstats$colorby <- paste(grstats$genotype,grstats$primed, grstats$grna)
    }
    
    ######## The actual plotting
    ggplotly(ggplot(grstats,aes(x=day,y=y, group=group, color=colorby, text=group)) + 
               geom_line()+
               xlab("Day")+
               ylab("Relative abundance")+
               ggtitle("") +
               scale_color_discrete(name = " ")
             , tooltip = c("x", "y", "color", "text", "group"))
  }
  
  
  
  
  ##############
  output$plot_grstats_tcplot <- renderPlotly({
    
    current_pool <- input$grstats_pool
    grstats <- all_timecourses[[current_pool]]
    
    ## Pick the right unit to show
    grstats <- grstats[[input$grstats_units]]  
    #print(head(grstats))
    
    plotTC(
      grstats,
      input$grstats_avg_grna, input$grstats_avg_mouse, FALSE, FALSE, #input$grstats_avg_genotype, input$grstats_avg_treatment,
      input$grstats_gene, input$grstats_colorby
    )
      
  })
  
  
  
    
}

