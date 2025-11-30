library(shiny)

if(FALSE){
  #install.packages("reactlog")

  #If you run this code, you will get improved debug info to help develop the code further
  library(reactlog)
  options(shiny.reactlog = TRUE)
  reactlog_enable()
}

runApp(".", port=8080, host="0.0.0.0")

if(FALSE){
  shiny::reactlogShow()
}
