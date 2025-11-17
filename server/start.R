library(shiny)

if(FALSE){
  #install.packages("reactlog")
  # cmd + Fn + 3
  library(reactlog)
  options(shiny.reactlog = TRUE)
  reactlog_enable()
}

runApp(".", port=10010, host="0.0.0.0")

if(FALSE){
  shiny::reactlogShow()
}