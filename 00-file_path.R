library(tools)

# support function which resolves filepath by extension
pp <- function(fname){
  fdir <- switch (file_ext(fname),
    "txt" = "ICELOGO",
    "pdf" = "FIGURES",
    "png" = "FIGURES",
    "csv" = "TABLES"
  )
  paste0("./",fdir,"/",fname)
}
