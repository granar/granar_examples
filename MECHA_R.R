
MECHA_R <- function(path = "granar_examples/MECHAv4_light.py"){

  tic <- Sys.time()
  message("Launch MECHA")
  MECHA <- try(source_python(path), silent = T)
  if(is.null(MECHA)){MECHA <- c(0)}
  if(str_detect(MECHA[1], "Error ")){
    warning("error NaN: kr_tot <- NA")
    kr_tot <- c(NA, NA, NA)
  }
  if(MECHA[1] == "Error in py_run_file_impl(file, local, convert) : \n  IndexError: index 25 is out of bounds for axis 0 with size 25\n"){
    warning("Undersized matrix of connected cell.")
    kr_tot <- c(NA, NA, NA)
  }
  print(Sys.time()-tic)
  return(kr_tot)
}