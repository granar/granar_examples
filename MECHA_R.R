
MECHA_R <- function(path = "granar_examples/MECHAv4_light.py", Kx = F){
  
  tic <- Sys.time()
  message("Launch MECHA")
  MECHA <- try(source_python(path), silent = T)
  if(is.null(MECHA)){MECHA <- c(0)}
  if(str_detect(MECHA[1], "Error ")){
    message("error NaN: kr_tot <- NA")
    kr_tot <- c(NA, NA, NA)
  }
  if(MECHA[1] == "Error in py_run_file_impl(file, local, convert) : \n  IndexError: index 25 is out of bounds for axis 0 with size 25\n"){
    message("Undersized matrix of connected cell.")
    kr_tot <- c(NA, NA, NA)
  }
  print(Sys.time()-tic)
  if(Kx){
    K <- data.frame(kr_tot,K_xyl_spec)
    return(K)
  }else{
    K <- kr_tot
    return(K)
    }
}