update_MECHA_inputs <- function(MECHA_inputs_list){
  # SET PATHS
  MECHA_inputs_list$Geometry$param$path$value                                     <- "Daniela/Arabido4bis.xml" # XML Cross section file
  # MECHA_inputs_list$General$properties$Output$path                                <- output_path               # Output folder 
  MECHA_inputs_list$Geometry$param$Plant$value                                    <- ""
  MECHA_inputs_list$BC$properties$path_scenarios$Output$path                      <- ""
  MECHA_inputs_list$Hydraulics$param$path_hydraulics$Output$path                  <- ""
  
  # ============================================================================ #
  # Set Transient as 0 so that we only simulate steady-state
  MECHA_inputs_list$Hormones$param$Time_information$Transient_sim$value           <- as.integer(0)
  # ============================================================================ #
  # Set observation point
  MECHA_inputs_list$Geometry$param$Printrange$Print_layer$value                   <- iLayer
  # ============================================================================ #
  # Paraview outputs
  MECHA_inputs_list$General$properties$Paraview$WallFlux                          <- 1
  MECHA_inputs_list$General$properties$Paraview$MembraneFlux                      <- 1
  MECHA_inputs_list$General$properties$Paraview$PlasmodesmataFlux                 <- 1
  # ============================================================================ #
  # PLASMODESMATA HYDRAULIC CONDUCTANCE
  MECHA_inputs_list$Hydraulics$param$Kplrange$Kpl$value                           <- 3.5e-5              
  MECHA_inputs_list$Hydraulics$param$Kplrange$Kpl$PCC_factor                      <- as.double(1)        
  MECHA_inputs_list$Hydraulics$param$Kplrange$Kpl$PPP_factor                      <- as.double(1)
  MECHA_inputs_list$Hydraulics$param$Kplrange$Kpl$PST_factor                      <- as.double(1)
  MECHA_inputs_list$Hydraulics$param$Kplrange$Kpl$endo_in_factor                  <- as.double(10)
  MECHA_inputs_list$Hydraulics$param$Kplrange$Kpl$endo_out_factor                 <- as.double(1e-4) 
  # ============================================================================ #
  # PLASMODESMATA FREQUENCY
  MECHA_inputs_list$Hydraulics$param$Fplxheight_stele_sieve$value                 <- as.double(35000)  # Default : 350000
  MECHA_inputs_list$Hydraulics$param$Fplxheight_stele_comp$value                  <- as.double(35000)  # Default : 350000
  # ============================================================================ #
  # CELL WALL CONDUCTIVITIES
  MECHA_inputs_list$Hydraulics$param$kw_barrier_range$kw_barrier$Casp             <- as.double(1e-8)
  MECHA_inputs_list$Hydraulics$param$kw_barrier_range$kw_barrier$Sub_out          <- as.double(1e-16)
  MECHA_inputs_list$Hydraulics$param$kw_barrier_range$kw_barrier$Sub_in           <- as.double(1e-16)
  # ============================================================================ #
  # SET BOUNDARY CONDITIONS
  MECHA_inputs_list$BC$properties$BC_xyl_range$BC_xyl[[2]][["flowrate_prox"]]     <- 0.005                # Xylem BC as flowrate 
  MECHA_inputs_list$BC$properties$BC_sieve_range$BC_sieve[[2]][["osmotic"]]       <- as.double(-1.26e4)   # Phloem BC as osmotic (-1.26e4) or not ('nan')
  MECHA_inputs_list$BC$properties$Psi_cell_range$Psi_cell[[2]][["Os_hetero"]]     <- as.integer(1)        # other cells as osmotic (1) or not (0)
  # ============================================================================ #
  # ELONGATION RATE
  MECHA_inputs_list$BC$properties$Elong_cell_range$Elong_cell[[2]][["midpoint_rate"]] <- as.double(0.034) # Ellongation rate. Set 0 if no.
  
  return(MECHA_inputs_list)
}


