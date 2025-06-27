
create_inputs_object <- function(){
  
  General_inputs <- list(
    
  )
  # ============== #
  BC_inputs <- list(
    "Q_xyl"                  = 0.003,     # Optimized value 
    "BC_sieve_pressure_prox" = "nan",     # Set 1.08E04 for putting phloem under pressure
    "BC_sieve_osmotic"       = "nan"      # Set -1.26E04 for putting phloem under osmotic conditions
  )
  # ============== #
  Geometry_inputs <- list(
    "PDS"                    =  7.45e-5   # From litterature 
  )
  # ============== #
  Hydraulics_inputs <- list(
    "Kw"                     = 1e-7,
    "Kpl"                    = 3.5e-5,
    "Kpl_fac_in"             = 1,         # Optimized = 10
    "Kpl_fac_out"            = 1,         # Optimized = 1e-4
    "Kw_sub_in"              = 1e-16,               
    "Kw_sub_out"             = 1e-16,
    "Kw_CS"                  = 1e-16,      # Optimized = 1e-8
    "kAQP"                   = 4.3e-4     # From Ehlert et al 2009
  )
  # ============== #
  Hormones_inputs <- list(
    "Diff_PD"                = 6.2       # Optimized in 2021
    ,"Diff_PW"               = 8.62      # Optimized in 2021
    ,"Diff_X"                = 10.5       # Optimized in 2021
    ,"transient"             = 1 
  )
  # ============== #
  inputs <- list(
    "General"   = General_inputs,
    "BC"        = BC_inputs,
    "Hydraulic" = Hydraulics_inputs,
    "Geometry"  = Geometry_inputs,
    "Hormones"  = Hormones_inputs
  )
  # ============== #
  return(inputs)
}

# ============================================================================ #

write.MECHA.param <- function(input_path, # path where MECHA param files are stored
                              output_path,
                              MECHA.inputs,
                              Barriers,
                              Contact,
                              UniXwalls = F,
                              cid,
                              iLayer,
                              lengths,
                              heights,
                              washin_cc = c(1, 1, 0, 0), 
                              washout_cc = c(0, 0, 0, 0),
                              Ini_cond = 1
                              
                              
){
  
  write_BC_xml(input_path, 
               output_path, 
               BC_inputs = MECHA.inputs$BC)
  
  write_Geom_xml(input_path, 
                 Barriers, 
                 lengths, 
                 heights, 
                 Geometry_inputs = MECHA.inputs$Geometry)
  
  write_Hydraulics_xml(input_path = input_path
                       ,output_path = output_path
                       ,Barriers = Barriers
                       ,Hydraulics_inputs = MECHA.inputs$Hydraulic
                       ,Contact = Contact
                       ,UniXwalls = UniXwalls)
  
  write_Hormones_xml(input_path = input_path
                     ,Hormones_inputs = MECHA.inputs$Hormones
                     ,cid = cid
                     ,iLayer = iLayer
                     ,lengths = lengths
                     ,heights = heights
                     ,washin_cc = washin_cc
                     ,washout_cc = washout_cc
                     ,Ini_cond = Ini_cond
                     )
  
}

# ============================================================================ #

write_BC_xml <- function(input_path, 
                         output_path, 
                         BC_inputs
                         ){
  #' @param input_path Path in where BC.xml will be written
  #' @param output_path Path to the folder of different scenarios
  #' @param inouts A list object containing parameters
  #' 
  
  Q_xyl                  <- BC_inputs$Q_xyl
  BC_sieve_pressure_prox <- BC_inputs$BC_sieve_pressure_prox
  BC_sieve_osmotic       <- BC_inputs$BC_sieve_osmotic
  
  
  lines <- c("<?xml version='1.0' encoding='utf-8'?>",
             "",
             "<!--#Boundary conditions: the same total number of scenarios should be prescribed for soil, xylem and cells",
             "The first scenario is currently 'forced' as it is used to calculate the kr and STF intrinsic properties-->",
             "<properties>",
             "    <!-- Path to the folders of different scenarios, note that the listed scenarios are gathered in the same folder -->",
             "    <path_scenarios>",
             paste0("        <Output path='", output_path,"'/>"),
             "    </path_scenarios>",
             "",
             "    <!-- Soil water potentials",
             "        Units: hPa",
             "        Left and right osmotic potential allow creating a continuous gradient of osmotic potentials accross root surface and cortex (if no Casparian strip was developped in the exodermis)",
             "        osmotic symmetry controls gradients of apoplastic osmotic potentials (1: linear left-right gradient; 2: linear central symmetry, where 'right' stands for the inner part of the cortex)-->",
             "    <Psi_soil_range>",
             "        <osmotic_convection flag='0'/>",
             "        <osmotic_diffusivity value='0'/>",
             "        <Psi_soil pressure_left='0.0'   pressure_right='0.0'     osmotic_left='0.0'	osmotic_right='0.0'	osmotic_symmetry='1'	osmotic_shape='1'/>",
             "        <Psi_soil pressure_left='0.0'   pressure_right='0.0'     osmotic_left='0.0'	osmotic_right='0.0'	osmotic_symmetry='1'	osmotic_shape='1'/>",
             "    </Psi_soil_range>",
             "",
             "    <!-- Xylem boundary conditions",
             "        Units: hPa (pressure) or cm^3/d (flow rate, positive when transpiring)",
             "        deltaP is a boundary pressure condition set as compared to the xylem equilibrium pressure (hPa). This boundary condition can only be used after a no-flow BC scenario.",
             "        In case of flux BC, there can either be xylem outflow on the proximal side only (then flowrate_dist should be null), or on both sides (positive being towards the shoot)",
             "        -1.5E3 osmotic: Reference treatment (Enns et al., 2000)",
             "        -2.2E3 osmotic: KNO3 treatment to help guttation (Enns et al., 2000)",
             "        100% transpiration (4 mm from barrier) = 0.00507024 * 15 * 30 mm / 21 mm = 0.108648",
             "        The xylem osmotic potential propagates throughout stelar walls -->",
             "    <BC_xyl_range>",
             "        <osmotic_convection flag='0'/>",
             "        <osmotic_diffusivity value='0'/>",
             "        <BC_xyl pressure_prox='-5.0E3'   deltaP_prox='nan'	pressure_dist='nan'	    flowrate_prox='nan'			flowrate_dist='nan'			osmotic_endo='0.0'	osmotic_xyl='0.0'	osmotic_symmetry='1'	osmotic_shape='1'/> <!-- Used to calculate kr -->",
             paste0("        <BC_xyl pressure_prox='nan'	  deltaP_prox='nan'	pressure_dist='nan'	    flowrate_prox='", Q_xyl, "'	flowrate_dist='0.0000000'	osmotic_endo='0.0'	osmotic_xyl='0.0'	osmotic_symmetry='1'	osmotic_shape='1'/> <!-- 4 mm default transpiration flowrate_prox='0.00507024'	flowrate_dist='0.00411957 -->"),
             "    </BC_xyl_range>",
             "",
             "    <!-- Phloem sieve tubes water potentials",
             "        Units: hPa (pressure) or cm^3/d (flow rate, positive when absorbing water from the cross-section into sieve elements)",
             "        -1.42E4 osmotic: Barley phloem (Pritchard, 1996)",
             "        1.1E4 pressure: Barley phloem, estimation (Pritchard, 1996)",
             "        -1.26E4 osmotic: Arabidopsis phloem (Hall et al., 2006) Melting point-depression",
             "        If nan pressure & deltaP & flowrate: no phloem BC",
             "        If both flowrate_*==0.0 and Barrier>0 (differentiation zone): no flux in phloem",
             "        If both flowrate_*==0.0 and Barrier==0 (elongation zone): phloem source equals elongation sink, no soil water involved -->",
             "    <BC_sieve_range>",
             "        <BC_sieve pressure_prox='0.8E4'  deltaP_prox='nan'	pressure_dist='nan'	flowrate_prox='nan'	flowrate_dist='nan'  osmotic='0.0'    />",
             # "        <BC_sieve pressure_prox='nan'    deltaP_prox='nan'	pressure_dist='nan'	flowrate_prox='nan'	flowrate_dist='nan'  osmotic='-8.0E3' />",
             paste0("        <BC_sieve pressure_prox='", BC_sieve_pressure_prox, "'	deltaP_prox='nan'	pressure_dist='nan'	flowrate_prox='nan'			flowrate_dist='nan'  osmotic='", BC_sieve_osmotic, "' />"),
             "    </BC_sieve_range>",
             "",
             "    <!-- Cellular protoplasts osmotic potentials",
             "        Units: hPa",
             "        Os_cortex is the cortical cell osmotic potential",
             "        -8.0E3: Reference treatment (Enns et al., 2000)",
             "        -1.26E4: KNO3 treatment to help guttation (Enns et al., 2000)",
             "        Os_hetero can generate heterogeneity in cell osmotic potential distribution",
             "        0: Uniform",
             "        1: non-uniform, reference treatment (Enns et al., 2000)",
             "        2: non-uniform, with KNO3 treatment to help guttation (Enns et al., 2000)",
             "        3: Simple test with uniform osmotic potential in pericycle and stele at -5000 hPa)",
             "        & Reflection coefficients",
             "            s_hetero (unitless)",
             "            0: Uniform",
             "            1: Non-uniform, stele twice more permeable to solute",
             "            2: Non-uniform, cortex twice more permeable to solute",
             "            s_factor  (unitless)",
             "            0->1: multiplies all reflection coefficients -->",
             "    <Psi_cell_range>",
             "        <Psi_cell Os_cortex='0.0'    Os_hetero='0' s_hetero='0' s_factor='1.0'/>",
             "        <Psi_cell Os_cortex='-8.0E3' Os_hetero='0' s_hetero='0' s_factor='1.0'/>",
             "    </Psi_cell_range>",
             "",
             "    <!-- Cell elongation rate",
             "        'midpoint_rate' units: cm/d",
             "        0.034 cm/d: Beemster and Baskin, 1998",
             "        'side_rate_difference' units: cm/d, positive if more elongation to the right of the cross-section, negative otherwise",
             "        0.0072 cm/d: Draft Dietrich et al.",
             "        'kappa_dot' units: radian/micron/d",
             "        0.000007 radian/micron/d: .......... -->",
             "    <Elong_cell_range>",
             "        <Elong_cell midpoint_rate='0' 	 side_rate_difference='0'  /> <!-- kappa_dot='0.000007' -->",
             "        <Elong_cell midpoint_rate='0' 	 side_rate_difference='0'  />",
             "    </Elong_cell_range>",
             "",
             "    <!-- Relative volume of water in symplast and in apoplast",
             "        Units: dimensionless",
             "        Cellulosic wall volumetric water content",
             "        0.69: Extract from eucalyptus leaves (Gaff and Carr, 1961)",
             "        Cell protoplast volumetric water content",
             "        0.7: Assumed to have the density of 1 molar sucrose (Gaff and Carr, 1961) + http://homepages.gac.edu/~cellab/chpts/chpt3/table3-2.html",
             "        -->",
             "    <Water_fractions Symplast='0.7' Apoplast='0.69'/>",
             "</properties>")
  
  writeLines(lines, paste0(input_path, 'BC_optim.xml'))
}

# ==============================================================================

write_Geom_xml <- function(input_path, 
                           Barriers, 
                           lengths, 
                           heights, 
                           Geometry_inputs){
  #' @param input_path Path where Geometry_optim.xml will be written
  #' @param Barriers vector containing Barriers types in MECHA format
  #' @param lengths vector containing lengths for each virtual root segment
  #' @param heights vector containing heights for each virtual root segment
  
  PDS <- Geometry_inputs$PDS
  
  # ====================== #
  
  lines <- c(
    '<?xml version="1.0" encoding="utf-8"?>',
    '',
    '<param>',
    '    <!-- Plant type -->',
    '    <Plant value=\'Arabido\' /> <!--#Maize / Arabido / Millet / Barley /-->',
    '',
    '    <!-- Image path and properties -->',
    '    <path value=\'Daniela/Arabido4bis.xml\' />',
    '    <im_scale value="0.16" /> <!--image scale (micron per pixel)-->',
    '',
    '    <!-- Maturity level',
    '    0: No apoplastic Barriers',
    '    1: Endodermal Casparian strip (radial walls)',
    '    2: Endodermal suberization except at passage cells',
    '    3: Endodermis full suberization',
    '    4: Endodermis full suberization and exodermal Casparian strip (radial walls)',
    '    5: Endodermal and exodermal Casparian strips (radial walls) -->',
    '    <Maturityrange> <!-- All the listed barrier types will be simulated and reported in separate files "***b1", "***b2", "***b3",... -->')
  
  for(i in seq(1, length(lengths))){
    lines <- c(lines, paste0('       <Maturity Barrier="', Barriers[i], '" height="', heights[i]*1000, '" Nlayers="', round(lengths[i]/heights[i]), '"/>'))}  

  lines <- c(lines,
    '    </Maturityrange>',
    '    <Printrange>',
    '        <Print_layer value="12"/>',
    '    </Printrange>',
    '    <Xwalls value="1" /> <!-- 0: No transverse walls in the 2D simulations; 1: Transverse walls included in 2D simulations -->',
    '    <PileUp value="2" /> <!-- 0: Simulating different levels of maturity separately (2D); 1: Simulating different levels of maturity separately but with continuous boundary conditions (3D); Simultaneously solving flow across all fully connected levels of maturity (full 3D) -->',
    '',
    '    <!-- Topological info (passage cells and intercellular spaces) -->',
    '    <passage_cell_range>',
    '        <passage_cell id="-1" /> <!-- ID number of passage cells in the endodermis, ideally in front of early metaxylem vessels -->',
    '    </passage_cell_range>',
    '    <aerenchyma_range>',
    '        <aerenchyma id="-1" /> <!-- ID number of aerenchyma "cells" in the cortex -->',
    '    </aerenchyma_range>',
    '    <InterC_perim_search value="0" /> <!-- 0: no search of intercellular spaces based on their perimeter; 1: search of intercellular spaces based on their perimeter -->',
    '    <InterC_perim1 value="0" /> <!-- (micrometers) threshold for maximum perimeter of intercellular spaces in first cortical cell layers (next to endodermis) -->',
    '    <InterC_perim2 value="0" /> <!-- (micrometers) threshold for maximum perimeter of intercellular spaces in 2nd cortical cell layers (next to endodermis) -->',
    '    <InterC_perim3 value="0" /> <!-- (micrometers) ... -->',
    '    <InterC_perim4 value="0" /> <!-- (micrometers) ... -->',
    '    <InterC_perim5 value="0" /> <!-- (micrometers) ... -->',
    '    <kInterC value="0.0" /> <!-- 0 for air-filled intercellular spaces / inf for water-filled intercellular spaces -->',
    '',
    '    <!-- Geometrical info -->',
    '    <cell_per_layer cortex="nan" stele="nan"/>',
    '    <diffusion_length cortex="nan" stele="nan"/>  <!-- (microns) -->',
    '',
    '    <!-- Cell wall properties -->',
    '    <thickness value="0.2" /> <!-- Double cell wall thickness (microns) Andeme-Onzighi et al. (2002)) -->',
    '',
    '    <!-- Plasmodesmata geometrical properties -->',
    paste0('    <PD_section value="', PDS, '" /> <!-- PD open cross-section (microns^2) (Ehlers and Bel, 1999; Ehlers and Kollmann, 2001) -->'),
    '',
    '    <!-- Were xylem vessels cut in smaller pieces (debug) -->',
    '    <Xylem_pieces flag="0" /> <!-- 0: No; 1: Yes -->',
    '</param>'
  )
  
  writeLines(lines, paste0(input_path, "Geometry_optim.xml"))
}

# ==============================================================================

write_Hydraulics_xml <- function(input_path, output_path, 
                                 Barriers,
                                 Hydraulics_inputs,
                                 Contact = c(0, 0, 9999, 0), 
                                 UniXwalls = F){
  #' @param input_path Path where Hydraulics_optim.xml will be written
  #' @param output_path Path to the folders of different hydraulic scenarios
  #' @param Barriers A 4 size vector containing which type of barrier is present in each root part.0: no apoplastic barrier nor mature xylem; 1: endodermal Casparian strip; 2: suberised endodermis with passage cells; 3: fully suberised endodermis
  
  Kw          = Hydraulics_inputs$Kw                   # Default value
  Kw_CS       = Hydraulics_inputs$Kw_CS 
  Kw_sub_out  = Hydraulics_inputs$Kw_sub_out
  Kw_sub_in   = Hydraulics_inputs$Kw_sub_in
  kAQP        = Hydraulics_inputs$kAQP
  Kpl         = Hydraulics_inputs$Kpl
  Kpl_fac_in  = Hydraulics_inputs$Kpl_fac_in
  Kpl_fac_out = Hydraulics_inputs$Kpl_fac_out
  
  
  lines <- c("<?xml version='1.0' encoding='utf-8'?>",
             "<!-- Hydraulic options -->",
             "<param>",
             "    <!-- Path to the folders of different hydraulic scenarios -->",
             "    <path_hydraulics> <!-- By order of selected property: Xcontact / leakiness / kw / Kpd / kAQP -->",
             paste0("        <Output path='", output_path, "'/>"),
             "    </path_hydraulics>",
             "    <!-- Cell wall hydraulic conductivity",
             "    Review in Transport in plants II: Part B Tissues and Organs, A. Lauchli: 1. Apoplasmic transport in tissues",
             "    Units: cm^2/hPa/d",
             "    6.6E-03: Soybean hypocotyl Steudle and Boyer (1985)",
             "    2.4E-04: Zhu and Steudle (1991)",
             "    1.2E-05: Nitella cell Tyree (1968)",
             "    6.0E-06: Nitella cell walls Zimmermann and Steudle (1975)",
             "    1.3E-7:  Cellulose wall Briggs (1967) for thickness of 1.5 micron",
             "    1.8E-9:  Maize root cell wall Tyree (1973) for thickness of 1.5 micron -->",
             "    <kwrange>",
             paste0("        <kw value='", Kw, "' /> <!-- hydraulic conductivity of standard walls -->"),
             "    </kwrange>",
             "    <!-- Apoplastic barrier conductivity",
             "    Units: cm^2/hPa/d",
             "    4.3E-07: 1 nanometer pores, Poiseuille law",
             "    1.0E-16: Virtually impermeable",
             "    Casp concerns the hydraulic conductivity of lignin in radial walls of 'gatekeeper cells' (endodermis/exodermis), while 'Sub' concerns the hydraulic conductivity of suberized secondary walls-->",
             "    <kw_barrier_range> <!-- hydraulic conductivity of suberised and lignified walls -->")
  
  # for(bi in seq(1, length(Barriers))){
  #   if(UniXwalls == T){
  #     if(Barriers[bi] == 1){ # If only CS
  #       lines <- c(lines, paste0("        <kw_barrier value='", kw,"' />"))
  #     }
  #     else{
  #       lines <- c(lines, paste0("        <kw_barrier value='", kwi,"' />"))
  #     }  
  #   }
  #   else{
  #     if(kwi == 1){
  #       lines <- c(lines, paste0("        <kw_barrier Casp='", kw,"' Sub='1.0E-16' />")) # CHANGER Sub ici !
  #     }
  #     else{
  #       lines <- c(lines, paste0("        <kw_barrier Casp='", kwi,"' Sub='1.0E-16'/>")) # CHANGER Sub ICI !
  #     }  
  #   }
  #   
  #   
  # }
  
  # for(bi in seq(1, length(Barriers))){
  #   if(UniXwalls){
  #     # If we only set xylem walls -> kw_CS is set to normal kw value
  #     lines <- c(lines, paste0("        <kw_barrier value='", kw,"' />"))
  #   }else{
  #     lines <- c(lines, paste0("        <kw_barrier Casp='", Kw_CS,"' Sub='", Kw_sub, "'/>")) # CHANGER Sub ICI !
  #   }
  # }
  
  lines <- c(lines, paste0(
    "        <kw_barrier Casp='", Kw_CS, "' Sub_out='", Kw_sub_out, "' Sub_in='", Kw_sub_in, "' />"
  ))
  
  lines <- c(lines,
             "    </kw_barrier_range>",
             "    <!-- Cell membrane permeability, with separate contribution of the biphospholipid layer (km) and AQP (kAQP) ",
             "    Attention, kaqp==0 is used as flag for each air-filled intercellular space, and cancels any flow across the membrane (as if km was equal to 0 too)",
             "    For knocked out AQPs, please use a very small value, such as 1.0E-16",
             "    Units: cm/hPa/d",
             "    6.8E-5 to 1.9E-4: Lp Steudle and Jeschke (1983)",
             "    4.3E-4: kAQP Elhert et al. (2009), difference in maize cell hydraulic conductivity between control and acid-treated cells (AQP closed)",
             "    3.0E-5: km after removal of kAQP and kpl from Elhert et al. (2009) and Bret-Hart and Silk (1994)-->",
             "    <km value='3.0E-5' />",
             "    <kAQPrange>",
             paste0("        <kAQP value='", kAQP, "'  epi_factor='1.0' exo_factor='1.0' cortex_factor='1.0' endo_factor='1.0' stele_factor='1.0' />"),
             "    </kAQPrange>",
             "    <ratio_cortex value='1.0'/> <!-- (-) ratio of cell membranes hydraulic conductivity between inner and outer cortex -->",
             "    <!-- Individual plasmodesma conductance and plasmodesmatal frequency",
             "    Kpl Units: cm^3/hPa/d/plasmodesmata",
             "    9.1E-13: Very low Kpl (extrapolation)",
             "    5.3E-12: Geometrical average from Bret-Hart and Silk (1994)",
             "    3.1E-11: Estimation from Ginsburg & Ginzburg (1970)",
             "    Frequency*height is the conserved quantity when height changes (default height in measured tissue assumed to be 100 microns)",
             "    Measurements reported below come from elongated maize cell (~100 microns)",
             "    Fpl by height Units: plasodesmata*cm/cm^2",
             "    3.5E5: Default average (0.35/micron^2 * 100 microns)",
             "    1.3E5: height * Plasmodesmal frequency between epidermis and cortex from Zhu et al. (1998)",
             "    5.1E5: height * Plasmodesmal frequency between cortex and cortex from Zhu et al. (1998)",
             "    2.6E5: height * Plasmodesmal frequency between cortex and endodermis from Zhu et al. (1998)",
             "    3.1E5: height * Plasmodesmal frequency between endodermis and endodermis from Zhu et al. (1998)",
             "    3.0E5: height * Plasmodesmal frequency between endodermis and pericycle from Zhu et al. (1998)",
             "    3.5E5: height * Plasmodesmal frequency between stele and stele (parenchyma) from Zhu et al. (1998) -->",
             "    <Kplrange>",
             paste0("        <Kpl value='", Kpl, "' cortex_factor='1.0' PCC_factor='1.0' PPP_factor='1.0' PST_factor='1.0' endo_in_factor='", Kpl_fac_in, "' endo_out_factor='", Kpl_fac_out, "' />"),
             "    </Kplrange>",
             "    <Freq_source value='1'/> <!-- 1: Hydraulics file (based on tissue type); 2: CellSet output file (wall per wall) -->",
             "    <Freq_factor value='1.0'/> <!-- Multiplies all PD frequencies -->",
             "    <Fplxheight value='3.5E5'/> <!-- Default height * plasmodesmatal frequency -->",
             "    <Fplxheight_epi_exo value='3.5E5'/> ",
             "    <Fplxheight_outer_cortex value='1.3E5'/>",
             "    <Fplxheight_cortex_cortex value='5.1E5'/>",
             "    <Fplxheight_cortex_endo value='2.6E5'/>",
             "    <Fplxheight_endo_endo value='3.1E5'/>",
             "    <Fplxheight_endo_peri value='3.0E5'/>",
             "    <Fplxheight_peri_peri value='3.5E5'/>",
             "    <Fplxheight_peri_stele value='3.5E5'/>",
             "    <Fplxheight_stele_stele value='3.5E5'/>",
             "    <Fplxheight_stele_comp value='3.5E5'/>",
             "    <Fplxheight_peri_comp value='3.5E5'/>",
             "    <Fplxheight_comp_comp value='3.5E5'/>",
             "    <Fplxheight_comp_sieve value='9.0E5'/>",
             "    <Fplxheight_peri_sieve value='3.5E5'/>",
             "    <Fplxheight_stele_sieve value='3.5E5'/>",
             "    <Defheight value='100'/> <!-- Default height (microns) -->",
             "    <!-- Axial conductances (cm^3/hPa/d) -->",
             "    <Kax_source value='1'/> <!-- 1: Poiseuille law (based on cross-section area); 2: Prescribed here below (for all sieve tubes, and vessel per vessel) -->",
             "    <K_sieve_range> <!-- Xylem vessel hydraulic conductance (cm^3/hPa/d) -->",
             "        <K_sieve id='38' value='4.0E-5'/> <!-- raw estimation -->",
             "        <K_sieve id='40' value='4.0E-5'/>",
             "        <K_sieve id='50' value='4.0E-5'/>",
             "        <K_sieve id='55' value='4.0E-5'/>",
             "    </K_sieve_range>",
             "    <K_xyl_range> <!-- Xylem vessel hydraulic conductance (cm^3/hPa/d) -->",
             "        <K_xyl id='24' value='1.1E-03'/> <!-- estimation based on transverse surface area -->",
             "        <K_xyl id='29' value='6.7E-03'/>",
             "        <K_xyl id='43' value='4.5E-03'/>",
             "        <K_xyl id='51' value='8.0E-03'/>",
             "        <K_xyl id='67' value='7.6E-04'/>",
             "    </K_xyl_range>",
             "    <!-- Xcontact (microns) is the X threshold coordinate of contact between soil and root (lower X not in contact with soil)  -->",
             "    <!-- Note that contact is also dealt with in the *Geometry.xml input file where the id's of cells in contact with the soil are listed  -->",
             "    <Xcontactrange>")
  
  for(i in seq(1, length(Contact))){
    lines <- c(lines, paste0("        <Xcontact value='", Contact[i], "'/> <!-- Arabido1: 120 microns includes 1 cell, 109 microns includes 5 cells to the right -->'"))
  }
  
  lines <- c(lines, "    </Xcontactrange>",
             "</param>")
  
  writeLines(lines, paste0(input_path, "Hydraulics_optim.xml"))
}

# ============================================================================ #

write_Hormones_xml <- function(input_path, 
                               Hormones_inputs,
                               cid, 
                               iLayer, 
                               lengths, 
                               heights,
                               washin_cc = c(1, 1, 0, 0), washout_cc = c(0, 0, 0, 0),
                               Ini_cond = 1){
  
  # Diff_PD = 6.2, Diff_PW = 8.62, Diff_X = 10.5, # Diffusivity PD, Primary Walls & Xylem
  Diff_PD   = Hormones_inputs$Diff_PD
  Diff_PW   = Hormones_inputs$Diff_PD
  Diff_X    = Hormones_inputs$Diff_X
  transient = Hormones_inputs$transient
  
  lines <- c(
    '<?xml version="1.0" encoding="utf-8"?>',
    ' ',
    '<param>',
    '    <!-- Hormone convection-diffusion-degradation (activated by setting the value of Sym_contagion and/or Apo_contagion to 2 in *_General.xml)',
    '         Hormone movement is calculated at steady-state within the symplast or apoplast.',
    '         Degradation',
    '         Units: mol degraded / mol-day',
    '         The diffusivity & degradation constant allows avoiding the infinite accumulation of solutes in locations that do not export the hormone.',
    '         24: average lifetime of an hour (Normanly 1997, Auxin metabolism)',
    '         Diffusivity',
    '         Units: cm^2 / day',
    '         0.14: value for carboxyfluorescein in symplast (Deinum et al., 2019) 162 micron^2 / s',
    '         0.3: value for auxin, scaled by using law of diffusivity inversely proportional to molecule size (Deinum et al., 2019)',
    '         0.00069: value for water in phospholipid bilayers (Khakimov et al. 2008) 8e-13 m^2/s',
    '         0.035: value for mannitol in cell walls (Knipfer et al., 2007)',
    '         1.2: value for K (Foster & Miklavcic, 2017)',
    '         2.0: Value for D2O in water (Spyrou et al., 2016)',
    '         Partition coefficient',
    '         Units: -',
    '         4.2e-5: Water in phospholipid bilayers (Khakimov et al. 2008)',
    '         Membrane thickness',
    '         Units: micron',
    '         3e-3: Phospholipid bilayer thickness (Reus et al., 2008)',
    '         The source cells are defined in *_Geometry*.xml, and the concentration in source cells is set to the constant 1 (relative concentration) -->',
    '    <Hormone_movement>',
    '        <Degradation_constant_H1 value="0.0"/>',
    '        <MB_H1 Diffusivity="0.00069" Partition_coef="4.2e-5" Thickness="3.0e-3"/>',
    paste0('        <Diffusivity_PD_H1 value="', Diff_PD, '"/>'), # *11e-9 ?
    paste0('        <Diffusivity_PW_H1 value="', Diff_PW, '"/>'), # *11e-9 ?
    paste0('        <Diffusivity_X_H1 value="',  Diff_X, '"/>'),  # *11e-9 ?
    '        <Dif_BC flag="1"/>',
    '        <H1_D2O flag="1"/>',
    '    </Hormone_movement>',
    '    ',
    '    <!-- Temporal information for transient simulations of tracer movement ',
    '         **Set Transient_sim to 1 in order to simulate tracer transient flow**',
    '         Units: day',
    '         Time steps are not adaptive',
    '         Tracer outputs are recorded at the selected print times, and at time 0',
    '         The simulation ends at the last print time -->',
    '    <Time_information>',
    paste0('        <Transient_sim value="', transient, '"/>'),
    '        <Switch_duration value="1.16e-9"/> <!-- Transition time during which the boundary condition progressively changes from ini to transi values -->',
    '        <Switch_fractionning value="4.0e1"/> <!-- Fractionning of the transition time into small time steps of duration "Switch_duration/Switch_fractionning" -->',
    '        <Time_step_transi value="1.0e-9"/> <!-- 3.48e-6 -->',
    '        <Time_step_adapt value="1" cc_change_low="0.1" cc_change_high="0.2" time_step_factor="1.4" time_step_min="1.0e-9" time_step_max="2.0e-4"/> <!-- value 0: Constant time step; value 1: Adaptive time step -->',
    '        <Observation_frequency value="1"/> <!-- How many time steps before saving tracer concentration at observation points -->',
    '        <Print_times>',
    '            <Print_time value="1.5e-2"/>',
    #'                <Print_time value="5.0e-5"/>',
    #'                <Print_time value="5.0e-4"/>',
    #'                <Print_time value="5.0e-3"/>',
    #'                <Print_time value="1.0e-2"/>',
    #'                <Print_time value="1.5e-2"/>',
    '        </Print_times>',
    '    </Time_information>',
    '    ',
    '    <!-- Observation points',
    '         Records tracer concentration at all times',
    '         id corresponds to the cell id number, starting at 0 (see cell_map from CellSet or Granar)',
    '         layer corresponds to the cross-section number, distal being 0 -->',
    '    <Observation_range> <!-- xylem parenchyma: 24, 29, 43, 51, 67   companion cells: 32, 42, 52, 56 -->')
  
  for(i in seq(1, length(cid))){
    for(j in seq(1, length(iLayer))){
      lines <- c(lines, paste0('        <Obs id="', cid[i], '"  layer="', iLayer[j], '"  precision="high" />'))
    }
  }
  
  lines <- c(lines, 
             '    </Observation_range>',
             '    ',
             '    <!-- Numerical method to solve transient tracer movement ',
             '         1: Explicit',
             '         2: Implicit 1',
             '         3: Implicit 2',
             '         -->',
             '    <Numerical_method>',
             '        <Method value="3" />',
             '    </Numerical_method>',
             '    ',
             '    <!-- Type of initial conditions for transient simulations of tracer movement ',
             '         1: From simulated steady state with prescrbed "concentration_ini" and "flowrate_ini", see Sym_contagion and/or Apo_contagion inputs',
             '         2: From prescribed concentrations in *3D_cc*.xml file',
             '         -->',
             '    <Initial_conditions>')
  
  if(Ini_cond==1){
    lines <- c(lines, paste0('        <Type value="1" path="none" />'))
  }
  else{
    lines <- c(lines, paste0('        <Type value="2" path="', Ini_cond, '"/>'))
  }
  
  lines <- c(lines,   
             '    </Initial_conditions>',
             '    ',
             '    <!-- Hormone active transport (activated by setting the value of Sym_contagion and Apo_contagion to 2 in *_General.xml)',
             '         Active transport is calculated assuming that hormone concentration is low enough to be proportional to concentration (linear part of the Michaelis-Menten curve)',
             '         Efflux carriers send the hormone from the symplast to the apoplast, while influx cariiers do the opposite',
             '         Active transport constants (=Vmax/KM)',
             '         Units: liter / day-micron^2',
             '         3.9E-16: ABA efflux carrier ABCG25 (Kuromori et al. 2010), assuming 55 kDa carrier, and 100 carriers per micron^2',
             '         5.8E-14: ABA influx carrier NFP4, previously AIT1 (Kanno et al., 2012), assuming same carrier density as in transformed yeast, and yeast spherical cell diameter of 6 microns',
             '         Direction: +1 is influx, -1 is efflux',
             '         Tissue: 1 		= exodermis',
             '                 2 		= epidermis',
             '                 3 		= endodermis',
             '                 4 		= cortex',
             '                 5 		= stele (excluding the pericycle)',
             '                 11, 23 = phloem',
             '                 12, 26 = companion cell',
             '                 13, 19 = protoxylem',
             '                 13, 20 = metaxylem',
             '                 16, 21 = pericycle (21 being the xylem pole pericycle)',
             '         The source cells are defined below, and the concentration in source cells is set to the constant 1 (relative concentration) -->',
             '    <Hormone_active_transport>',
             '        <carrier_range>',
             '            <carrier tissue="-1"  constant="7.9E-11" direction="-1"/> <!-- xylem parenchyma: 24, 29, 43, 51, 67   companion cells: 32, 42, 52, 56 -->',
             '        </carrier_range>',
             '    </Hormone_active_transport>',
             '    ',
             '    <!-- Symplastic contagion inputs',
             '         If flow_BC is set to 0: cc_dist is the prescribed static relative concentration in the cell (cc_prox is not used). Units: relative concentration',
             '         If flow_BC is set to 1: The concentration will be multiplied by the incoming water flow rate in order to get the concentration flux arriving into the cell. Units: relative concentration to be multiplied by the associated convective flux in cm^3/d',
             '         In case of solute flow BC, concentrations in the distal and proximal ends can be prescribed, though only the concentration upstream of water flow will be considered',
             '         In case of transient simulation, steady-state BC generate initial conditions from the steady-state solution (except if cc prescribed in input file). Transient BC define the new boundary condition at time 0.',
             '         The current code does not allow to have both pressure- and flux-type boundary conditions at the same node',
             '         Set Transient_sim to "nan" in Time_information in order not to have a transient simulation',
             '         -->',
             '    <Sym_Contagion>',
             '        <source_range>',
             '            <Steady-state>',
             '                <source id="-1" layer_range_min="0" layer_range_max="0" cc_dist="1.000000" cc_prox="0.000100" flow_BC="1"/> <!-- xylem parenchyma: 24, 29, 43, 51, 67   companion cells: 32, 42, 52, 56 -->',
             '            </Steady-state>',
             '            <Transient>',
             '                <source id="-1" layer_range_min="0" layer_range_max="0" cc_dist="0.000100" cc_prox="0.000100" flow_BC="1"/> <!-- xylem parenchyma: 24, 29, 43, 51, 67   companion cells: 32, 42, 52, 56 -->',
             '            </Transient>',
             '        </source_range>',
             '        <target_range>',
             '            <target id="-1"/> <!-- Low water potential side cortex Hydrotropism: 78, 79 -->',
             '        </target_range>',
             '        <immune_range>',
             '            <immune id="-1"/>',
             '        </immune_range>',
             '    </Sym_Contagion>',
             '    ',
             '    <!-- Apoplastic contagion inputs ',
             '         Cell walls can either be selected through the id number of their cell (see Cell_map*.png) or through their tissue identifier (see above)',
             '         For Dirichlet concentration BC at root surface, select tissue "2", prescribe the concentration in cc_soil, and set flow_BC to 0',
             '         For Neumann solute flux BC at root surface, select tissue "2", prescribe the concentration in cc_soil, and set flow_BC to 1		 -->',
             '    <Apo_Contagion>',
             '        <source_range>',
             '            <Steady-state>')
  
  temp <- 0
  
  for(i in seq(1, length(washout_cc))){
    lines <- c(lines, 
               paste0('                <source tissue="2" id="-1" layer_range_min="', temp, '" layer_range_max="', temp+(lengths[i]/heights[i])-1, '" cc_dist="nan" cc_prox="nan" cc_soil="', washin_cc[i], '" flow_BC="1"/> <!-- xylem parenchyma: 24, 29, 43, 51, 67 -->')
    )  
    temp <- temp+(lengths[i]/heights[i])
  }
  
  # ['                <source tissue="2" id="-1" layer_range_min="' num2str(sum(lengths(1:2)./heights(1:2))) '" layer_range_max="' num2str(sum(lengths./heights)-1) '" cc_dist="nan" cc_prox="nan" cc_soil="0.000000" flow_BC="1"/>'],
  lines <- c(lines, 
             '            </Steady-state>',
             '            <Transient>')
  
  temp <- 0
  
  for(i in seq(1, length(washout_cc))){
    lines <- c(lines, 
               paste0('                <source tissue="2" id="-1" layer_range_min="', temp, '" layer_range_max="', temp+(lengths[i]/heights[i])-1, '" cc_dist="nan" cc_prox="nan" cc_soil="', washout_cc[i], '" flow_BC="1"/> <!-- xylem parenchyma: 24, 29, 43, 51, 67 -->'))
    temp <- temp+(lengths[i]/heights[i])
  }
  
  #['                <source tissue="2" id="-1" layer_range_min="' num2str(sum(lengths(1:2)./heights(1:2))) '" layer_range_max="' num2str(sum(lengths./heights)-1) '" cc_dist="nan" cc_prox="nan" cc_soil="0.000000" flow_BC="1"/>'],
  
  lines <- c(lines, 
             '            </Transient>',
             '        </source_range>',
             '        <target_range>',
             '            <target id="-1"/> <!-- Low water potential side cortex Hydrotropism: 78, 79 -->',
             '        </target_range>',
             '        <immune_range>',
             '            <immune id="-1"/>',
             '        </immune_range>',
             '    </Apo_Contagion>',
             '    ',
             '    <!-- Contact (cell id) is the list of cells in contact with the water meniscus at the interface with the soil  -->',
             '    <Contactrange>',
             '        <Contact id="-1"/>',
             '    </Contactrange>',
             '</param>')
  
  # Write the vector of lines to a file
  writeLines(lines, paste0(input_path, "Hormones_optim.xml"))
}

