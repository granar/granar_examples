load_MECHA_out <- function(path = "MECHA/Projects/granar/out/Tomato/Root/Project_Test/results/"){
  ls_files <- list.files(path, pattern = ".txt")
  conds <- tibble()
  
  for(i in seq(1, length(ls_files))){
    
    if(ls_files[i] == "position.txt"){
      next
    }
    
    out.temp <- unlist(read.delim(paste0(path, ls_files[i])))
    kx.temp <- unlist(strsplit(out.temp[7], " "))[5]  # cm^4/hPa/d
    kr.temp <- unlist(strsplit(out.temp[8], " "))[4]  # cm/hPa/d
    peri.temp <- unlist(strsplit(out.temp[6], " "))[3]  # cm
    
    df.temp <- tibble("ID" = i,
                      "Kx" = kx.temp,
                      "kr" = kr.temp,
                      "perimeter" = peri.temp)
    
    conds <- rbind(conds, df.temp)
  }
  
  return(conds)
}

ApoSymp <- function(path = "MECHA/Projects/granar/out/Root/Project_Test/baseline/Macro_prop_1,0.txt"){
  
  coef_width_symplast=4/5
  mpercm=0.01
  dpersec=1/3600/24
  
  mecha_out_0 <- readLines(path)
  pos_dist <- which(mecha_out_0 =="Radial distance from stele centre (microns): ")
  pos_STF_up <- which(mecha_out_0 =="Standard Transmembrane uptake Fractions (%): ")
  pos_STF_rel <- which(mecha_out_0 =="Standard Transmembrane release Fractions (%): ")
  pos_end <- which(mecha_out_0 == "Scenario 1 ")
  pos_UP <- which(mecha_out_0 =="Stele, cortex, and epidermis uptake distribution cm^3/d: ")
  pos_UN <- which(mecha_out_0 =="Stele, cortex, and epidermis release distribution cm^3/d: ")
  pos_QXyl <- which(mecha_out_0 =="Xylem uptake distribution cm^3/d: ")
  pos_QPhlo <- which(mecha_out_0 =="Phloem uptake distribution cm^3/d: ")
  pos_disc <- which(mecha_out_0 =="Number of radial discretization boxes: ")
  
  q_tot <- parse_number( mecha_out_0[grepl("q_tot", mecha_out_0)]) # cm2 d-1
  c_height <- parse_number(unlist(str_split( mecha_out_0[grepl("height:", mecha_out_0)], " ")))
  c_height <- c_height[!is.na(c_height)]
  
  Q_tot = (q_tot*c_height)
  
  
  disc <- as.numeric(unlist(str_split(mecha_out_0[pos_disc+1], " ")))
  disco <- tibble(nlayer = disc[!is.na(disc)][-1], type = c("stele", "pericycle", "endodermis", "cortex", "exodermis", "epidermis"))
  
  
  dist <- parse_number(mecha_out_0[seq(pos_dist+1,pos_STF_up-2,1) ])
  STF_up <- parse_number(mecha_out_0[seq(pos_STF_up+1,pos_STF_rel-2,1) ])
  STF_down <- parse_number(mecha_out_0[seq(pos_STF_rel+1,pos_end-3,1) ])
  
  UP1 <- parse_number(mecha_out_0[seq(pos_UP[1]+1,pos_UN[1]-2,1) ])
  UN1 <- parse_number(mecha_out_0[seq(pos_UN[1]+1,pos_QXyl[1]-2,1) ])
  Qxyl1 <- parse_number(mecha_out_0[seq(pos_QXyl[1]+1,pos_QPhlo[1]-2,1) ])
  
  flux <- tibble(dist = dist, STF_up, STF_down, UP1,UN1, Qxyl1)
  
  Type <- NULL
  for (k in 1:nrow(disco) ){
    tmp <- rep(disco$type[k], disco$nlayer[k])
    Type <- c(Type, tmp)
  }
  flux$type <- Type
  
  flux$passage <- 0
  flux$passage[flux$type == "endodermis"] <- c(0,1,1,0)
  flux$intern <- 0
  flux$intern[flux$type == "endodermis"] <- c("in","in","out","out")
  
  flux_passage <- flux
  flux <- flux[flux$passage == 0 , ]
  
  flux$UP1[flux$intern == "in"] <- sum(flux_passage$UP1[flux_passage$intern == "in"])
  flux$UN1[flux$intern == "in"] <- sum(flux_passage$UN1[flux_passage$intern == "in"])
  flux$UP1[flux$intern == "out"] <- sum(flux_passage$UP1[flux_passage$intern == "out"])
  flux$UN1[flux$intern == "out"] <- sum(flux_passage$UN1[flux_passage$intern == "out"])
  
  endo <- flux$dist[flux$type == "endodermis"][1]
  flux <- flux %>%# filter(!duplicated(dist))%>%
    mutate(r_coord_layer = dist - endo)
  
  Flux <- NULL
  for(h in 2:nrow(flux)){
    tmp_m <- flux[c(h-1,h),]
    if((unique(tmp_m$type)[1] %in% c("exodermis", "endodermis") | unique(tmp_m$type)[2] %in% c("exodermis", "endodermis")) & length(unique(tmp_m$type)) == 1){
      tmp <- tmp_m%>%mutate(label = "gate")
    }else{
      
      mid_wall <- mean(tmp_m$r_coord_layer)
      memb = (mid_wall-tmp_m$r_coord_layer[1])*coef_width_symplast
      tmp_m$label <- "mid_cell"
      
      if(tmp_m$Qxyl1[2] != 0){
        Addon <- tmp_m[2,]%>%mutate(r_coord_layer = mid_wall,
                                    label = "mid_wall",
                                    Qxyl1 = 0)
      }else{
        Addon <- NULL
      }
      
      tmp <- rbind(tmp_m[1,], # mid cell
                   tmp_m[1,]%>%mutate(r_coord_layer = r_coord_layer+memb,
                                      label = "membrane_in"), # membrane
                   tmp_m[1,]%>%mutate(r_coord_layer = r_coord_layer+memb,
                                      label = "membrane_out"), # wall
                   Addon,
                   tmp_m[2,]%>%mutate(r_coord_layer = mid_wall,
                                      label = "mid_wall"), # mid wall
                   
                   tmp_m[2,]%>%mutate(r_coord_layer = r_coord_layer-memb,
                                      label = "membrane_out"), # wall
                   tmp_m[2,]%>%mutate(r_coord_layer = r_coord_layer-memb,
                                      label = "membrane_in"), # membrane
                   tmp_m[2,]
      )
    }
    Flux <- rbind(Flux, tmp)
    
  }
  
  Flux <- Flux %>%mutate(r_coord_layer = -r_coord_layer)%>%arrange (r_coord_layer)
  
  Flux$type[Flux$label == "mid_wall"] <- "cell_wall"
  Flux$id <- 1:nrow(Flux)
  
  Flux$Sympl <- NA
  Flux$Apo <- NA
  
  Flux$Sympl[1:2] <- Flux$UP1[1:2]
  Flux$Sympl[3] <- Flux$Sympl[2]+Flux$UN1[3]
  Flux <- rbind(Flux[1,]%>%mutate(Sympl = 0,
                                  Apo = 0,
                                  type = "soil",
                                  r_coord_layer = r_coord_layer-diff(Flux$r_coord_layer[c(1,3)])/5),
                Flux[1,]%>%mutate(Sympl = 0,
                                  Apo = Q_tot,
                                  type = "soil",
                                  r_coord_layer = r_coord_layer-diff(Flux$r_coord_layer[c(1,3)])/5),
                Flux[1,]%>%mutate(Sympl = 0, Apo = Q_tot),
                Flux)
  
  com <- Flux[grepl("membrane", Flux$label),]
  
  for(i in 3:nrow(com)){
    if(i %% 2 == 0){
      if(com$type[i] %in% c("exodermis", "endodemris") ){
        com$Sympl[i] <- com$Sympl[i-1]+com$UN1[i]+com$UP1[i]
      }else{
        if (length(com$label[i][grepl("out", com$label[i])])>0){
          com$Sympl[i] <- com$Sympl[i-1]+com$UP1[i]
        }else{
          com$Sympl[i] <- com$Sympl[i-1]+com$UN1[i]
        }
      }
    }else{
      com$Sympl[i] <- com$Sympl[i-1]
    }
    Flux$Sympl[Flux$id == com$id[i]] <- com$Sympl[i]
  }
  
  
  cxyl <- Flux[grepl("mid_wall", Flux$label),]
  cxyl$Apo[1] <- Q_tot
  
  for(i in 2:nrow(cxyl)){
    cxyl$Apo[i] <- cxyl$Apo[i-1]+cxyl$Qxyl1[i]
    Flux$Apo[Flux$id == cxyl$id[i]] <- cxyl$Apo[i]
  }
  
  return(Flux)
}

plot_water_flux <- function(Flux, apobar = 1){
  Q_tot = Flux$Apo[2]
  r_0 <- Flux$r_coord_layer[1]
  
  colo <- c("Apoplastic compartment" = "burlywood3", "Symplastic compartment"= "dodgerblue2")
  
  pl <- Flux%>%
    filter(!is.na(Sympl) )%>%
    ggplot()+
    geom_polygon(aes(r_coord_layer, Apo, fill = "Apoplastic compartment"), data = Flux%>%filter(!is.na(Apo) ))+
    geom_polygon(aes(r_coord_layer, Sympl, fill = "Symplastic compartment") )+
    geom_line(aes(r_coord_layer, Sympl))+
    geom_vline(aes(xintercept = r_coord_layer), alpha = 0.2, linetype = 2, data = Flux%>%filter(label == "mid_wall"))+
    geom_vline(aes(xintercept = r_coord_layer), alpha = 0.5, linetype = 2, data = Flux%>%filter(type == "soil"))+
    geom_hline(yintercept = Q_tot)+
    geom_hline(yintercept = 0)+
    xlim(r_0-2, 90)+
    theme_classic()+
    ylab("Partition of radial water flow rates \n between compartments [cm3/d]")+
    xlab("Distance from the endodermis layer [um]")+
    labs(fill = "Compartment")+
    scale_fill_manual(values = colo)
  
  if(apobar == 1){
    pl <- pl + geom_vline(aes(xintercept = 0), alpha = 0.6, size = 2, color = "red")
  }
  if(apobar == 2){
    pl <- pl + geom_vline(aes(xintercept = 0), alpha = 0.6, size = 2, color = "red")+
      geom_vline(aes(xintercept = r_coord_layer), alpha = 0.6, size = 2, color = "red",
                 data = Flux%>%filter(id %in% c(48,57)))
  }
  if(apobar == 3){
    pl <- pl + geom_vline(aes(xintercept = 0), alpha = 0.6, size = 2, color = "red")+
      geom_vline(aes(xintercept = r_coord_layer), alpha = 0.6, size = 2, color = "red",
                 data = Flux%>%filter(id %in% c(48,57)))+
      geom_vline(aes(xintercept = r_coord_layer), alpha = 0.6, size = 2, color = "red",
                 data = Flux%>%filter(label == "mid_cell",
                                      type == "exodermis"))
    
  }
  
  
  print(pl)
}
