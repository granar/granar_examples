
library(gganimate)


setwd("./granar/")
document()
setwd("..")
install("granar")

params <- read_param_xml("input_granar.xml")

nodes <- NULL
for(i in c(2:6)){
  params$value[params$name == "xylem" & type == "n_files"] <- i
  sim <- create_anatomy(parameters = params, verbatim=T)
  temp <- sim$nodes
  temp$rand <- i
  nodes <- rbind(nodes, temp)
}

pl <- nodes %>% 
  ggplot(aes(x, y, group=id_cell, fill=type, frame = rand)) + 
  geom_polygon(colour="white") + 
  #geom_text(data=sim$cells, aes(x, y, label=id_cell-1), size=3) + 
  theme_classic() + 
  coord_fixed() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) 

gganimate(pl)
