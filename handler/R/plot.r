library(ggplot2)
library(bio3d)
library(reshape)

read_score_file <- function(score_file_path){
  read.delim(score_file_path, header=T)
}

pairwise_rmsd <- function(protein.list){
  # compare (heatmap of RMSD of all a collection or protiens)
  # x and y are the same here need to make a data farme in
  # the format var 1, var 2, value
  # http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
  pairwise.rmsd <- list()
  c <- 1
  for (i in 1:length(protein.list)){
    protein.a <- protein.list[[i]]
    for (j in 1:length(protein.list)){
      protein.b <- protein.list[[j]]
      pairwise.rmsd[[c]] <- c(name a, name b, rmsd(cords a, cords b))
    }
  }
}

style_pairwise_heatmap <- function(heatmap){
  # http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
  heatmap + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5)) 
  
  
}


