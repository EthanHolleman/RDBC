library(ggplot2)
library(bio3d)
library(reshape)
library(RColorBrewer)

score_file_path = '/home/ethan/Documents/agg_results.tsv'
pdb_dir = '/home/ethan/Documents/gino/shc_conformation_comparison/structures/top_15_repacked'

read_score_file <- function(score_file_path){
  read.delim(score_file_path, header=T)
}

total_score_vs_interface_delta_by_protein <- function(df){
  ggplot(df, aes(x=total_score, y=interface_delta_X, color=ligand_name)) + 
    labs(x="Total Score", y="Interface Delta") + geom_point() + 
    facet_wrap(~target_protein) + theme_minimal()
}

save_plot <- function(plot, path){
  
  ggsave
  
}

read_pdb_dir_as_list <- function(pdb.dir){
  files <- list.files(pdb.dir, pattern = "\\.pdb$")
  pdb.list <- list()
  for (f in files){
    f.abs <- paste(pdb.dir, f, sep='/')
    pdb.list[[f]] <- read.pdb(f.abs)
  }
  return(pdb.list)
}

pairwise_rmsd <- function(protein.list){
  # compare (heatmap of RMSD of all a collection or protiens)
  # x and y are the same here need to make a data farme in
  # the format var 1, var 2, value
  # http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
  pairwise.rmsd <- list()
  c <- 1
  for (name_a in names(protein.list)){
    protein.a <- protein.list[[name_a]]
    for (name_b in names(protein.list)){
      protein.b <- protein.list[[name_b]]
      coords.a <- protein.a$xyz
      coords.b <- protein.b$xyz
      print(c)
      pairwise.rmsd[[c]] <- c(name_a, name_b, rmsd(coords.a, coords.b))
      c <- c + 1
    }
  }
  pairwise.rmsd.df <- do.call(rbind, pairwise.rmsd)
  # pairwise.rmsd.df <- pairwise.rmsd.df[lower.tri(pairwise.rmsd.df)] <- NA
  pairwise.rmsd.df <- as.data.frame(pairwise.rmsd.df)
  colnames(pairwise.rmsd.df) <- c('protein_a', 'protein_b', 'rmsd')
  return(pairwise.rmsd.df)
}

protein_rmsd_dendrogram <- function(protein.dir){
  
  pdb.list <- read_pdb_dir_as_list(protein.dir)
  rmsd.df <- pairwise_rmsd(pdb.list)
  rmsd.df.wide <- as.data.frame(pivot_wider(rmsd.df, names_from = protein_a, values_from = rmsd))
  rownames(rmsd.df.wide) <- rmsd.df.wide[, 1]
  rmsd.df.wide <- rmsd.df.wide[, -1]
  hclust(dist(rmsd.df.wide))
  
  
}

pairwise_rmsd_heatmap <- function(protein.dir){
  pdb.list <- read_pdb_dir_as_list(protein.dir)
  rmsd.df <- pairwise_rmsd(pdb.list)
  rmsd.df$rmsd <- as.numeric(rmsd.df$rmsd)
  ggplot(rmsd.df, aes(x=protein_a, y=protein_b, fill=rmsd)) + geom_tile()
}

# Want a heatmap of best scores arent chaning that much 

style_pairwise_heatmap <- function(heatmap){
  # http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
  heatmap + geom_text(aes(protein_a, protein_b, label = rmsd), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank()) + 
    scale_fill_distiller(palette = "Spectral") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    theme(axis.text=element_text(size=12),
           axis.title=element_text(size=14,face="bold"))
      #legend.justification = c(1, 0),
      # legend.position = c(0.6, 0.7),
      # legend.position = 'bottom',
      # legend.direction = "horizontal")
    # guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
    #                              title.position = "top", title.hjust = 0.5)) 
}

score_function <- function(ligand.ts, ligand.idx, min.ts, max.ts, min.idx, max.idx){
  
  ((2 * ((ligand.ts - min.ts) / (max.ts - min.ts)))-1) + ((2 * ((ligand.idx - min.idx) / (max.idx - min.idx)))-1)
  
}

score_ligand <- function(agg.df, ligand.name, protein.name){
  agg.df.ligand <- subset(agg.df, ligand_name==ligand.name & target_protein==protein.name)
  if (nrow(agg.df.ligand) > 0){
    min.ts <- min(agg.df.ligand$total_score)
    max.ts <- max(agg.df.ligand$total_score)
    min.idx <- min(agg.df.ligand$interface_delta_X)
    max.idx <- max(agg.df.ligand$interface_delta_X)
    
    scores <- c()
    for (i in 1:nrow(agg.df.ligand)){
      ligand.ts <- agg.df.ligand[i, ]$total_score
      ligand.idx <- agg.df.ligand[i, ]$interface_delta_X
      scores <- c(scores, score_function(ligand.ts, ligand.idx, min.ts, 
                                         max.ts, min.idx, max.idx))
    }
    print(length(scores))
    agg.df.ligand.scores <- cbind(agg.df.ligand, scores)
    colnames(agg.df.ligand.scores) <- c(colnames(agg.df.ligand), 'composite_score')
    return(agg.df.ligand.scores)
  }else{
    return(data.frame())
  }
  }

add_scores_to_agg_df <- function(agg.df){
  dfs <- list()
  i = 1
  for (ligand in unique(agg.df$ligand_name)){
    for (protein in unique(agg.df$target_protein)){
      print(paste(ligand, protein))
      scores.df <- score_ligand(agg.df, ligand, protein)
      if (nrow(scores.df) > 0){
        dfs[[i]] <- scores.df
      }
      print(i)
      i <- i + 1
    }
  }
  do.call(rbind, dfs)
}

get_best_poses_from_scored_agg_df <- function(agg.df.scored){
  
  
  
}



rmsd_ligands_from_pdb_list <- function(pdb_list){
  
  # rmsd.list <- list()
  # i <- 1
  # for 
  
}


rank_ligands_from_agg_df <- function(agg.df){
  
}

# HCLUSTERING 
# pivot_wider(rmsd.df, names_from = protein_a, values_from = rmsd)
# rownames(g) <- g[, 1]
# g <- g[, -1]





