# dendrogram creation and plotting 
library(parallel)
library(bio3d)
library(factoextra)
library(argparse)
#library(plotly)
library(ggplot2)

cl <- parallel::makeCluster(4)

get_arg_parser <- function(){
  parser <- ArgumentParser(description='Analyize Rosetta ligand RMSD based heirarhical clustering.')
  parser$add_argument('-r', '--rds', help='Path to RDS file containing a RMSD distance matrix for all ligand poses.')
  parser$add_argument('-s', '--score_file', help='Path to score file containing results of all ligand docking simulations contained by the file specified with the -r argument. This should be formated as a tsv file from the original Rosetta output.')    
  parser$add_argument('-p', '--pdb_dir', default="", help='Path to directory containing pdb results of ligand docking simulations.')
  parser$add_argument('-m', '--max_poses', default=500, help='Max number of poses to evaluate. Default 500. Increasing this greatly increases computation time.')
  parser$add_argument('-o', '--output_dir', default='', help='Folder to save plots to.')
  parser$add_argument('-k', '--max_clusters', default=25, help='Max number of clusters to use.')
  return(parser)
  # args <- parser$parse_args()
  # args$max_clusters <- as.numeric(args$max_clusters)
  # args$max_poses <- as.numeric(args$max_poses)
  # return(args)
}

process_parsed_args <- function(args){
  args$max_clusters <- as.numeric(args$max_clusters)
  args$max_poses <- as.numeric(args$max_poses)
  return(args)

}


read_rds <- function(rds.path, max_poses=500){
    rmsd.df <- readRDS(rds.path)
    if (nrow(rmsd.df) > max_poses){
        rmsd.df <- rmsd.df[1:max_poses, 1:max_poses]
    }
    rownames(rmsd.df) <- names(rmsd.df) # should be a square distance matrix
    return(rmsd.df)
}

read_score_file <- function(score_file.path){
    return(read.csv(score_file.path))
}

hcut_rmsd_df <- function(rmsd.df, k=25){
    gap <- fviz_nbclust(rmsd.df, hcut, 'gap_stat', k.max=k)
    message(1)
    n_cuts <- match(max(gap$data$gap), gap$data$gap)
    message(2)
    hcut <- hcut(rmsd.df, k=n_cuts, hc_method='complete')
    message(3)
    return(hcut)
}

# Scorefile description field will not have the .pdb extension but this
# extension should be present as the row / col names of the distance
# matrix read using the read_rds file. In order to match these two data
# frames up we add the the description + /pdb extension to a new field
# called "description_pdb"
add_filepath_to_scorefile_description <- function(score.df){
    descrip_with_ext <- lapply(score.df$description, function(s){paste(s, '.pdb', sep='')})
    descrip_with_ext <- unlist(descrip_with_ext)
    score.df$description_pdb <- descrip_with_ext
    return(score.df)
}


select_scorefile_rows_present_in_dist_matrix <- function(score.df, rmsd.df){
    score.df.subset <- subset(score.df, description_pdb %in% rownames(rmsd.df))
    return(score.df.subset)
}


add_cluster_field_to_score_df <- function(score.df, hcut){
    score.df <- score.df[order(names(hcut$cluster)), ]
    score.df$cluster <- hcut$cluster
    return(score.df)
}


remove_uncontacted_ligands <- function(score.df){
  s <-  score.df
    score.df.contact <- subset(score.df, ligand_is_touching_X==1)
    return(score.df.contact)
}


interface_delta_by_cluster_boxplot <- function(score.df){
    ggplot(score.df, aes(x=as.factor(cluster), y=interface_delta_X)) + 
    geom_boxplot(fill='#8b21a3') + 
    theme_minimal() + labs(x='Cluster', y='Interface Delta')
}

total_score_by_cluster_boxplot <- function(score.df){
    ggplot(score.df, aes(x=as.factor(cluster), y=total_score)) + 
    geom_boxplot(fill='#8b21a3') + 
    theme_minimal() + labs(x='Cluster', y='Total Score')
}

# Calculate the average coordinates of a ligand chain (X) from a pdb object
average_coord_of_chain <- function(pdb.object){
    chain_coords <- subset(pdb.object$atom, chain=='X', select=c('x', 'y', 'z'))
    mean_coords <- c(mean(chain_coords$x), mean(chain_coords$y), mean(chain_coords$z))
    return(mean_coords)
}

# Read pdb files from a directory (should be same directory that the score
# file is from) and calculate the average coordinates of the ligand (specified
# by chain id 'X') and add to the score file dataframe.
add_average_chain_coords_to_score_df <- function(pdb.dir, score.df){
    files.all <- list.files(pdb.dir, patter='\\.pdb$', full.names = T)
    files.select <- list()
    for (f in files.all){  # get files actually in the 
    if (basename(f) %in% score.df$description_pdb){
      files.select[[basename(f)]] <- f
    }
  }
  pdb.objects <- parLapply(cl, files.select, read.pdb)
  mean.coords <- parLapply(cl, pdb.objects, average_coord_of_chain)
  mean.coords <- t(mean.coords)
  mean.coords.df <- do.call(rbind, mean.coords)
  rownames(mean.coords.df) <- names(mean.coords)
  colnames(mean.coords.df) <- c('mean_ligand_x', 'mean_ligand_y', 'mean_ligand_z')
  cbind(score.df, mean.coords.df)
}


distance_to_binding_pocket <- function(score.df, binding_coords=c(-2.32, -13.109,-5.0907)){
  distance_to_pocket <- c()
  for (i in 1:nrow(score.df)){
    d <- dist.xyz(
      c(score.df[i, ]$mean_ligand_x, score.df[i, ]$mean_ligand_y, score.df[i, ]$mean_ligand_z), 
      binding_coords
      )
    distance_to_pocket <- c(distance_to_pocket, d)
  }
  score.df$distance_to_pocket <- distance_to_pocket
  return(score.df)
}


average_ligand_coords_by_cluster_plot <- function(score.df){
    fig <- plot_ly(mean.coords, x = ~mean_ligand_x, y = ~mean_ligand_y, z = ~mean_ligand_z, color = ~as.factor(cluster))
    fig %>%
    add_trace(
        marker = list(
        size = 8,
        line = list(
            color = 'rgb(0, 0, 0)',
            width = 2
        )
        ),
        showlegend = F
    )
    return(fig)
}

get_top_n_percent_by_metric <- function(score.df, metric, percent=0.25){
  score.df.sorted <- score.df[order(score.df[, metric]), ]
  top_n_rows <- round(nrow(score.df.sorted) * percent)
  return(score.df.sorted[1:top_n_rows, ])
}

average_ligand_distance_to_point_vs_delta_x <- function(score.df, x_lab=''){
    ggplot(score.df, aes(distance_to_pocket, y=interface_delta_X,
                         text = paste(
                           "Distance to PTB: ", distance_to_pocket, "\n",
                           # "Interface Delta: ", interface_delta_X, "\n",
                           "Name: ", description, "\n",
                           sep = ""
                           ))) + 
    geom_smooth(method='lm') + geom_point(aes(color=as.factor(cluster))) + 
    theme_minimal() + labs(x=x_lab, y='Interface Delta', color='Cluster') + 
    theme(text = element_text(size=14))
}


average_ligand_distance_to_point_vs_total_score <- function(score.df, x_lab=''){
    ggplot(score.df, aes(distance_to_pocket, y=total_score,
                         text = paste(
                            "Distance to PTB: ", distance_to_pocket, "\n",
                            "Interface Delta: ", interface_delta_X, "\n",
                            "Name: ", description, "\n",
                            sep = ""
                          ))) + 
    geom_smooth(method='lm') + geom_point(aes(color=as.factor(cluster))) + 
    theme_minimal() + labs(x=x_lab, y='Total Score', color='Cluster') + 
    theme(text = element_text(size=14)) 
}

save_ggplot_to_directory <- function(plot, dir, filename){
  full_path <- file.path(dir, filename)
  full_path.html <- paste(full_path, '.html', sep='')
  #print(full_path)
  ggsave(full_path, plot, dpi=500)
  #plot.plotly <- ggplotly(plot)
  #htmlwidgets::saveWidget(as_widget(plot.plotly), full_path.html)
  return(full_path)
}


main <- function(args){
  print(args)
    #if (args$rds & args$score_file){
        message('Reading rds file')
        rmsd.df <- read_rds(args$rds, args$max_poses)
        
        message('Reading score file')
        score.df <- read_score_file(args$score_file)
        
        message('Clustering ligands')
        hcut <- hcut_rmsd_df(rmsd.df, args$max_clusters)
        
        message('Formating score file')
        score.df <- add_filepath_to_scorefile_description(score.df)
        score.df <- select_scorefile_rows_present_in_dist_matrix(score.df, rmsd.df)
        score.df <- add_cluster_field_to_score_df(score.df, hcut)
        score.df <- remove_uncontacted_ligands(score.df)
        saveRDS(score.df, file.path(args$output_dir, 'agg_results.rds'))
        
        if (dir.exists(args$pdb_dir)){
          score.df <- add_average_chain_coords_to_score_df(args$pdb_dir, score.df)
          score.df <- distance_to_binding_pocket(score.df)
          
          
          ts_dist_plot <- average_ligand_distance_to_point_vs_total_score(
            get_top_n_percent_by_metric(score.df, "total_score")
            )
          save_ggplot_to_directory(ts_dist_plot, args$output_dir, 'ts_dist_plot.png')
          id_dist_plot <- average_ligand_distance_to_point_vs_delta_x(
            get_top_n_percent_by_metric(score.df, "interface_delta_X")
            )
          save_ggplot_to_directory(id_dist_plot, args$output_dir, 'id_dist_plot.png')
        }
        
        
        
        message('Plotting results')
        id_boxplot <- interface_delta_by_cluster_boxplot(score.df)
        ts_boxplot <- total_score_by_cluster_boxplot(score.df)
        id_boxplot_best <- interface_delta_by_cluster_boxplot(
          get_top_n_percent_by_metric(score.df, "total_score")
        )
        ts_boxplot_best <- interface_delta_by_cluster_boxplot(
          get_top_n_percent_by_metric(score.df, "interface_delta_X")
        )
        save_ggplot_to_directory(id_boxplot, args$output_dir, 'id_boxplot.png')
        save_ggplot_to_directory(ts_boxplot, args$output_dir, 'ts_boxplot.png')
        save_ggplot_to_directory(id_boxplot_best, args$output_dir, 'id_boxplot_best.png')
        save_ggplot_to_directory(ts_boxplot_best, args$output_dir, 'id_boxplot_best.png')
        #ggsave(id_boxplot, file.path(args$output_dir, 'id_boxplot.png'))
        #ggsave(ts_boxplot, file.path(args$output_dir, 'ts_boxplot.png'))
        #ggsave(id_boxplot_best, file.path(args$output_dir, 't'))
        # need to add command line args to allow specifying an expected
        # binding site (x, y z) coordinates in order to make the rest of the 
        # plots

    #}else{
    #    message('Must specify both a rds file and the corresponding score file!')
    #}
}


if (!interactive()){
    parser = get_arg_parser()
    parsed_args <- parser$parse_args()
    args <- process_parsed_args(parsed_args)
    main(args)
}
# parser <- get_arg_parser()
# args <- c(
#   "-r", "/home/ethan/data/igg/rotations/Gino/corrected_hclust/0061_results/rmsd_dfs/RTX60933293", 
#   "-s", "/home/ethan/data/igg/rotations/Gino/random_docking/RTX60933293/results/score.sc",  
#   "-o",  "/home/ethan/data/igg/rotations/Gino/clust_test", 
#   "-k", "25", 
#   "-m", "1000",
#   '-p', "/home/ethan/data/igg/rotations/Gino/random_docking/RTX60933293/results")
# args <- process_parsed_args(parser$parse_args(args))
# main(args)

parallel::stopCluster(cl)


