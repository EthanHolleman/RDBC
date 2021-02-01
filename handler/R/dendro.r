# dendrogram creation and plotting 
library(parallel)
library(bio3d)
library(factoextra)
library(argparse)
library(plotly)

get_args <- function(){
  parser <- ArgumentParser(description='Analyize Rosetta ligand RMSD based heirarhical clustering.')
  parser$add_argument('-r', '--rds', help='Path to RDS file containing a RMSD distance matrix for all ligand poses.')
  parser$add_argument('-s', '--score_file', help='Path to score file containing results of all ligand docking simulations contained by the file specified with the -r argument. This should be formated as a tsv file from the original Rosetta output.')    
  parser$add_argument('-m', '--max_poses', default=500, help='Max number of poses to evaluate. Default 500. Increasing this greatly increases computation time.')
  parser$add_argument('-o', '--output_dir', default='', help='Folder to save plots to.')
  parser$parse_args()
}

cl <- parallel::makeCluster(4)

read_rds <- function(rds.path, max_poses=500){
    rmsd.df <- readRDS(rds.path)
    if (nrow(rmsd_dfs) > max_poses){
        rmsd.df <- rmsd.df[1:max_poses, 1:max_poses]
    }
    rownames(rmsd.df <- names(rmsd.df))  # should be a square distance matrix
    return(rmsd.df)
}

read_score_file <- function(score_file.path){
    return(read.delim(score_file.path))
}

hcut_rmsd_df <- function(rmsd.df){
    gap <- fviz_nbclust(rmsd.df, hcut, 'gap_stat', k=25)
    n_cuts <- match(max(gap$data$gap), gap$data$gap)
    hcut <- hcut(rmsd.df, k=n_cuts, hc_method='complete')
    return(hcut)
}

# Scorefile description field will not have the .pdb extension but this
# extension should be present as the row / col names of the distance
# matrix read using the read_rds file. In order to match these two data
# frames up we add the the description + /pdb extension to a new field
# called "description_pdb"
add_filepath_to_scorefile_description <- function(score.df){
    descrip_with_ext <- lapply(score_file$description, function(s){paste(s, '.pdb', sep='')})
    descrip_with_ext <- unlist(descrip_with_ext)
    score.df$description_pdb <- descrip_with_ext
    return(score.df)
}


select_scorefile_rows_present_in_dist_matrix <- function(score.df, rmsd.df){
    score.df.subset <- subset(score.df, description_pdb %in% rownames(rmsd.df))
    return(score.df.subset)
}


add_cluster_field_to_score_df <- function(score.df, hcut){
    score.df <- score.df[order(names(cut$cluster)), ]
    score.df$cluster <- hcut$cluster
}


remove_uncontacted_ligands <- function(score.df){
    score.df.contact <- subset(sample.subset, ligand_is_touching_X==1)
    return(score.df.contact)
}


interface_delta_by_cluster_boxplot <- function(score.df){
    ggplot(sample.subset, aes(x=as.factor(cluster), y=interface_delta_X)) + 
    geom_boxplot(fill='#8b21a3') + 
    theme_minimal() + labs(x='Cluster', y='Interface Delta')
}

total_score_by_cluster_boxplot <- function(score.df){
    ggplot(sample.subset, aes(x=as.factor(cluster), y=total_score)) + 
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
    if (basename(f) %in% rmsd.scores$description_pdb){
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


average_ligand_distance_to_point_vs_delta_x <- function(score.df, x_lab=''){
    ggplot(score.df, aes(distance_to_pocket, y=interface_delta_X)) + 
    geom_smooth(method='lm') + geom_point(aes(color=as.factor(cluster))) + 
    theme_minimal() + labs(x=x_lab, y='Interface Delta', color='Cluster') + 
    theme(text = element_text(size=14))
}

average_ligand_distance_to_point_vs_total_score <- function(score.df, x_lab=''){
    ggplot(score.df, aes(distance_to_pocket, y=total_score)) + 
    geom_smooth(method='lm') + geom_point(aes(color=as.factor(cluster))) + 
    theme_minimal() + labs(x=x_lab, y='Total Score', color='Cluster') + 
    theme(text = element_text(size=14)) 
}


main <- function(args){
    if (args$rds & args$score_file){
        rmsd.df <- read_rds(args$rds, args$max_poses)
        score.df <- read_score_file(args$score_file)
        hcut <- hcut_rmsd_df(rmsd.df)
        score.df <- add_filepath_to_scorefile_description(score.df)
        score.df <- select_scorefile_rows_present_in_dist_matrix(score.df, rmsd.df)
        score.df <- add_cluster_field_to_score_df(score.df, hcut)
        score.df <- remove_uncontacted_ligands(score.df)
        id_boxplot <- interface_delta_by_cluster_boxplot(score.df)
        ts_boxplot <- total_score_by_cluster_boxplot(score.df)
        ggsave(id_boxplot, file.path(args$output_dir, 'id_boxplot.png'))
        ggsave(ts_boxplot, file.path(args$output_dir, 'ts_boxplot.png'))
        # need to add command line args to allow specifying an expected
        # binding site (x, y z) coordinates in order to make the rest of the 
        # plots

    }else{
        message('Must specify both a rds file and the corresponding score file!')
    }


}


if (!interactive()){
    main(get_args())
}


parallel::stopCluster(cl)


