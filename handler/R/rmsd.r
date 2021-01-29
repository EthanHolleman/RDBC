library(bio3d)
library(argparse)
library(factoextra)
library(reshape)
library(tidyr)
library(parallel)

# number_cores <- detectCores()

get_args <- function(){
  parser <- ArgumentParser(description='Calculate RMSD across ligand poses in a pdb file collection.')
  parser$add_argument('-r', '--results_dir', help='Path to results directory produced by RDBC program.')
  parser$add_argument('-o', '--output_dir', help='Path to output directory to write all results to.')
  parser$parse_args()
}


read_pdb_dir_as_list <- function(pdb.dir){
  files <- list.files(pdb.dir, pattern = "\\.pdb$")
  pdb.list <- list()
  i <- 1
  for (f in files){
    f.abs <- paste(pdb.dir, f, sep='/')
    pdb.list[[f]] <- read.pdb(f.abs)
    message(paste('Pdb files read:', i))
    i <- i + 1
  }
  return(pdb.list)
}

extract_ligand_coords_from_pdb_object <- function(pdb.object, chain_name='X'){
  ligand.df <- subset(pdb.object$atom, chain==chain_name, c('x', 'y', 'z'))
}

xyz_coord_df_to_vector <- function(xyz.df){
  coords <- c()
  for (i in 1:nrow(xyz.df)){
    for (j in 1:3){
      coords <- c(coords, xyz.df[i, j])
    }
  }
  assertthat::are_equal(length(coords), nrow(xyz.df)*ncol(xyz.df))
  coords
}

get_ligand_coord_vector_from_pdb <- function(pdb.object, chain_name='X'){
  ligand.df <- extract_ligand_coords_from_pdb_object(pdb.object, chain_name)
  xyz_coord_df_to_vector(ligand.df)
}

pairwise_ligand_rmsd <- function(protein.list){
  
  c <- 1
  pairwise.rmsd <- list()
  for (name_a in names(protein.list)){
    message(name_a)
    pairwise.rmsd[[name_a]] <- c()
    for (name_b in names(protein.list)){
      # Need to compare only the ligand coords
      ligand.a <- protein.list[[name_a]]
      ligand.b <- protein.list[[name_b]]
      ligand.a.coords <- get_ligand_coord_vector_from_pdb(ligand.a)
      ligand.b.coords <- get_ligand_coord_vector_from_pdb(ligand.b)
      rmsd <- rmsd(ligand.a.coords, ligand.b.coords)
      pairwise.rmsd[[name_a]] <- c(pairwise.rmsd[[name_a]], rmsd)
    }
  }
  pairwise.rmsd.df <- do.call(cbind, pairwise.rmsd)
  pairwise.rmsd.df <- as.data.frame(pairwise.rmsd.df)
  rownames(pairwise.rmsd.df) <- names(protein.list)
  return(pairwise.rmsd.df)
}

write_pairwise_rmsd_df_to_csv <- function(pairwise.rmsd.df, filepath){
  print(filepath)
  write.csv(filepath, pairwise.rmsd.df, quote = FALSE)
}

calculate_rmsd_for_run <- function(results_dir, output_dir, ligand_dir_to_results='results'){
  ligand_dirs <- list.dirs(path = results_dir, full.names = TRUE, recursive = FALSE)
  for (ligand_dir in ligand_dirs){
    print(ligand_dir)
    ligand_dir.results <- paste(ligand_dir, ligand_dir_to_results, sep='/')
    message('Reading PDB list')
    pdb.list <- read_pdb_dir_as_list(ligand_dir.results)
    print(length(pdb.list))
    message('Calculating RMSD')
    pairwise.rmsd.df <- pairwise_ligand_rmsd(pdb.list)
    output.path <- paste(output_dir, basename(ligand_dir), sep = '')
    output.path <- paste(output.path, '.csv', sep = '')
    write.csv(pairwise.rmsd.df, output.path, , quote = FALSE)
    message(output.path)
  }
}

read_rmsd_file <- function(filepath){
  as.data.frame(read.csv(filepath, row.names = 1))
}

rmsd_df_to_hcut <- function(rmsd.df){
  # get optimal number of groups based on gap statistic
  gap <- fviz_nbclust(rmsd.df, hcut, 'gap_stat')
  n_cuts <- match(max(gap$data$gap), gap$data$gap)
  hcut <- hcut(rmsd.df, k=n_cuts, stand = TRUE)
}

main <- function(command_args){
  calculate_rmsd_for_run(command_args$results_dir, command_args$output_dir)
}

if (!(interactive())){
  main(get_args())
}

