library(parallel)
library(bio3d)
library(factoextra)
library(argparse)

get_pdb_dir_from_results_dir <- function(results_dir, 
                                         ligand_dir_to_results='results'){
    ligand_dirs <- list.dirs(path = results_dir, full.names = TRUE, recursive = FALSE)
    results.dirs <- c()
    for (ligand_dir in ligand_dirs){
        results_dirs <- c(results_dirs, file.path(ligand_dir, ligand_dir_to_results))    
    results_dirs
    }
}


get_args <- function(){
  parser <- ArgumentParser(description='Calculate RMSD across ligand poses in a pdb file collection.')
  parser$add_argument('-r', '--results_dir', help='Path to results directory produced by RDBC program.')
  parser$add_argument('-o', '--output_dir', help='Path to output directory to write all results to.')
  parser$parse_args()
}

get_pdb_paths <- function(pdb.dir){
    pdb.list <- list()
    pdb.files <- list.files(pdb.dir, pattern='\\.pdb$')
    for (f in pdb.files){
        f.abs <- file.path(pdb.dir, f)
        pdb.list[[f]] <- f.abs
    }
    return(pdb.list)
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

get_ligand_coords_from_pdb_path <- function(pdb.filepath){

    pdb.object <- bio3d::read.pdb(pdb.filepath)
    ligand.df <- subset(pdb.object$atom, chain=='X', c('x', 'y', 'z'))
    xyz_coord_df_to_vector(ligand.df)
}

get_ligand_coords_from_pdb_files <- function(pdb.filepaths.list){
    parLapply(cl, pdb.filepaths.list, get_ligand_coords_from_pdb_path)
}

calculate_pairwise_rmsd <- function(coords.list){
    ligand.rmsd.list <- list()
    for (src_coord.name in names(coords.list)){
        print(src_coord.name)
        src_coord <- coords.list[[src_coord.name]]
        rmsd.values <- parLapply(cl, coords.list, bio3d::rmsd, a=src_coord)
        rmsd.values <- as.numeric(rmsd.values)
        ligand.rmsd.list[[src_coord.name]] <- rmsd.values
    }
    df <- as.data.frame(do.call(cbind, ligand.rmsd.list))
    df
}


rmsd_for_all_ligands <- function(results_dir, output_dir, ligand_dir_to_results='results'){
    print(results_dir)
    ligand_dirs <- list.dirs(results_dir, full.names=TRUE, recursive=FALSE)
    rmsd.list <- list()
    output_dir.rmsd <- file.path(output_dir, 'rmsd_dfs')
    if (! dir.exists(output_dir.rmsd)){
      dir.create(output_dir.rmsd)
    }
    for (dir in ligand_dirs){
        dir.pdb <- file.path(dir, ligand_dir_to_results)
        output_prefix <- file.path(output_dir.rmsd, basename(dir))
        pdb.paths <- get_pdb_paths(dir.pdb)
        if (length(pdb.paths) > 0){
        ligand.coords <- get_ligand_coords_from_pdb_files(pdb.paths)
        ligand.coords.rmsd <- calculate_pairwise_rmsd(ligand.coords)
        rmsd.list[[basename(dir)]] <- ligand.coords.rmsd
        saveRDS(ligand.coords.rmsd, output_prefix)
        }
    }
    rmsd.list
}

hclust_from_rmsd_df <- function(rmsd.df, output_dir, ligand.name){
    
    gap <- fviz_nbclust(rmsd.df, hcut, 'gap_stat', k.max=nrow(rmsd.df)-1)
    n_cuts <- match(max(gap$data$gap), gap$data$gap)
    hcut <- hcut(rmsd.df, k=n_cuts, stand = TRUE, hc_method='complete')
    saveRDS(hcut, file.path(output_dir, ligand.name))
}

hclust_all_ligands <- function(rmsd.list, output_dir){
    hclust_dir <- file.path(output_dir, 'hclust')
    if (!dir.exists(hclust_dir)){
      dir.create(hclust_dir)
    }
    hclust.rep <- rep(hclust_dir, length(rmsd.list))
    #parLapply(cl, rmsd.list, hclust_from_rmsd_df, ligand.name=names(rmsd.list))
              #output_path=hclust.rep)
    for (ligand.name in names(rmsd.list)){
      rmsd.df <- rmsd.list[[ligand.name]]
      hcut <- hclust_from_rmsd_df(rmsd.df, hclust_dir, ligand.name)
    }
}

rmsd_hclust <- function(rmsd.df){
    nbcl <- fviz_nbclust(rmsd.df, hcut, k.max = nrow(rmsd.df)-1)
    
}

main <- function(args){
    results_dir <- args$results_dir
    output_dir <- args$output_dir
    rmsd.list <- rmsd_for_all_ligands(results_dir, output_dir)
    hclust_all_ligands(rmsd.list, output_dir)
}

cl <- makeCluster(detectCores(), type='FORK')

if (!interactive()){
    main(get_args())
    stopCluster(cl)
}