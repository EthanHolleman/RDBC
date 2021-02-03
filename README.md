# RDBC
Rosetta ligand docking batch job submission control and organizer. Created in order to make submitting large numbers of 
Rosetta simulations to SLURM workload manager easier and more organized.

RDBC uses templates for Rossetta xml scripts and options files as well as for batch scripts which makes running many similar 
simulations easier as new templates can be created quickly. Templates are filled where `{}` are found from top to 
bottom, so currently the order of `{}` matters. Use the templates in any of the `template` directories as examples.

## Docking simulations

`python3 rh.py --help`

```
usage: rh.py [-h] [-l LIGANDS] [-p PROTEIN] [-o PARENT] [-e EXE] [-m] [-d] [-i NUM_ITERS] [-x XML_TEMPLATE] [-op OPTIONS_TEMPLATE] [-b BATCH_TEMPLATE] [-a AGGREGATE_RESULTS_PATH]
             [-f AGGREGATED_FILEPATH] [-pd TARGET_PROTEIN_DIR]

optional arguments:
  -h, --help            show this help message and exit
  -l LIGANDS, --ligands LIGANDS
                        Path to ligand directory.
  -p PROTEIN, --protein PROTEIN
                        Path to protein (pdb) file (or directory) without ligand.
  -o PARENT, --parent PARENT
                        Parent directory that all results are written under.
  -e EXE, --exe EXE     Path to RosettaScripts exe
  -m, --moist           Moist run. Do everything normally but only submit the first job.
  -d, --dry             Dry run. Do everything except submit jobs.
  -i NUM_ITERS, --num_iters NUM_ITERS
                        Number of iterations to run for each structure. Default 2000.
  -x XML_TEMPLATE, --xml_template XML_TEMPLATE
                        Path to Rosetta script XML template.
  -op OPTIONS_TEMPLATE, --options_template OPTIONS_TEMPLATE
                        Path to Rosetta options script to use for all jobs.
  -b BATCH_TEMPLATE, --batch_template BATCH_TEMPLATE
                        Path to batch script template to use.
  -a AGGREGATE_RESULTS_PATH, --aggregate_results_path AGGREGATE_RESULTS_PATH
                        Path to results from a completed job that you want to aggregate the results from.
  -f AGGREGATED_FILEPATH, --aggregated_filepath AGGREGATED_FILEPATH
                        If using -a, path to file to write aggregated results to.
  -pd TARGET_PROTEIN_DIR, --target_protein_dir TARGET_PROTEIN_DIR
                        If using -a, path to directory containing all target proteins. Including this option will write the target protein path as a field in the aggregated results.
```

## R scripts

The repo also contains a few R scripts used for analysis and graphing after docking simulations are completed.
If they have a command line interface they will be listed below, otherwise the script is intended to be
used interactively. Scripts with command line interfaces can still be used interactively (usually easier to
troubleshoot).

### rmsd_p.r

```
usage: rmsd_p.r
       [-h] [-r RESULTS_DIR] [-o OUTPUT_DIR]

Calculate RMSD across ligand poses in a pdb file collection making use of all
available processor cores. Save RMSD dataframe as rds file and preform
heirarchial clustering using the factoextra package.

optional arguments:
  -h, --help            show this help message and exit
  -r RESULTS_DIR, --results_dir RESULTS_DIR
                        Path to results directory produced by RDBC program.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Path to output directory to write all results to.
 ```
 
 ### dendro.r 
 
 Work in progress, still need to add some stuff to produce all plots and test.
 
 ```
 usage: dendro.r [-h] [-r RDS] [-s SCORE_FILE] [-m MAX_POSES] [-o OUTPUT_DIR]

Analyize Rosetta ligand RMSD based heirarhical clustering.

optional arguments:
  -h, --help            show this help message and exit
  -r RDS, --rds RDS     Path to RDS file containing a RMSD distance matrix for
                        all ligand poses.
  -s SCORE_FILE, --score_file SCORE_FILE
                        Path to score file containing results of all ligand
                        docking simulations contained by the file specified
                        with the -r argument. This should be formated as a tsv
                        file from the original Rosetta output.
  -m MAX_POSES, --max_poses MAX_POSES
                        Max number of poses to evaluate. Default 500.
                        Increasing this greatly increases computation time.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Folder to save plots to.
  ```
  
 
 

