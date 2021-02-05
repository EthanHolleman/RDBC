import argparse
from pathlib import Path
import sys
from handler import DEFAULT_BATCH, DEFAULT_OPTIONS, DEFAULT_XML


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--ligands', help='Path to ligand directory.')
    parser.add_argument('-p', '--protein', help='Path to protein (pdb) file (or directory) without ligand.')
    parser.add_argument('-o', '--parent', help='Parent directory that all results are written under.')
    parser.add_argument('-e', '--exe', help='Path to RosettaScripts exe')
    parser.add_argument('-m', '--moist', action='store_true', help='Moist run. Do everything normally but only submit the first job.')
    parser.add_argument('-d', '--dry', action='store_true', help='Dry run. Do everything except submit jobs.')
    parser.add_argument('-i', '--num_iters', default=2000, type=int, help='Number of iterations to run for each structure. Default 2000.')
    parser.add_argument('-x', '--xml_template', default=DEFAULT_XML, help='Path to Rosetta script XML template.')
    parser.add_argument('-op', '--options_template', default=DEFAULT_OPTIONS, help='Path to Rosetta options script to use for all jobs.')
    parser.add_argument('-b', '--batch_template', default=DEFAULT_BATCH, help='Path to batch script template to use.')
    parser.add_argument('-a', '--aggregate_results_path', default=False, help='Path to results from a completed job that you want to aggregate the results from.')
    parser.add_argument('-f', '--aggregated_filepath', default='agg_results.tsv', help='If using -a, path to file to write aggregated results to.')
    parser.add_argument('-pd', '--target_protein_dir', help='If using -a, path to directory containing all target proteins. Including this option will write the target protein path as a field in the aggregated results.')
    parser.add_argument('-mi', '--multi_iterations', type=int, help='Run the same experiment multiple times. Used for random docking to get more samples (poses).')
    args = parser.parse_args()

    if args.aggregate_results_path:
        return args
    
    if not args.parent:
        print('Please specify a parent directory to write results to (-o)')
        sys.exit()

    args.parent = Path(args.parent)
    args.protein = Path(args.protein)
    args.exe = Path(args.exe)

    if not args.parent.is_dir():
        args.parent.mkdir()

    return args
