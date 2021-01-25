import argparse
from pathlib import Path
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
    a = parser.parse_args()
    a.parent = Path(a.parent)
    a.protein = Path(a.protein)
    a.exe = Path(a.exe)
    if not a.parent.is_dir():
        a.parent.mkdir()
    return a