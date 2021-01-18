import argparse
from pathlib import Path


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--ligands', help='Path to ligand directory.')
    parser.add_argument('-p', '--protein', help='Path to protein (pdb) file without ligand.')
    parser.add_argument('-o', '--parent', help='Parent directory that all results are written under.')
    parser.add_argument('-e', '--exe', help='Path to RosettaScripts exe')
    parser.add_argument('-m', '--moist', action='store_true', help='Moist run. Do everything normally but only submit the first job.')
    parser.add_argument('-d', '--dry', action='store_true', help='Dry run. Do everything except submit jobs.')
    parser.add_argument('-i', '--num_iters', default=2000, type=int, help='Number of iterations to run for each structure. Default 2000.')

    a = parser.parse_args()
    a.parent = Path(a.parent)
    a.protein = Path(a.protein)
    a.exe = Path(a.exe)
    if not a.parent.is_dir():
        a.parent.mkdir()
    return a