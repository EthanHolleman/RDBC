import argparse
from pathlib import Path


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--ligands', help='Path to ligand directory.')
    parser.add_argument('-p', '--protein', help='Path to protein (pdb) file without ligand.')
    parser.add_argument('-o', '--parent', help='Parent directory that all results are written under.')

    a = parser.parse_args()
    a.parent = Path(a.parent)
    return a