# RDBC
Rosetta ligand docking batch job submission control and organizer

## Usage guide

`python3 rh.py --help`

```
usage: rh.py [-h] [-l LIGANDS] [-p PROTEIN] [-o PARENT] [-e EXE] [-m] [-d]
             [-i NUM_ITERS] [-x XML_TEMPLATE]

optional arguments:
  -h, --help            show this help message and exit
  -l LIGANDS, --ligands LIGANDS
                        Path to ligand directory.
  -p PROTEIN, --protein PROTEIN
                        Path to protein (pdb) file without ligand.
  -o PARENT, --parent PARENT
                        Parent directory that all results are written under.
  -e EXE, --exe EXE     Path to RosettaScripts exe
  -m, --moist           Moist run. Do everything normally but only submit the
                        first job.
  -d, --dry             Dry run. Do everything except submit jobs.
  -i NUM_ITERS, --num_iters NUM_ITERS
                        Number of iterations to run for each structure.
                        Default 2000.
  -x XML_TEMPLATE, --xml_template XML_TEMPLATE
                        Path to Rosetta script XML template.
```

