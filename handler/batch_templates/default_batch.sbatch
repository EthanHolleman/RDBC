#!/bin/bash -l
#SBATCH -J {}  # job name goes here
#SBATCH --partition=med
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -o {}/out_%j.txt
#SBATCH -e {}/error_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=etholleman@ucdavis.edu
#SBATCH --time=12:00:00
#SBATCH --mem 24000

{}  # command goes here
