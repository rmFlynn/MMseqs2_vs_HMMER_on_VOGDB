#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=300gb
#SBATCH --time=1-00:00:00
#SBATCH --job-name=mmseqs_sweep
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rmflynn@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/DRAM/hmmer_mmseqs2_testing_take_2/job_scripts/MMseqs2_scripts
~/miniconda3/envs/DRAM/bin/python mmseqs_from_setup_to_sweep.py
