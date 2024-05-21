#!/bin/bash
#SBATCH --job-name=ray_gen_halo_667
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=1000mb
#SBATCH --partition=main
#SBATCH --qos=main

module load python/python3/3.9.6
module load cmake/3.26.3_gcc_11.1.0

python spectrauploader.py
