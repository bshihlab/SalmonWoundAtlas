#!/bin/sh
#SBATCH -J install_scanpy              
#SBATCH --time=10:00:00 
#SBATCH --mem=20G
# Initialise the environment modules
source /etc/profile

# input

# Run the program
module add anaconda3/2022.05

cd "/mmfs1/scratch/hpc/43/shihb/2024/20240110_rose_salmon_skin"

source activate conda/scanpy

pip install scanpy
conda env export | grep -v "^prefix: " > analysis/scanpy_environment.yml

conda deactivate

