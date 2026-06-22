#!/bin/sh
#SBATCH -J scanpy        
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=20:00:00
#SBATCH --output=log/scanpy.out
#SBATCH --error=log/scanpy.err
# Initialise the environment modules
source /etc/profile

# input
cd /mmfs1/scratch/hpc/43/shihb/2024/20240110_rose_salmon_skin

# load anaconda
module add anaconda3/2022.05

# Run the program
source activate /mmfs1/storage/users/shihb//20240110_rose_salmon_skin/conda/scanpy

python script/scanpy_toh5ad.py