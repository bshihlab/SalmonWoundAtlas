#!/bin/sh
#SBATCH -J conda_env_create              
#SBATCH --time=01:00:00 
#SBATCH --mem=20G
# Initialise the environment modules
source /etc/profile

# input
WORKING_DIR=${working_dir}
YML_FP=${yml_fp}
CONDA_ENVNAME=${conda_envname}
CONDA_DIR=${conda_dir}



# load anaconda
module add anaconda3/2022.05

# Run the program
conda env create --prefix $CONDA_DIR/$CONDA_ENVNAME -f $YML_FP

