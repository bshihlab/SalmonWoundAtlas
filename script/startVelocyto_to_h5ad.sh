#!/bin/sh
#SBATCH -J starVelocyto_to_h5d   
#SBATCH --output="log/starVelocyto_to_h5d.out"
#SBATCH --error="log/starVelocyto_to_h5d.err"
#SBATCH --time=20:00:00 
#SBATCH --mem=50G
# Initialise the environment modules
source /etc/profile



# input
WORKING_DIR=${working_dir}
SAMPLE_LIST=${in_filelist}
CONDA_DIR=${conda_dir}
STAR_DIR=${star_dir}
VELOCYTO_DIR=${velocyto_dir}

CURRENT_SAMPLE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $SAMPLE_LIST)

cd $WORKING_DIR

# load anaconda
module add anaconda3/2022.05


# Run the program
source activate $CONDA_DIR/scanpy/

python script/starVelocyto_to_Anndata.py $CURRENT_SAMPLE $STAR_DIR $VELOCYTO_DIR


conda deactivate
