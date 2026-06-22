#!/bin/sh
#SBATCH -J cellbender        
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=20:00:00
#SBATCH --output=log/cellbender.out
#SBATCH --error=log/cellbender.err
# Initialise the environment modules
source /etc/profile

# input
WORKING_DIR=${working_dir}
CONDA_DIR=${conda_dir}
IN_FILELIST=${in_filelist}
IN_DIR=${in_dir}
OUT_DIR=${out_dir}


# load anaconda
module add anaconda3/2022.05

# current ID is $SLURM_ARRAY_TASK_ID 
SAMPLE_LINE=($(awk "NR==$SLURM_ARRAY_TASK_ID" $IN_FILELIST))
SAMPLE_ID=${SAMPLE_LINE[0]}
ESTIMATED_CELL_NUM=${SAMPLE_LINE[1]}
DROPLET_INCLUDE=${SAMPLE_LINE[2]}


#IN_DIR=$IN_DIR/$SAMPLE_ID/Solo.out/Gene/raw/
OUT_DIR=$OUT_DIR/${SAMPLE_ID}
mkdir -p $OUT_DIR
cd $OUT_DIR

IN_FP=${IN_DIR}/${SAMPLE_ID}.h5ad
OUT_FP=${SAMPLE_ID}.h5


# Run the program
source activate $CONDA_DIR/cellbender

#cellbender remove-background \
#                 --input $IN_DIR \
#                 --output $OUT_FP \
#				 --cpu-threads 4

cellbender remove-background --cpu-threads 8 \
	--input $IN_FP --output $OUT_FP --expected-cells $ESTIMATED_CELL_NUM 
	--total-droplets-included $DROPLET_INCLUDE
