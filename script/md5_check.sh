#!/bin/sh
#SBATCH -J md5sum        
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=20:00:00
#SBATCH --output=log/md5sum.out
#SBATCH --error=log/md5sum.err
# Initialise the environment modules
source /etc/profile

# input
IN_FILELIST=${in_filelist}
WORKING_DIR=${working_dir}

mkdir -p md5_results
IN_FP=$(awk "NR==$SLURM_ARRAY_TASK_ID" $IN_FILELIST)
OUT_FP=md5_results/${SLURM_ARRAY_TASK_ID}.txt

cd $WORKING_DIR

md5sum $IN_FP > $OUT_FP
