#!/bin/sh
#SBATCH -J sra_download                
#SBATCH --output="log/sra_download.out"
#SBATCH --error="log/sra_download.err"
#SBATCH --time=20:00:00 
#SBATCH --mem=8G
# Initialise the environment modules
source /etc/profile



# input
WORKING_DIR=${working_dir}
CONDA_DIR=${conda_dir}
FASTQ_DIR=${fastq_dir}
SAMPLE_LIST=${in_filelist}

CURRENT_SAMPLE=$(awk "NR==$SLURM_ARRAY_TASK_ID" $SAMPLE_LIST)

cd $WORKING_DIR

# load anaconda
module add anaconda3/2022.05


# Run the program
source activate /mmfs1/storage/users/shihb/20231213_dalmatian_genetics/conda/sra_tools/

mkdir -p $(dirname $FASTQ_DIR) $FASTQ_DIR $FASTQ_DIR/$CURRENT_SAMPLE

if [ ! -f "${FASTQ_DIR}/${CURRENT_SAMPLE}_1.fastq.gz" ]; then
        prefetch --max-size 100g $CURRENT_SAMPLE && vdb-validate $CURRENT_SAMPLE && fastq-dump --split-files --origfmt --gzip $CURRENT_SAMPLE -O $FASTQ_DIR/$CURRENT_SAMPLE
fi

conda deactivate
