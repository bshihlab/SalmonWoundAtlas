#!/bin/sh
#SBATCH -J star_genome_index              
#SBATCH --time=20:00:00 
#SBATCH --mem=80G
#SBATCH --output=log/star_genome_index.out
#SBATCH --error=log/star_genome_index.err
# Initialise the environment modules
source /etc/profile

# input
WORKING_DIR=${working_dir}
CONDA_DIR=${conda_dir}
DATA_DIR=${data_dir}
REF_FP=${ref_fp}
GTF_FP=${gtf_fp}


# load anaconda
module add anaconda3/2022.05

# Run the program
source activate $CONDA_DIR/star

STAR --runMode genomeGenerate \
		--runThreadN 1 \
		--genomeDir $DATA_DIR/ref/idx/ \
		--genomeFastaFiles $REF_FP \
		--sjdbGTFfile $GTF_FP

