#!/bin/sh
#SBATCH -J cellranger_align        
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=20:00:00
#SBATCH --output=log/cellranger_align.out
#SBATCH --error=log/cellranger_align.err
# Initialise the environment modules
source /etc/profile

# input
WORKING_DIR=${working_dir}
CONDA_DIR=${conda_dir}
FASTQ_DIR=${fastq_dir}
GENOME_IDX=${genome_idx}
IN_FILELIST=${in_filelist}
GTF_FP=${gtf_fp}
OUT_DIR=${out_dir}

# load anaconda

# current ID is $SLURM_ARRAY_TASK_ID 
export MRO_DISK_SPACE_CHECK=disable

SAMPLE_ID=$(awk "NR==$SLURM_ARRAY_TASK_ID" $IN_FILELIST)

cd $OUT_DIR

# Run the program
PATH=$PATH:/mmfs1/scratch/hpc/43/shihb/2024/20240110_rose_salmon_skin/cellranger/cellranger-10.0.0/bin

cellranger count --id=$SAMPLE_ID \
   --fastqs=$WORKING_DIR/$FASTQ_DIR/$SAMPLE_ID \
   --sample=$SAMPLE_ID \
   --create-bam=false \
   --transcriptome=Salmo_salar_115 
   