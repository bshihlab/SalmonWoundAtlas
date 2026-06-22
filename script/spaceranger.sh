#!/bin/bash
#SBATCH -J spaceranger
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=20:00:00
#SBATCH --output=log/spaceranger.out
#SBATCH --error=log/spaceranger.err
# Initialise the environment modules
source /etc/profile

# input
WORKING_DIR=${working_dir}
FASTQ_DIR=${fastq_dir}
GENOME_IDX=${genome_idx}
IN_FILELIST=${in_filelist}
OUT_DIR=${out_dir}


# current ID is $SLURM_ARRAY_TASK_ID 
SAMPLE_LINE=($(awk "NR==$SLURM_ARRAY_TASK_ID" $IN_FILELIST))
SPATIAL_DIR="data/spatial/"
SAMPLE_ID=${SAMPLE_LINE[0]}
FASTQ_DIR=${FASTQ_DIR}/${SAMPLE_LINE[0]}
IMAGE=${SPATIAL_DIR}/image/${SAMPLE_LINE[1]}
SLIDE=${SAMPLE_LINE[2]}
SLIDE_FILE=${SPATIAL_DIR}/gpr/${SAMPLE_LINE[3]}
AREA=${SAMPLE_LINE[4]}

OUT_DIR=$OUT_DIR/$SAMPLE_ID/
mkdir -p $OUT_DIR


# Run the program
export MRO_DISK_SPACE_CHECK=disable

#module add spaceranger/2.0.1
PATH=$PATH:/mmfs1/scratch/hpc/43/shihb/2024/20240110_rose_salmon_skin/spaceranger/spaceranger-4.1.0/bin


spaceranger count \
   --id=$SAMPLE_ID \
   --transcriptome=Salmo_salar_Ssal_v3_1 \
   --fastqs=$FASTQ_DIR \
   --sample=$SAMPLE_ID \
   --image=$IMAGE \
   --slide=$SLIDE \
   --create-bam=false \
   --slidefile=$SLIDE_FILE \
   --area=$AREA \
   --reorient-images true \
   --output-dir $OUT_DIR
