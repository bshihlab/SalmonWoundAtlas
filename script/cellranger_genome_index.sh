#!/bin/bash
#SBATCH -J cellranger_index
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=20:00:00
#SBATCH --output=log/spaceranger_idx.out
#SBATCH --error=log/spaceranger_idx.err
# Initialise the environment modules
source /etc/profile

# input
WORKING_DIR=${working_dir}
CONDA_DIR=${conda_dir}
DATA_DIR=${data_dir}
REF_FP=${ref_fp}
GTF_FP=${gtf_fp}

mkdir -p $DATA_DIR/ref_cellranger

#module add spaceranger/2.0.1
PATH=$PATH:/mmfs1/scratch/hpc/43/shihb/2024/20240110_rose_salmon_skin/cellranger/cellranger-10.0.0/bin


# Run the program
cellranger mkgtf \
    $GTF_FP \
    ${GTF_FP/.gtf/_filtered.gtf} \
	--attribute=gene_biotype:protein_coding \
	--attribute=gene_biotype:lincRNA \
	--attribute=gene_biotype:antisense
	
	
cellranger mkref \
    --genome=Salmo_salar_115 \
    --fasta=$REF_FP \
    --genes=${GTF_FP/.gtf/_filtered.gtf} \
	--nthreads 8