#!/bin/bash
#SBATCH -J spaceranger
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
GENOME_FP=${genome_fp}
GTF_FP=${gtf_fp}

# load anaconda
export MRO_DISK_SPACE_CHECK=disable

#module add spaceranger/2.0.1
PATH=$PATH:/mmfs1/scratch/hpc/43/shihb/2024/20240110_rose_salmon_skin/spaceranger/spaceranger-4.1.0/bin

spaceranger mkgtf $GTF_FP ${GTF_FP/.gtf/_filtered.gtf} \
                   --attribute=gene_biotype:protein_coding \
                   --attribute=gene_biotype:lincRNA \
                   --attribute=gene_biotype:antisense 
				   
spaceranger mkref --genome=Salmo_salar_Ssal_v3_1 --fasta=$GENOME_FP --genes=${GTF_FP/.gtf/_filtered.gtf} --nthreads 8