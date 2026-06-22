#!/bin/sh
#SBATCH -J star_align        
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=20:00:00
#SBATCH --output=log/star_align.out
#SBATCH --error=log/star_align.err
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
OUT_DIR=${out_dir}
BARCODE_FP="data/barcode/3M-february-2018.txt"
R1_LENGTH=${R1_length}

# load anaconda
module add anaconda3/2022.05

# current ID is $SLURM_ARRAY_TASK_ID 
SAMPLE_ID=$(awk "NR==$SLURM_ARRAY_TASK_ID" $IN_FILELIST)

CURRENT_FASTQ_DIR=$DATA_DIR/$SAMPLE_ID
R1=$(ls -m $FASTQ_DIR/${SAMPLE_ID}/*_R1*)
R2=$(ls -m $FASTQ_DIR/${SAMPLE_ID}/*_R2*)
R1=$(echo $R1 | tr -d ' ')
R2=$(echo $R2 | tr -d ' ')

OUT_DIR=$OUT_DIR/$SAMPLE_ID/
mkdir -p $(dirname $OUT_DIR) $OUT_DIR


# Run the program
source activate $CONDA_DIR/star

STAR --genomeDir $GENOME_IDX \
		--quantMode GeneCounts \
		--soloType CB_UMI_Simple \
		--sjdbGTFfile $GTF_FP \
		--soloBarcodeReadLength 0 \
		--soloCBstart 1 \
		--soloCBlen 16 \
		--soloUMIstart 17 \
		--soloUMIlen 12 \
		--runThreadN 4 \
		--soloCBwhitelist $BARCODE_FP \
		--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
		--soloUMIfiltering MultiGeneUMI_CR \
		--soloUMIdedup 1MM_CR \
		--clipAdapterType CellRanger4 \
		--outFilterScoreMin 30 \
		--outSAMtype BAM SortedByCoordinate \
		--readFilesCommand gunzip -c \
		--limitBAMsortRAM 137338953472 \
		--outSAMattributes CR UR CY UY CB UB \
		--soloFeatures Gene GeneFull Velocyto \
		--outFileNamePrefix $OUT_DIR \
		--readFilesIn $R2 $R1 \
		--soloCellReadStats Standard \
		--limitOutSJcollapsed 2000000 \
		--outFilterMatchNmin 66 \
		--outFilterMatchNminOverLread 0.4 \
		--outFilterScoreMinOverLread 0.4 \
		--soloMultiMappers EM \
		--soloCellFilter EmptyDrops_CR 		