#!/bin/sh

#### Settings ####
STORAGE="/mmfs1/storage/users/shihb/"
WORKING_DIR="/mmfs1/scratch/hpc/43/shihb/2024/20240110_rose_salmon_skin"
DATA_DIR="data"
FASTQ_DIR="$DATA_DIR/fastq"
SCRIPT_DIR="${WORKING_DIR}/script"
ANALYSIS_DIR="${WORKING_DIR}/analysis"
FILELIST_DIR=$DATA_DIR/filelist
FILELIST_FP=$DATA_DIR/filelist/sample.txt

CONDA_DIR="${STORAGE}/$(basename $WORKING_DIR)/conda"
LOG_DIR="${WORKING_DIR}/log"
TOOLS_DIR="${WORKING_DIR}/tools"
VELOCYTO_DIR=$ANALYSIS_DIR/velocyto
SCANPY_DIR=$ANALYSIS_DIR/scanpy

REF_GENOME_DIR="${DATA_DIR}/ref/"
REF_GENOME_IDX_DIR="${DATA_DIR}/ref/idx"

REF_STORAGE_DIR=$STORAGE/ref

GENOME_FTP="http://ftp.ensembl.org/pub/release-115/fasta/salmo_salar/dna/Salmo_salar.Ssal_v3.1.dna.toplevel.fa.gz"
GTF_FTP="http://ftp.ensembl.org/pub/release-115/gtf/salmo_salar/Salmo_salar.Ssal_v3.1.115.gtf.gz"

FILE_LIST_FP=$FILELIST_DIR/all.txt 


# derived settings
REF_FILE=$(basename $GENOME_FTP)
REF_FP=${REF_GENOME_DIR}/${REF_FILE/.gz/}

GTF_FILE=$(basename $GTF_FTP)
GTF_FP=${REF_GENOME_DIR}/${GTF_FILE/.gz/}


#### Settings

cd $WORKING_DIR

mkdir -p $STORAGE $DATA_DIR $FASTQ_DIR $REF_DIR $SCRIPT_DIR $ANALYSIS_DIR \
	$CONDA_DIR $LOG_DIR $REF_GENOME_DIR $REF_GENOME_IDX_DIR $FILELIST_DIR $VELOCYTO_DIR \
	$ANALYSIS_DIR/spaceranger $ANALYSIS_DIR/cellranger $SCANPY_DIR



