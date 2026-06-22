##### scRNAseq salmon wound healing
WORKING_DIR="/mmfs1/scratch/hpc/43/shihb/2024/20240110_rose_salmon_skin"
mkdir -p $WORKING_DIR
cd $WORKING_DIR
source ${WORKING_DIR}/script/_full_process_settings.sh
module add anaconda3/2022.05
module add vim/8.1


#### Set up tools and environment ####
# conda env create from yml files
COND_ENVS=("biomart" "spaceranger")
for CURRENT_ENV in ${COND_ENVS[@]}; do 
sbatch \
--output=${LOG_DIR}/env_create_${CURRENT_ENV}.out \
--error=${LOG_DIR}/env_create_${CURRENT_ENV}.err \
--export=working_dir=$WORKING_DIR,\
yml_fp=$SCRIPT_DIR/yml/${CURRENT_ENV}.yml,\
conda_envname=${CURRENT_ENV},\
conda_dir=$CONDA_DIR \
$SCRIPT_DIR/conda_env_create.sh
done

# yml didn't work, Used pip install instead on 22/03/24
module add anaconda3/2022.05
conda create -p $STORAGE/20240110_rose_salmon_skin/conda/cellbender python=3.7
source activate $STORAGE/20240110_rose_salmon_skin/conda/cellbender 
pip install cellbender 
conda deactivate

# download spaceranger
curl -o spaceranger-4.1.0.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-4.1.0.tar.gz?Expires=1777592899&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=R2NBRE0W-F0Ich9eK0qVIHEdIlAQTyGYeK~8XbVgXfZOQDCnXG97bc5BFhKEM20jlklDugpSLenQR5rgDos3Ndz3u63HyUNhO779ZtpI0bDMmppaqdtiAiTPfgrE7p8DBDCGclNTL1BYUhpE6tjokbtwNfCn5ykfgFZsnbh8nI2qkmCJrD1WuFDNZARrj15Xhy0OPDDt4ConqjBc9e3EWBt0Zz4u~qdg3IXXJ1eOL87WWaDJZX~DwiMzXFwimeaHkh6jojErhz3ZzXEDd59oT-3EPznj7kIRTbZG9Q3TtPDMiJ5ypbczM3nmD3n1xiG5s2WEtLjo00BP3iCGhi15EA__"
tar -xzf spaceranger-4.1.0.tar.gz
mkdir spaceranger
cp -r spaceranger-4.1.0 spaceranger/  # HEC delete the file if the file date is not within 30 days

# download cellranger
curl -o cellranger-10.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz?Expires=1777591716&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=jTleOGxHKqteASE0-bUmNEn4aI8B88g6k27yXmgf3PDVPC3yzPRLhYQzbT66zZ3PBDRuGl-xQuOvdRGPJ9h66x~w7BrSdckSzEV5-eKph3Z4Hc27QBYTWCGF3ecD66JE6s5qz~cbDmUd-kZ2w3xR2-Ba44ABsai3i55QD1c5nv8-77tC6TMVThzor2J4K~DWoFM7V5dX3UfaYt5MBHjCScHz4HRbtkRbLDI~1BwndBWmIdIwybKRpA5VWlOFxaPIKi8DviGwGoRLMWIUavLFoS6B9Cc8cIHHvUctvLXinBwM9c-w0YMgm9qzjCpFUZbGurrxZyAW91vk~RazOBFS7A__"
tar -xzf cellranger-10.0.0.tar.gz
mkdir cellranger
cp -r cellranger-10.0.0 cellranger/  # HEC delete the file if the file date is not within 30 days



# Used pip install instead on 08/11/24
module add anaconda3/2022.05
conda create -p conda/scanpy python=3.11
pip install -U loompy
pip install -U scvelo
conda install -c conda-forge scikit-image # version 
conda deactivate

# or use 
qsub script/scanpy_install.sh


#### Download reference
curl $GENOME_FTP --output $REF_GENOME_DIR/$(basename $GENOME_FTP)
curl $GTF_FTP --output $REF_GENOME_DIR/$(basename $GTF_FTP)
cd $REF_GENOME_DIR
gunzip $(basename $GENOME_FTP)
gunzip $(basename $GTF_FTP)


#### Download samples
# Copy over samples from LUNA

#### Download samples from PRJNA970284: SRR24471129, SRR24471130, SRR24471131, SRR24471132
echo SRR24471129 > data/sra_list.txt
echo SRR24471130 >> data/sra_list.txt

sbatch \
--array=1-2 \
--output=${LOG_DIR}/download_sra.out \
--error=${LOG_DIR}/download_sra.err \
--export=working_dir=$WORKING_DIR,\
conda_dir=$CONDA_DIR,\
fastq_dir=$FASTQ_DIR,\
in_filelist=data/sra_list.txt \
$SCRIPT_DIR/download_sra.sh

# change file name so it follows the same format as the other files (file 1 and 2 are index)
for file in data/fastq/SRR*/*_3.fastq.gz; do
mv $file ${file/_3.fastq.gz/_R1_001.fastq.gz}
done
for file in data/fastq/SRR*/*_4.fastq.gz; do
mv $file ${file/_4.fastq.gz/_R2_001.fastq.gz}
done


#### Generate genome index
sbatch --export=working_dir=$WORKING_DIR,\
	data_dir=$DATA_DIR,\
	ref_fp=$REF_FP,\
	gtf_fp=$GTF_FP,\
	conda_dir=$CONDA_DIR \
	$SCRIPT_DIR/star_genome_index.sh


#### Align with STARsolo
sbatch \
--array=1-7 \
--export=working_dir=$WORKING_DIR,\
fastq_dir=$DATA_DIR/fastq,\
genome_idx=$DATA_DIR/ref/idx,\
in_filelist=filelist/snRNAseq_filelist.txt,\
out_dir=$ANALYSIS_DIR/star,\
gtf_fp=$GTF_FP,\
conda_dir=$CONDA_DIR \
$SCRIPT_DIR/star_align.sh


#### gzip star solo outputs 
cp -r analysis/star ./
for CURRENT_DIR in star/*; do
echo $CURRENT_DIR
cd $CURRENT_DIR/Solo.out/GeneFull/raw/
for FILE in ./*.*; do
echo $FILE
gzip $FILE 
done
cd -
done


#### Run scanpy
sbatch script/scanpy_run.sh


#### Remove background with cellbender
sbatch \
--array=1-7 \
--export=working_dir=$WORKING_DIR,\
in_filelist=filelist/snRNAseq_filelist_cellbender.txt,\
in_dir=$ANALYSIS_DIR/scanpy,\
out_dir=$ANALYSIS_DIR/cellbender,\
conda_dir=$CONDA_DIR \
$SCRIPT_DIR/cellbender.sh


#### Use cellranger to find barcodes that passed cellranger QC on default settings
# Index genome
sbatch \
--export=working_dir=$WORKING_DIR,\
data_dir=$DATA_DIR,\
ref_fp=$REF_FP,\
gtf_fp=$GTF_FP \
$SCRIPT_DIR/cellranger_genome_index.sh


# 10x need the input fastq in a specific format 
# copy data folder
cp -r data/fastq/ data/fastq_10x
# rename files to cellranger compatible file name

# Align
sbatch \
--array=1-7 \
--export=working_dir=$WORKING_DIR,\
fastq_dir=$DATA_DIR/fastq_10x,\
in_filelist=filelist/snRNAseq_filelist.txt,\
out_dir=$ANALYSIS_DIR/cellranger,\
gtf_fp=$GTF_FP,\
conda_dir=$CONDA_DIR \
$SCRIPT_DIR/cellranger_count.sh



#### Spatial transcriptomics
# Upload spatial data into a folder data/spatial
# md5sum
ls data/spatial/*/*/*.fastq.gz > data/spatial_files.txt

# Index genome
sbatch \
--export=working_dir=$WORKING_DIR,\
data_dir=$DATA_DIR,\
genome_fp=$REF_FP,\
gtf_fp=$GTF_FP,\
conda_dir=$CONDA_DIR \
$SCRIPT_DIR/spaceranger_index.sh


# Align
SPATIAL_METADATA="data/spatial/metadata.txt"
sbatch \
--array=1-4 \
--export=working_dir=$WORKING_DIR,\
fastq_dir=$DATA_DIR/fastq,\
out_dir=$ANALYSIS_DIR/spaceranger,\
in_filelist=filelist/spatial_metadata.txt \
$SCRIPT_DIR/spaceranger.sh

