# Salmon Wound Atlas
This repository contains the code for analysing single-nuclei sequencing (snRNA-Seq) and spatial transcriptomics (ST) data on the wound healing of Atlantic salmon. The processed data can also be visualised on our webapp at https://bshihlab.shinyapps.io/SalmonWoundAtlas/.
## Data availability
Please first download raw data from NCBI Sequence Read Archive (SRA) or NCBI Gene Expression Omnibus (GEO). 

SRA (snRNA-Seq): SRR29294027, SRR24471130, SRR38814853, SRR38814854, SRR38814855, SRR38814856.

SRA (ST): SRR27965648, SRR27965649, SRR30953451, and SRR30953452

GEO (snRNA-Seq): GSE33584

## Processing pipelines
Please use the script/_full_process.sh for processing the raw sequencing files. The parameters may need to be changed in the __full_process_settings.sh file to match your local settings. The scripts were written for running on a slurm high-performance computing cluster.
In brief, the reads were aligned to Salmo_salar.Ssal_v3.1 using STAR, and ambient RNA contamination were removed using cellbender. Cellranger was also used for alignment, but only for the purpose for identifying cells that passes its QC checks for filtering in the downstream QC filtering in R.
After Processing steps associated with R are detailed in script/scRNAseq.r, whereas those for ST are detailed in script/spatial.r.
Processed data generated from R are used to make the shinyapp.
