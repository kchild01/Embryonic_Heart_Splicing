#!/bin/bash
#SBATCH -J rMATS_runCS13to23FQ
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=child@uchc.edu
#SBATCH --mem=128G
#SBATCH -o rMATS_runCS13to23FQ_%j.out
#SBATCH -e rMATS_runCS13to23FQ_%j.err
export PATH="/home/FCAM/jcotney/miniconda3/bin:$PATH"
cd ~/.conda/envs
conda init bash
source ~/.bashrc
conda activate rmats_4.1.1
export OUTPUT=/home/FCAM/kchild/ANALYSIS/RNASEQ/Human_embryonic_heart_and_H9/rMATS/output_rMATSFQ_GAv25_firststrand/13to23
export GTF_LOCATION=/home/FCAM/kchild/DATA/RNASEQ/Human_embryonic_heart_and_H9/Gencode_Human_annotations
export GENOMEDIR=/home/FCAM/kchild/GENOME/Homo_sapiens/UCSC/hg38/Sequence/STARIndex2.7.9a
export rMATS=/home/FCAM/kchild/ANALYSIS/RNASEQ/Human_embryonic_heart_and_H9/rMATS
cd /home/FCAM/kchild/DATA/RNASEQ/Human_embryonic_heart_and_H9/STAR_aligned_reads/rMATS
python ~/.conda/envs/rmats/rMATS/rmats.py \
--s1 $rMATS/s3_1.txt \
--gtf $GTF_LOCATION/gencode.v25.annotation.gtf --bi $GENOMEDIR \
--libType fr-firststrand -t paired --nthread 16 \
--readLength 75 --variable-read-length --novelSS --od $OUTPUT --tmp $OUTPUT/tmp --task prep
conda deactivate
