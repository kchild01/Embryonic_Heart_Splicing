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
source ~/.bashrc
conda activate rmats_4.1.1
export GTF_LOCATION=/home/FCAM/kchild/DATA/RNASEQ/Human_embryonic_heart_and_H9/Gencode_Human_annotations
export GENOMEDIR=/home/FCAM/kchild/GENOME/Homo_sapiens/UCSC/hg38/Sequence/STARIndex2.7.9a
python ~/.conda/envs/rmats/rMATS/rmats.py \
--b1 b1_24.txt \
--gtf $GTF_LOCATION/gencode.v25.annotation.gtf --bi $GENOMEDIR \
--libType fr-firststrand -t paired --nthread 16 \
--readLength 75 --variable-read-length --novelSS --od ./ --tmp ../prep --task post --statoff
conda deactivate
