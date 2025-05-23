##DESCRIPTION: This is a loop to concatenate files from sequencing experiments.
##USAGE: ./MakeFileConcatenate_rbg.sh
##INSTRUCTION: User must define variables up to end of user input. Go into folder containing with script to run. Make script executeable before running.

#Create variable for list of files to be processed
export email="child@uchc.edu"
export threads=16
export memory="128G"
export bamdir=./
export samplistPATH=./
export bash_location=/home/FCAM/kchild/.bashrc
export outPATH=./
export slurmDir=./
export GENOMEDIR=/home/FCAM/kchild/GENOME/Homo_sapiens/UCSC/hg38/Sequence/STARIndex2.7.9a
export GTF_LOCATION=/home/FCAM/kchild/DATA/RNASEQ/Human_embryonic_heart_and_H9/Gencode_Human_annotations
export rMATS=/home/FCAM/kchild/ANALYSIS/RNASEQ/Human_embryonic_heart_and_H9/rMATS
export abbrv="rMATS_prep"

##################
##END USER INPUT##
##################

cat ${samplistPATH}/sample_list.txt | while read i
do
echo "#! /usr/bin/bash
#SBATCH --job-name=fileConcat_${i}_${abbrv}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c ${threads}
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mem=${memory}
#SBATCH --mail-type=END
#SBATCH -o sbatch_%x_%j.out
#SBATCH -e sbatch_%x_%j.err
#List hostname & date for troubleshooting
hostname
date
source ~/.bashrc
conda activate rmats_4.1.1
python ~/.conda/envs/rmats/rMATS/rmats.py \
--b1 ${i} \
--gtf ${GTF_LOCATION}/gencode.v25.annotation.gtf --bi ${GENOMEDIR} \
--libType fr-firststrand -t paired --nthread 16 \
--readLength 75 --variable-read-length --novelSS --od post --tmp prep --task prep
conda deactivate" > ${i}_${abbrv}.slurm
done
