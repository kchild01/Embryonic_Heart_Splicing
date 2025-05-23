#!/bin/bash
#SBATCH -J rmats2sashimiplot_run
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=child@uchc.edu
#SBATCH --mem=64G
#SBATCH -o rmats2sashimiplot_runSE_%j.out
#SBATCH -e rmats2sashimiplot_runSE_%j.err
source ~/.bashrc
conda activate rmats_4.1.1
export OUTPUT=/home/FCAM/kchild/ANALYSIS/RNASEQ/Human_embryonic_heart_and_H9/rMATS/output_rMATSFQ_GAv25_firststrand/13to23/sashimi/SE/13v23
export rMATS=/home/FCAM/kchild/ANALYSIS/RNASEQ/Human_embryonic_heart_and_H9/rMATS/output_rMATSFQ_GAv25_firststrand/13to23/sashimi/SE
export BAM=/home/FCAM/kchild/ANALYSIS/RNASEQ/Human_embryonic_heart_and_H9/rMATS/output_rMATSFQ_GAv25_firststrand/13to23/tmp
cd /home/FCAM/kchild/ANALYSIS/RNASEQ/Human_embryonic_heart_and_H9/rMATS
rmats2sashimiplot \
--b1 $BAM/bam1_1/Aligned.sortedByCoord.out.bam,$BAM/bam1_2/Aligned.sortedByCoord.out.bam,$BAM/bam1_3/Aligned.sortedByCoord.out.bam \
--b2 $BAM/bam1_22/Aligned.sortedByCoord.out.bam,$BAM/bam1_23/Aligned.sortedByCoord.out.bam,$BAM/bam1_24/Aligned.sortedByCoord.out.bam \
-t SE -e SE.MATS.JC.txt \
--l1 CS13 --l2 CS23 --exon_s 1 --intron_s 5 \
-o $OUTPUT --group-info $rMATS/grouping.gf
conda deactivate