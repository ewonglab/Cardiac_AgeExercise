#!/bin/bash
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Run CellBender remove-background for primary snRNA-seq matrices.
# Steps: 1) Locate raw Cell Ranger outputs; 2) invoke configured CellBender settings; 3) write corrected matrices for QC.
# =======================================
# Source reformatted from snRNA_AgeEx/RunCellbedner.sh.
# Set the required variables in config/.env before running this script.

#PBS -P zk16
#PBS -N cellbender
#PBS -l ngpus=1
#PBS -l ncpus=12
#PBS -l mem=64GB
#PBS -l walltime=48:00:00
#PBS -q gpuvolta
#PBS -l wd
#PBS -l storage=scratch/zk16+gdata/zk16+gdata/gd7
#PBS -j oe


conda activate CellBender
Runs=("MAE2_1" "MAE4_3" "MAS1_4" "MAS3_2" "YE2_1" "YE4_3" "YS1_4" "YS3_2")
for Run in "${Runs[@]}"
do
	echo $Run
    fname=${HEART_SNRNA_ROOT:?Set HEART_SNRNA_ROOT in your environment}/data/CellRanger/${Run}/outs/raw_feature_bc_matrix.h5
    Out=${HEART_SNRNA_ROOT:?Set HEART_SNRNA_ROOT in your environment}/data/2022-01-23.Cellbender/${Run}_cellbender.h5
    echo $fname
    echo $Out
    cellbender remove-background \
                 --input ${fname} \
                 --output ${Out} \
                 --cuda \
                 --expected-cells 20000 \
                 --total-droplets-included 30000 \
                 --fpr 0.01 \
                 --epochs 150
done

