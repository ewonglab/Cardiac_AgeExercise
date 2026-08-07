#!/bin/bash
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Run HOMER known-motif enrichment for 10-week c-FOS peak sets.
# Steps: 1) Locate cleaned peak BED files; 2) run HOMER against mm10 for each set; 3) record job completion and logs.
# =======================================
# Source reformatted from 10_week_exercsie_Cut&Tag/40_homer_std_bg_10wk(1).pbs.
# Set CUTTAG_10_WEEK_ROOT and HOMER_BIN in config/.env before submission.
#PBS -P zk16
#PBS -q normalbw
#PBS -l walltime=08:00:00
#PBS -l mem=128GB
#PBS -l ncpus=28
#PBS -l storage=gdata/zk16+scratch/zk16
#PBS -N ct10wk_homer_std
#PBS -j oe
#PBS -l wd

set -euo pipefail
module purge
export PATH="${HOMER_BIN:?Set HOMER_BIN before submitting this job}:$PATH"

RUN="${CUTTAG_10_WEEK_ROOT:?Set CUTTAG_10_WEEK_ROOT before submitting this job}"
QC=$RUN/post_analysis_20260217_v1/qc_10wk
PEAKS=$QC/peaks_clean
OUT=$RUN/post_analysis_20260217_v1/motifs_10wk/homer_standard_genome_bg
mkdir -p "$OUT"
THREADS=${PBS_NCPUS:-28}

for peaks in $PEAKS/*C-FOS.clean.bed; do
  [ -f "$peaks" ] || continue
  group=$(basename "$peaks" .clean.bed)
  outdir=$OUT/$group
  mkdir -p "$outdir"
  findMotifsGenome.pl "$peaks" mm10 "$outdir" -size given -nomotif -p "$THREADS" > "$outdir/homer.stdout.log" 2> "$outdir/homer.stderr.log"
done

echo DONE > "$OUT/DONE.flag"
