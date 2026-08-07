#!/bin/bash
# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Create CellBender-corrected matrices for AAV dnAP-1 snRNA-seq libraries.
# Steps: 1) Locate raw matrices; 2) run remove-background with configured settings; 3) write corrected matrices for R processing.
# =======================================
# Source reformatted from dnAP-1_snRNA/CellBender_create.sh.
# Set the required variables in config/.env before running this script.

###############################################################################
# SETUP: create one separate CellBender PBS job per OCM sample
###############################################################################

# Run:
#   bash make_cellbender_ocm_individual_pbs.sh
#
# Then submit all jobs:
#   bash $DNAP1_SNRNA_ROOT/data/CellBender_OCM/scripts/submit_all_cellbender_ocm_jobs.sh
#
# Or submit one sample manually:
#   qsub $DNAP1_SNRNA_ROOT/data/CellBender_OCM/scripts/cellbender_001_POOL_SAMPLE.pbs

set -euo pipefail

###############################################################################
# USER SETTINGS
###############################################################################

BASE_DIR=${DNAP1_SNRNA_ROOT:?Set DNAP1_SNRNA_ROOT in your environment}/data/CellRanger
POOLS=("Pool_MS" "Pool_SW" "Pool_VW")

OUT_DIR=${DNAP1_SNRNA_ROOT:?Set DNAP1_SNRNA_ROOT in your environment}/data/CellBender_OCM
SCRIPT_DIR="${OUT_DIR}/scripts"
LOG_DIR="${OUT_DIR}/logs"

PROJECT="zk16"
QUEUE="gpuvolta"
STORAGE="gdata/zk16+scratch/zk16"

WALLTIME="24:00:00"
NCPUS="12"
MEM="96GB"
NGPUS="1"
JOBFS="200GB"

# CellBender command.
# Edit this if your CellBender binary is somewhere specific.
# Example:
#   CELLBENDER_BIN=${CELLBENDER_BIN:?Set CELLBENDER_BIN in your environment}
CELLBENDER_BIN="${CELLBENDER_BIN:-cellbender}"

# Default CellBender parameters.
# Edit the manifest after generation if some samples need different values.
DEFAULT_EXPECTED_CELLS=8000
DEFAULT_TOTAL_DROPLETS_INCLUDED=25000
DEFAULT_FPR=0.01
DEFAULT_EPOCHS=150
DEFAULT_LOW_COUNT_THRESHOLD=100

mkdir -p "${OUT_DIR}" "${SCRIPT_DIR}" "${LOG_DIR}"

MANIFEST="${SCRIPT_DIR}/cellbender_ocm_manifest.tsv"
SUBMIT_ALL_SCRIPT="${SCRIPT_DIR}/submit_all_cellbender_ocm_jobs.sh"
SUBMIT_MISSING_SCRIPT="${SCRIPT_DIR}/submit_missing_cellbender_ocm_jobs.sh"
MISSING_SCRIPT="${SCRIPT_DIR}/find_missing_cellbender_outputs.sh"

###############################################################################
# HELPERS
###############################################################################

sanitize_name() {
    echo "$1" | sed 's/[^A-Za-z0-9_]/_/g' | sed 's/__*/_/g' | sed 's/^_//; s/_$//'
}

shell_quote() {
    printf "%q" "$1"
}

###############################################################################
# CREATE MANIFEST
###############################################################################

echo -e "idx\tpool\tsample\tinput_h5\toutput_h5\texpected_cells\ttotal_droplets_included\tfpr\tepochs\tlow_count_threshold\tpbs_script" > "${MANIFEST}"

idx=0

for pool in "${POOLS[@]}"; do
    SAMPLE_OUTS_DIR="${BASE_DIR}/${pool}/outs/per_sample_outs"

    if [[ ! -d "${SAMPLE_OUTS_DIR}" ]]; then
        echo "WARNING: Missing pool directory: ${SAMPLE_OUTS_DIR}" >&2
        continue
    fi

    while IFS= read -r sample_dir; do
        sample="$(basename "${sample_dir}")"

        candidates=(
            "${sample_dir}/sample_raw_feature_bc_matrix.h5"
            "${sample_dir}/count/sample_raw_feature_bc_matrix.h5"
            "${sample_dir}/raw_feature_bc_matrix.h5"
            "${sample_dir}/count/raw_feature_bc_matrix.h5"
            "${sample_dir}/sample_feature_bc_matrix.h5"
            "${sample_dir}/count/sample_feature_bc_matrix.h5"
        )

        input_h5=""

        for candidate in "${candidates[@]}"; do
            if [[ -f "${candidate}" ]]; then
                input_h5="${candidate}"
                break
            fi
        done

        if [[ -z "${input_h5}" ]]; then
            echo "WARNING: No raw/sample H5 matrix found for ${pool}/${sample}" >&2
            echo "         Checked:" >&2
            printf '         %s\n' "${candidates[@]}" >&2
            continue
        fi

        idx=$((idx + 1))

        sample_out_dir="${OUT_DIR}/${pool}/${sample}"
        mkdir -p "${sample_out_dir}"

        output_h5="${sample_out_dir}/${sample}_cellbender_filtered.h5"

        idx_pad="$(printf "%03d" "${idx}")"
        safe_pool="$(sanitize_name "${pool}")"
        safe_sample="$(sanitize_name "${sample}")"
        job_name="cb_${idx_pad}_${safe_pool}_${safe_sample}"

        # PBS job names can get ugly if very long.
        # Keep file name descriptive, but keep PBS -N shorter.
        short_job_name="$(echo "${job_name}" | cut -c1-45)"

        PBS_SCRIPT="${SCRIPT_DIR}/cellbender_${idx_pad}_${safe_pool}_${safe_sample}.pbs"

        echo -e "${idx}\t${pool}\t${sample}\t${input_h5}\t${output_h5}\t${DEFAULT_EXPECTED_CELLS}\t${DEFAULT_TOTAL_DROPLETS_INCLUDED}\t${DEFAULT_FPR}\t${DEFAULT_EPOCHS}\t${DEFAULT_LOW_COUNT_THRESHOLD}\t${PBS_SCRIPT}" >> "${MANIFEST}"

        POOL_Q="$(shell_quote "${pool}")"
        SAMPLE_Q="$(shell_quote "${sample}")"
        INPUT_H5_Q="$(shell_quote "${input_h5}")"
        OUTPUT_H5_Q="$(shell_quote "${output_h5}")"
        CELLBENDER_BIN_Q="$(shell_quote "${CELLBENDER_BIN}")"

        cat > "${PBS_SCRIPT}" <<EOF
#!/bin/bash
#PBS -P ${PROJECT}
#PBS -q ${QUEUE}
#PBS -l walltime=${WALLTIME}
#PBS -l ncpus=${NCPUS}
#PBS -l mem=${MEM}
#PBS -l ngpus=${NGPUS}
#PBS -l jobfs=${JOBFS}
#PBS -l storage=${STORAGE}
#PBS -N ${short_job_name}
#PBS -o ${LOG_DIR}/${job_name}.o
#PBS -e ${LOG_DIR}/${job_name}.e

set -eo pipefail

###############################################################################
# SAMPLE-SPECIFIC SETTINGS
###############################################################################

IDX="${idx}"
POOL=${POOL_Q}
SAMPLE=${SAMPLE_Q}
INPUT_H5=${INPUT_H5_Q}
OUTPUT_H5=${OUTPUT_H5_Q}

EXPECTED_CELLS="${DEFAULT_EXPECTED_CELLS}"
TOTAL_DROPLETS_INCLUDED="${DEFAULT_TOTAL_DROPLETS_INCLUDED}"
FPR="${DEFAULT_FPR}"
EPOCHS="${DEFAULT_EPOCHS}"
LOW_COUNT_THRESHOLD="${DEFAULT_LOW_COUNT_THRESHOLD}"

CELLBENDER_BIN=${CELLBENDER_BIN_Q}

###############################################################################
# START
###############################################################################

cd "\${PBS_O_WORKDIR}"

echo "============================================================"
echo "CellBender OCM individual PBS job"
echo "PBS job ID: \${PBS_JOBID}"
echo "Host: \$(hostname)"
echo "Start: \$(date)"
echo "============================================================"
echo "Index: \${IDX}"
echo "Pool: \${POOL}"
echo "Sample: \${SAMPLE}"
echo "Input: \${INPUT_H5}"
echo "Output: \${OUTPUT_H5}"
echo "Expected cells: \${EXPECTED_CELLS}"
echo "Total droplets included: \${TOTAL_DROPLETS_INCLUDED}"
echo "FPR: \${FPR}"
echo "Epochs: \${EPOCHS}"
echo "Low-count threshold: \${LOW_COUNT_THRESHOLD}"

###############################################################################
# LOAD CELLBENDER ENVIRONMENT
###############################################################################

# Edit this section if CellBender needs a conda/module environment.
#
# Example:
# module load python3
# source $MINIFORGE_ROOT/etc/profile.d/conda.sh
source activate CellBender

echo "Using CellBender command: \${CELLBENDER_BIN}"
which "\${CELLBENDER_BIN}" || true

###############################################################################
# CHECK INPUTS
###############################################################################

if [[ ! -f "\${INPUT_H5}" ]]; then
    echo "ERROR: Input file does not exist: \${INPUT_H5}" >&2
    exit 2
fi

OUTPUT_DIR="\$(dirname "\${OUTPUT_H5}")"
mkdir -p "\${OUTPUT_DIR}"

if [[ -f "\${OUTPUT_H5}" ]]; then
    echo "Output already exists, skipping:"
    echo "\${OUTPUT_H5}"
    echo "End: \$(date)"
    echo "============================================================"
    exit 0
fi



cd "\${OUTPUT_DIR}"
rm -f ckpt.tar.gz.tmp

"${CELLBENDER_BIN}" remove-background \
    --input "\${INPUT_H5}" \
    --output "\${OUTPUT_H5}" \
    --cuda \
    --expected-cells "\${EXPECTED_CELLS}" \
    --total-droplets-included "\${TOTAL_DROPLETS_INCLUDED}" \
    --fpr "\${FPR}" \
    --epochs "\${EPOCHS}" \
    --low-count-threshold "\${LOW_COUNT_THRESHOLD}"

###############################################################################
# CHECK OUTPUT
###############################################################################

if [[ ! -f "\${OUTPUT_H5}" ]]; then
    echo "ERROR: CellBender finished but output file was not created: \${OUTPUT_H5}" >&2
    exit 4
fi

echo "Finished sample: \${SAMPLE}"
echo "Output: \${OUTPUT_H5}"
echo "End: \$(date)"
echo "============================================================"
EOF

        chmod +x "${PBS_SCRIPT}"

    done < <(find "${SAMPLE_OUTS_DIR}" -mindepth 1 -maxdepth 1 -type d | sort)
done

N_JOBS="$(($(wc -l < "${MANIFEST}") - 1))"

if [[ "${N_JOBS}" -lt 1 ]]; then
    echo "ERROR: No samples found. Manifest contains no jobs: ${MANIFEST}" >&2
    exit 2
fi


echo ""
echo "DONE"
echo "Created one PBS script per sample."
echo ""
echo "Manifest:"
echo "  ${MANIFEST}"
echo ""
echo "Number of CellBender jobs:"
echo "  ${N_JOBS}"
echo ""
echo "PBS scripts:"
echo "  ${SCRIPT_DIR}/cellbender_*.pbs"
echo ""
echo "Submit all jobs:"
echo "  bash ${SUBMIT_ALL_SCRIPT}"
echo ""
echo "Check missing outputs:"
echo "  bash ${MISSING_SCRIPT}"
echo ""
echo "Submit only missing outputs:"
echo "  bash ${SUBMIT_MISSING_SCRIPT}"
echo ""
echo "Before submitting, check/edit expected cell numbers here:"
echo "  ${MANIFEST}"
