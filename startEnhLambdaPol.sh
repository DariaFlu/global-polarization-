#!/bin/bash
#SBATCH -D "/scratch3/dflusova/afterburner/slurm"
#SBATCH -J polLambdaAnal
#SBATCH -p nica
#SBATCH -a 1-100
#SBATCH --requeue
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --exclude=ncx112,ncx115,ncx117,ncx121,ncx127,ncx158,ncx159,ncx171,ncx181,ncx207,ncx214,ncx216,ncx222,ncx223,ncx224,ncx225,ncx227

# Load necessary environment
source /cvmfs/nica.jinr.ru/sw/os/login.sh latest
module add mpddev
source /lhep/users/dflusova/mpdroot/install/config/env.sh

# Set job identifiers
export JOB_ID=${SLURM_ARRAY_JOB_ID}
export TASK_ID=${SLURM_ARRAY_TASK_ID}

# Define paths (modify these as needed)
INPUT_DIR="/eos/nica/mpd/users/parfenov/SimData/UrQMD/xexe_2.87gev_mf/6195240/files/mcini/" # dir with input UrQMD data
CONFIG_DIR="/lhep/users/dflusova/lambda/afterburner/v.7/global-polarization-/" #path to afterburner
SRC_DIR="/lhep/users/dflusova/lambda/afterburner/global-polarization-/src/" # path to afterburner/src
OUTPUT_DIR="/scratch3/dflusova/afterburner/v.7/global-polarization-/out/" # path to afterburner/out
BUILD_DIR="/lhep/users/dflusova/lambda/afterburner/v.7/global-polarization-/build" # path to afterburner/build
INCLUDE_DIR="/lhep/users/dflusova/lambda/afterburner/v.7/global-polarization-/include" # path to afterburner/include

# Ensure dirs exist
mkdir -p "${OUTPUT_DIR}"

# Wait for CVMFS and EOS to be available
for i in {1..5}; do
    ls /cvmfs/ >/dev/null 2>&1 && ls /eos/nica/ >/dev/null 2>&1 && break
    sleep 15
done

# Build file names
INPUT_FILE="urqmd_xexe_2.87gev_mf_6195240_${TASK_ID}.mcini.root"
CONFIG_FILE="qa_out_xexe.root"
OUTPUT_FILE="result_urqmd_xexe_2.87gev_mf_6195240_${TASK_ID}.mcini.root"
# POLARIZATION_OUTPUT="result_global_polarization_urqmd_xexe_2.87gev_mf_6195240_${TASK_ID}.mcini.root"

#cd "${CONFIG_DIR}"

# Run the analysis with timing measurement
echo "Starting analysis for task ${TASK_ID} at $(date)"
START_TIME=$(date +%s)

# Make compiled lib visible (so gSystem->Load can find it)
export LD_LIBRARY_PATH="${BUILD_DIR}:${LD_LIBRARY_PATH}"

# ----- Run ROOT non-interactively (use compiled plugin; no JIT, no AutoDict_*) -----

echo "ROOT start"

root -b <<EOF
.x ${CONFIG_DIR}rootlogon.C
add_enhanced_lambda(
  "${INPUT_DIR}${INPUT_FILE}",
  "${OUTPUT_DIR}${OUTPUT_FILE}",
  "${CONFIG_DIR}${CONFIG_FILE}",
  2
)
.q
EOF

echo "ROOT stop"

END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "Analysis for task ${TASK_ID} completed at $(date)"
echo "Time taken: $((ELAPSED_TIME / 3600)) hours, $(( (ELAPSED_TIME % 3600) / 60 )) minutes, $((ELAPSED_TIME % 60)) seconds"

# ----- Verify output -----
if [[ -f "${OUTPUT_DIR}${OUTPUT_FILE}" ]]; then
  echo "Output OK: ${OUTPUT_DIR}${OUTPUT_FILE}"
else
  echo "Error: Output file was not created!"
  exit 1
fi

echo "Job ${TASK_ID} completed successfully"