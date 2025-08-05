#!/bin/bash
#SBATCH -D "/scratch3/dflusova/afterburner/slurm"
#SBATCH -J polLambdaAnal
#SBATCH -p nica
#SBATCH -a 1-2000
#SBATCH --requeue
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --exclude=ncx112,ncx115,ncx117,ncx121,ncx147,ncx153,ncx156,ncx158,ncx159,ncx171,ncx181,ncx207,ncx214,ncx216,ncx222,ncx223,ncx224,ncx225,ncx227

# Load necessary environment
source /cvmfs/nica.jinr.ru/sw/os/login.sh latest
module add mpddev
source /lhep/users/dflusova/mpdroot/install/config/env.sh

# Set job identifiers
export JOB_ID=${SLURM_ARRAY_JOB_ID}
export TASK_ID=${SLURM_ARRAY_TASK_ID}

# Define paths (modify these as needed)
INPUT_DIR="/eos/nica/mpd/users/parfenov/SimData/UrQMD/xexe_2.87gev_mf/6195240/files/mcini/"
CONFIG_DIR="/lhep/users/dflusova/lambda/afterburner/v.6/"
# OUTPUT_DIR="/lhep/users/dflusova/lambda/afterburner/v.6/out/"
OUTPUT_DIR="/scratch3/dflusova/afterburner/out/"

# Wait for CVMFS and EOS to be available
for i in {1..5}; do
    ls /cvmfs/ >/dev/null 2>&1 && ls /eos/nica/ >/dev/null 2>&1 && break
    sleep 15
done

# Create unique working directory
WORK_DIR="${CONFIG_DIR}lambda_analysis"
mkdir -p "${WORK_DIR}" || { echo "Failed to create working directory"; exit 1; }
cd "${WORK_DIR}" || { echo "Failed to enter working directory"; exit 1; }

# Build file names
INPUT_FILE="urqmd_xexe_2.87gev_mf_6195240_${TASK_ID}.mcini.root"
CONFIG_FILE="qa_out_xexe.root"
OUTPUT_FILE="result_urqmd_xexe_2.87gev_mf_6195240_${TASK_ID}.mcini.root"
# POLARIZATION_OUTPUT="result_global_polarization_urqmd_xexe_2.87gev_mf_6195240_${TASK_ID}.mcini.root"

#cd "${CONFIG_DIR}"

# Run the analysis with timing measurement
echo "Starting analysis for task ${TASK_ID} at $(date)"
START_TIME=$(date +%s)


# Run the analysis
root -l -b <<EOF
gSystem->Load("${CONFIG_DIR}AutoDict_vector_TVector3__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_vector_UParticle__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_vector_TVector3__cxx.so")
gSystem->Load("${CONFIG_DIR}AutoDict_vector_UParticle__cxx.so")
.L ${CONFIG_DIR}read_unigen_root.cpp
simulate_lambda_decays("${INPUT_DIR}${INPUT_FILE}", "${OUTPUT_DIR}${OUTPUT_FILE}", "${CONFIG_DIR}${CONFIG_FILE}", 50)
.q
EOF

END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "Analysis for task ${TASK_ID} completed at $(date)"
echo "Time taken: $((ELAPSED_TIME / 3600)) hours, $(( (ELAPSED_TIME % 3600) / 60 )) minutes, $((ELAPSED_TIME % 60)) seconds"

# Verify and copy results
if [[ -f "${OUTPUT_FILE}" ]]; then
    cp "${OUTPUT_FILE}" "${OUTPUT_DIR}"/
    # cp "${POLARIZATION_OUTPUT}" "${OUTPUT_DIR}"/
else
    echo "Error: Output file was not created!"
    exit 1
fi

# Cleanup
cd "${TMP}" && rm -rf "${WORK_DIR}"

echo "Job ${TASK_ID} completed successfully"