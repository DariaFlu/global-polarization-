# Load necessary environment
source /cvmfs/nica.jinr.ru/sw/os/login.sh latest
module add mpddev
source /lhep/users/dflusova/mpdroot/install/config/env.sh


# Define paths (modify these as needed)
BUILD_DIR="/lhep/users/dflusova/lambda/afterburner/release/build" # path to afterburner/build
INCLUDE_DIR="/lhep/users/dflusova/lambda/afterburner/release/include" # path to afterburner/include
INPUT_DIR="/scratch3/dflusova/afterburner/out/" # path to afterburner/out
OUTPUT_DIR="/lhep/users/dflusova/lambda/afterburner/release/out" # path to afterburner/out

# Build file names
INPUT_FILE="result_urqmd_xexe_2.87gev_mf_6195240"
OUTPUT_FILE="result_polarization_urqmd_xexe_2.87gev_mf_6195240.mcini.root"

# POLARIZATION_OUTPUT="result_global_polarization_urqmd_xexe_2.87gev_mf_6195240_${TASK_ID}.mcini.root"

#cd "${CONFIG_DIR}"

# Run the analysis with timing measurement
echo "Starting analysis at $(date)"
START_TIME=$(date +%s)

# Make compiled lib visible (so gSystem->Load can find it)
export LD_LIBRARY_PATH="${BUILD_DIR}:${LD_LIBRARY_PATH}"


echo "ROOT start"

root -b <<EOF
#include "${INCLUDE_DIR}/calc_global_polarization.hpp"
gSystem->Load("${BUILD_DIR}/libgp_macros");
gSystem->Load("${BUILD_DIR}/libgp_dict");
calc_global_polarization(
                        ${INPUT_DIR},
                        ${INPUT_FILE},
                        100,
                        ${OUTPUT_DIR},
                        ${OUTPUT_FILE},
                        2
                        )
.q
EOF

echo "ROOT stop"

END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "Analysis completed at $(date)"
echo "Time taken: $((ELAPSED_TIME / 3600)) hours, $(( (ELAPSED_TIME % 3600) / 60 )) minutes, $((ELAPSED_TIME % 60)) seconds"


echo "Completed successfully"