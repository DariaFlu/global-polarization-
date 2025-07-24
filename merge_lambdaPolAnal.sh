#!/bin/bash
#SBATCH -J polLambdaMerge
#SBATCH --output=/lhep/users/dflusova/lambda/afterburner/v.6/logs/merge_%j.log
#SBATCH --error=/lhep/users/dflusova/lambda/afterburner/v.6/logs/merge_%j.err
#SBATCH --requeue
#SBATCH --time=02:00:00

# Load necessary environment
source /cvmfs/nica.jinr.ru/sw/os/login.sh latest
module add mpddev
source /lhep/users/dflusova/mpdroot/install/config/env.sh

# Set paths
INPUT_DIR="/lhep/users/dflusova/lambda/afterburner/v.6/out"
OUTPUT_DIR="/lhep/users/dflusova/lambda/afterburner/v.6/lambda_analysis"
LOG_DIR="/lhep/users/dflusova/lambda/afterburner/v.6/logs"
OUTPUT_FILE="${OUTPUT_DIR}/merged_urqmd_xexe_2.87gev_$(date +%Y%m%d_%H%M).root"

# Create directories
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${INPUT_DIR}/processed"

# Function to verify ROOT files
verify_root_file() {
    local file=$1
    root -l -b -q <<EOF
{
    TFile* f = TFile::Open("${file}");
    if (!f || f->IsZombie()) {
        std::cerr << "Invalid ROOT file: ${file}" << std::endl;
        exit(1);
    }
    if (!f->GetListOfKeys()->Contains("decayes")) {
        std::cerr << "Missing 'decayes' tree in: ${file}" << std::endl;
        exit(1);
    }
    f->Close();
    exit(0);
}
EOF
}

# Create and verify file list
echo "Creating verified file list..."
rm -f "${LOG_DIR}/filelist.txt" "${LOG_DIR}/missing_files.log"
for i in {1..50}; do
    file="${INPUT_DIR}/result_urqmd_xexe_2.87gev_mf_6195240_${i}.mcini.root"
    if [[ -f "${file}" ]]; then
        if verify_root_file "${file}"; then
            echo "${file}" >> "${LOG_DIR}/filelist.txt"
        else
            echo "Corrupted file: ${file}" >> "${LOG_DIR}/missing_files.log"
        fi
    else
        echo "Missing file: ${file}" >> "${LOG_DIR}/missing_files.log"
    fi
done

# Check if we have files to merge
if [[ ! -s "${LOG_DIR}/filelist.txt" ]]; then
    echo "Error: No valid input files found!" >&2
    exit 1
fi

# Merge using hadd (faster than TChain for simple merges)
echo "Merging files using hadd..."
hadd -f "${OUTPUT_FILE}" $(cat "${LOG_DIR}/filelist.txt") 2>&1 | tee "${LOG_DIR}/merge_progress.log"

# Verify merged output
verify_root_file "${OUTPUT_FILE}"

# Count entries in merged file
ENTRIES=$(root -l -b -q <<EOF
{
    TFile* f = TFile::Open("${OUTPUT_FILE}");
    TTree* t = (TTree*)f->Get("decayes");
    if (!t) {
        std::cerr << "Error reading merged tree" << std::endl;
        exit(1);
    }
    std::cout << t->GetEntries();
    f->Close();
    exit(0);
}
EOF
)

echo "Merged successfully. Total entries: ${ENTRIES}"

# Archive processed files
while read -r file; do
    mv "${file}" "${INPUT_DIR}/processed/"
done < "${LOG_DIR}/filelist.txt"

# Generate checksum
md5sum "${OUTPUT_FILE}" > "${OUTPUT_FILE}.md5"

# Final status
echo "Merge completed at $(date)"
echo "Output file: ${OUTPUT_FILE}"
echo "Size: $(du -h "${OUTPUT_FILE}" | cut -f1)"
echo "MD5: $(cat "${OUTPUT_FILE}.md5")"
