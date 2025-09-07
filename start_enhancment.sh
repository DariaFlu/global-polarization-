#!/bin/bash

# Load necessary environment
source /cvmfs/nica.jinr.ru/sw/os/login.sh latest
module add mpddev
source /lhep/users/dflusova/mpdroot/install/config/env.sh

echo "=== SETTING ENVIRONMENT ==="

# Define the directory name
DIR_NAME="out"

# Check if directory exists and clear it
if [ -d "$DIR_NAME" ]; then
    echo "Clearing existing directory: $DIR_NAME"
    # Remove all contents of the directory
    rm -rf "${DIR_NAME:?}/"*
else
    echo "Creating directory: $DIR_NAME"
    # Create the directory
    mkdir -p "$DIR_NAME"
fi

echo "Operation completed - '$DIR_NAME' directory is ready"

echo "=== ENHANCMENT STARTS ==="


echo "=== ENHANCMENT DONE ==="

