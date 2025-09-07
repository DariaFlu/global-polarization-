# Global Polarization of Λ Hyperons at the MPD Experiment

The afterburner for hyperon global polarization setting and analysis. Utilizes UrQMD modeling data.

## Description

This project consists a set of scripts:

1. `void add_enhanced_lambda(TString inputFile, TString outputFile, TString confInFile, Int_t enhanceStat)`

   - Script for statistical lambda enhancement.
   - Create XYZ polarization vector for all paricles in event. Only lambdas have non-null vector.
   - Output files contains **event** in UEvent format and **XYZ-vector** for particles polarization. 

2. `void simulate_lambda_decays(TString inputFile, TString outputFile, TString confInFile, Int_t enhanceStat)`

   - Same as `add_enhanced_lambda()` with additional functions.
   - Simulate Λ decay into proton and pion
   - Macro produces secondary particles and then puts polarization in proton

3. `void calc_global_polarization(TString fileIn , Int_t NFiles, TString OutFileName, Int_t enhancedFlag = -100)`
   - Measures global polarization with macro that provides transverse momentum and rapidity binning
   - Fits resulting angular distribution to get the value of global polarization
   - Global polarization is measured with proton polarization

4. `rootlogon.C` loads needed librarirs. 

5. Repository contais some additional .cpp and .sh scripts.

6. Examples for slurm script are presented in `startEnhLambdaPol.sh`, `startLambdaPolAnal.sh` and `startLambdaPolCalc.sh`,
---

## Requirements

The afterburner is intended to work with UrQMD generated data in unigen format.

## Usage

change path to mpdroot in CMakeList.txt (line 117, 120, 147, 155)
``` bash
source /path-to-mpdroot/thisroot.sh
mkdir build && cd build
cmake .. && make
cd ../
root .x rootlogon.C
add_enhanced_lambda(...) #or other script
````

