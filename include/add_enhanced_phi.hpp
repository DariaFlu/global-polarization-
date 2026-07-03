#pragma once

#include "headers.hpp"
#include "set_phi_parameterization.hpp"
#include "set_lambda_parameterization.hpp"

int add_enhanced_phi(const TString& outputFile,
                    Long64_t fileStart,
                    Long64_t fileCount,
                    const TString& pattern,
                    const TString& qaFile,
                    Int_t enhanceStat,
                    UInt_t seed);