#pragma once


#include "headers.hpp"

// Main function used by the generator.
//  - Lambda_yield : file with 2D pT-y yield histograms for lambda hyperon
//  - fBVal     : impact parameter
//  - ULambda      : output UParticle with updated 4-momentum
//  - rand      : TRandom3 instance for random sampling
void set_lambda_parameterization(TFile* Lambda_yield,
                                    Double_t fBVal, 
                                    UParticle &ULambda, 
                                    TRandom3 &rand); 