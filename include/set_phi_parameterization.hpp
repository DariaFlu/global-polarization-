#pragma once

#include "TF1.h"
#include "TF2.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include "headers.hpp"

// Main function used by the generator.
//  - Phi_yield : file with 2D pT-y yield histograms for phi mesons
//  - fBVal     : impact parameter
//  - UPhi      : output UParticle with updated 4-momentum
//  - rand      : TRandom3 instance for random sampling
void set_phi_parameterization(TFile *Phi_yield, 
                                Double_t fBVal, 
                                UParticle &UPhi, 
                                TRandom3 &rand);
