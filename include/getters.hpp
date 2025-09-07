#pragma once

#include "headers.hpp"

// #include "TRandom3.h"
// #include "TMath.h"
// #include "TF1.h"

// #include <Math/RotationZ.h>
// #include <Math/Vector3D.h>

Double_t get_mean_polarization(Double_t sNN, Double_t centrality);// Value in % 
Double_t get_random_value(Double_t fMean, Double_t fSigma);
Double_t get_positive_phi(const Double_t& phi);
Double_t get_centrality  (Double_t fBVal);
Double_t get_costh(Double_t alpha, Double_t pol = 0.6);
Double_t get_V2   (Double_t sNN, Double_t centrality, Double_t lambda_pT, Double_t lambda_y);
Int_t    get_number_of_bin(Double_t fValue, Double_t fMinValue, Double_t fMaxValue, Int_t NBins); 
ROOT::Math::XYZVector get_pol_lambda(UParticle& lambda, Double_t _fpoly, Double_t _fSigmaPol = 0.3);