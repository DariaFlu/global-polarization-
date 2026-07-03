// post_phi_flow.h
// MC flow post-analyzer for phi(1020) in UniGen (UEvent/UParticle).
//
// Measures vn w.r.t. TRUE reaction plane PsiRP = event->GetPhi():
//   v_n(bin) = < cos(n*(phi - PsiRP)) > in bins of (cent, y, pT)
//
// Run:
//   root -l -b -q 'post_phi_flow.C+("out_files*.root","phi_flow_QA.root",0)'
//   root -l -b -q 'post_phi_flow.C+("out_files*.root","phi_flow_QA_enh.root",1)'

#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include "UEvent.h"
#include "UParticle.h"

Double_t get_centrality(Double_t fBVal);

static int centClass(Double_t cent);
static const char* centName(int ic);

static inline double wrapToPi(double x);

int post_phi_flow(const char* inPattern = "out_files*.root",
                  const char* outFile   = "phi_flow_QA.root",
                  int enhancedFactor    = 10); 