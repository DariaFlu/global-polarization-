#pragma once


// ROOT Core includes
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TTree.h"
#include "TBranch.h"
#include "TStyle.h"
#include "TSystem.h"

#include "TObjArray.h"

// Physics-specific includes
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TParticle.h"
#include <Math/RotationZ.h>
#include <Math/Vector3D.h>

// Random number generation
#include "TRandom3.h"

// Histogramming
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"

#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"


// Utility includes
#include "TMath.h"
#include "TString.h"

// MPDRoot objects
#include "URun.h"
#include "UEvent.h"
#include "UParticle.h"

// Standard C++ includes
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <map>

#include "getters.hpp"

class TString;
class TClonesArray;
class UParticle;