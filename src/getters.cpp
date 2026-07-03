#include "getters.hpp"


Double_t get_random_value(Double_t fMean, Double_t fSigma)
{
    TRandom3* rand = new TRandom3(0);  // Seed with 0 for reproducibility
    Double_t randomVal = rand->Gaus(fMean, fSigma);
    delete rand;
    return randomVal;
}

Int_t get_number_of_bin(Double_t fValue, Double_t fMinValue, Double_t fMaxValue, Int_t NBins) { //Function return number of bin
    if (NBins <= 0) return -1;
    if (fValue < fMinValue) return -1;  // Underflow
    if (fValue >= fMaxValue) return -1; // Overflow
    
    Double_t binWidth = (fMaxValue - fMinValue)/NBins;
    //std::cout << "bin_value" <<  (fValue - fMinValue)/binWidth << std::endl;
    Int_t bin = static_cast<Int_t>((fValue - fMinValue)/binWidth);
    
    // Protect against edge case where fValue == fMaxValue
    return (bin < NBins) ? bin : NBins-1;
}


Double_t get_positive_phi(const Double_t& phi){
    Double_t phi_rot = phi;
    if (phi < 0) phi_rot += 2.*TMath::Pi();
    return phi_rot;
}


Double_t get_costh(Double_t alpha, Double_t pol) {

    // TF1* randomCosTheta = new TF1("randomCosTheta", "1+(TMath::Pi()*[0]*[1])/8*x", -1.0, 1.0);
    // TF1* randomCosTheta = new TF1("randomCosTheta", "1+[0]*[1]*x", -1.0, 1.0);
    TF1* randomCosTheta = new TF1("randomCosTheta", "(1+[0]*x)", -1.0, 1.0);
    randomCosTheta->SetNpx(10000);

    // Set the parameters
    // randomCosTheta->SetParameters(alpha);

    // randomCosTheta->SetParameters(alpha, pol);
    randomCosTheta->SetParameters(alpha);

    // Generate a random value from the distribution
    Double_t value = randomCosTheta->GetRandom();
    // std::cout<<"random costh = "<<value<<std::endl;
    // Clean up
    delete randomCosTheta;
    
    return value;

}

Double_t get_V2(Double_t sNN, Double_t centrality, Double_t lambda_pT, Double_t lambda_y) {

	// Double_t v2 = 1.05*(0.8132-11.4/TMath::Power(sNN,2.1)) 
    //                   *((-0.01105*centrality + 0.000162*TMath::Power(centrality,2) )/(-0.01105*25 + 0.000162*TMath::Power(25,2) )) 
    //                   *(-1) * ( (0.32172*TMath::Power(lambda_y,1.805)-0.1578)*lambda_pT + (0.05838*TMath::Power(lambda_y,3.2172) - 0.0179)
    //                   *TMath::Exp(lambda_pT*(1.4562*TMath::Power(lambda_y,3.28814) + 0.22912))
    //                   *TMath::Sin(lambda_pT*(-1.004*TMath::Power(lambda_y,4.6398) + 3.88097)+(-3.6075*TMath::Power(lambda_y,4.6582) + 3.35966))); //v2 snn-cent-pT-y func
    // if (sNN <= 0 || centrality < 0 || lambda_pT < 0 || lambda_y < 0) return 0.;  // or throw an error
    lambda_y = TMath::Abs(lambda_y);
    
    // =============================================
    // 1. Energy-dependent term (depends on sNN)
    // =============================================
    Double_t energy_term = 1.05 * (0.8132 - 11.4 / TMath::Power(sNN, 2.1));

    // =============================================
    // 2. Centrality-dependent term (normalized at centrality=25)
    // =============================================
    Double_t centrality_numerator = 
        -0.01105 * centrality + 0.000162 * TMath::Power(centrality, 2);
    Double_t centrality_denominator = 
        -0.01105 * 25 + 0.000162 * TMath::Power(25, 2);
    Double_t centrality_term = centrality_numerator / centrality_denominator;

    // =============================================
    // 3. Sign flip (-1 factor)
    // =============================================
    Double_t sign_flip = -1;

    // =============================================
    // 4. Rapidity (lambda_y) and pT (lambda_pT) dependent term
    // =============================================
    // 4a. Linear term in lambda_pT (depends on lambda_y)
    Double_t linear_pT_term = 
        (0.32172 * TMath::Power(lambda_y, 1.805) - 0.1578) * lambda_pT;

    // 4b. Nonlinear term (exponential + sinusoidal dependence)
    Double_t exp_argument = 
        lambda_pT * (1.4562 * TMath::Power(lambda_y, 3.28814) + 0.22912);
    Double_t sin_phase = 
        lambda_pT * (-1.004 * TMath::Power(lambda_y, 4.6398) + 3.88097) 
        + (-3.6075 * TMath::Power(lambda_y, 4.6582) + 3.35966);
    
    Double_t nonlinear_term = 
        (0.05838 * TMath::Power(lambda_y, 3.2172) - 0.0179) 
        * TMath::Exp(exp_argument) 
        * TMath::Sin(sin_phase);

    // Combine all terms into final v2 expression
    Double_t v2 = 
        energy_term 
        * centrality_term 
        * sign_flip 
        * (linear_pT_term + nonlinear_term);

    // std::cout << "ERROR: NaN in nonlinear_term. Args: "
    //     << "exp_argument=" << exp_argument << ", sin_phase=" << sin_phase <<" , centrality_denominator = " << centrality_denominator << std::endl;

    // if (TMath::Abs(centrality_denominator) < 1e-10) {
    //     return 0.;  // Avoid division by zero
    // }

    // if (TMath::IsNaN(nonlinear_term)) {
    //     std::cout << "ERROR: NaN in nonlinear_term. Args: "
    //                 << "exp_argument=" << exp_argument << ", sin_phase=" << sin_phase << std::endl;
    // }
    // if (!TMath::Finite(sin_phase)) {
    //     return 0.;
    // }

    return v2;

}


Double_t get_mean_polarization(Double_t sNN, Double_t centrality){ return (2.8569/ TMath::Power(sNN,0.955513) ) * (2.4702 - 0.0461*centrality + 0.0042 * TMath::Power(centrality, 2)); }// Value in % 

// =============================================================================
// Centrality mapping 
// =============================================================================
Double_t get_centrality(Double_t fBVal){
    if(fBVal<3.44) {return 5;}
    else if(fBVal<4.88) {return 15;}
    else if(fBVal<5.84) {return 25;}
    else if(fBVal<6.64) {return 35;}
    else if(fBVal<7.44) {return 45;}
    else if(fBVal<8.08) {return 55;}
    else if(fBVal<8.72) {return 65;}
    else if(fBVal<9.36) {return 75;}
    else if(fBVal<9.84) {return 85;}
    else {return 95;}
}


ROOT::Math::XYZVector get_pol_lambda(UParticle& lambda, Double_t _fpoly, Double_t _fSigmaPol){ 
    
    Double_t fPhi = 0;

    Double_t fpolx  = 0; 
    Double_t fpolSx = 0.07;  
    Double_t fpoly  = _fpoly;
    Double_t fpolSy = 0.07; 
    Double_t fpolz  = 0;
    Double_t fpolSz = 0.07;

    TRandom3* rand = new TRandom3(0);  // Seed with 0 for reproducibility
    Double_t polX = rand->Gaus(fpolx,fpolSx);
    Double_t polY = rand->Gaus(fpoly,fpolSy);
    Double_t polZ = rand->Gaus(fpolz,fpolSz);

    // Double_t polX = gRandom->Gaus(fpolx,fpolSx);  // generate polarization direction
    // Double_t polY = gRandom->Gaus(fpoly,fpolSy);  // generate polarization direction /100
    // Double_t polZ = gRandom->Gaus(fpolz,fpolSz);  // generate polarization direction
    ROOT::Math::XYZVector polarizationVec = ROOT::Math::XYZVector(polX, polY, polZ);

    // std::cout<<"TVector3(polX, polY, polZ) "<<polX<<" "<<polY<<" "<<polZ<<std::endl;

    ROOT::Math::RotationZ rotateRP(fPhi);  // set rotation transformation
    polarizationVec = rotateRP*polarizationVec;  // rotate reaction plane

    // Float_t polmag = TMath::Sqrt(polarizationVec.mag());
    Float_t polmag = TMath::Sqrt(polarizationVec.mag2());

    // special case, overpolarized -> 100% polarized
    if (polmag>1.) {
        polarizationVec *= 1./polmag;  // scale to 1
        polmag = 1.;
    }
    // special case, no polarization -> generate random direction
    if (polmag == 0) {
    // random unitary vector, see https://mathworld.wolfram.com/SpherePointPicking.html
        Double_t x1,x2,R2,R;
        do {
            // x1 = 1-2*gRandom->Rndm();
            // x2 = 1-2*gRandom->Rndm();
            x1 = 1-2*rand->Rndm();
            x2 = 1-2*rand->Rndm();
            R2 = x1*x1+x2*x2;
        } while(R2 >= 1);
        R=2*TMath::Sqrt(1-R2);
        polarizationVec.SetXYZ(x1*R, x2*R, 1-2*R2);  // random unitary vector - background for signal
    }

    // Float_t xxx = gRandom->Rndm();
    Float_t xxx = rand->Rndm();
    if (xxx < 1./2.*(1.-polmag)) {  // probability of spin flip
        polarizationVec *= -1.;  // spin flip according to mean polarization
    }

    // std::cout<<"polmag =  "<<polmag<<std::endl;

    //std::cout<<"LambdaPolarization generated (X,Y,Z) "<<polarizationVec.X()<<" "<<polarizationVec.Y()<<" "<<polarizationVec.Z()<<std::endl;
    delete rand;
    return polarizationVec;
}


// =============================================================================
// RANDOM HELPERS
// =============================================================================

// --- FIX: no per-call TRandom3(0) allocations ---
// OLD was:
//   TRandom3 *rand = new TRandom3(0); // Seed with 0 for reproducibility
//   Double_t randomVal = rand->Gaus(fMean, fSigma);
//   delete rand;
// That repeats the same values and is expensive.
Double_t get_random_value_rand(TRandom3 &rand, Double_t fMean, Double_t fSigma)
{
  return rand.Gaus(fMean, fSigma);
}

// keep Δ in [0, 2π)
Double_t wrap_delta_0_2pi(Double_t d)
{
  while (d < 0) d += TMath::TwoPi();
  while (d >= TMath::TwoPi()) d -= TMath::TwoPi();
  return d;
}


// =============================================================================
// Build chain subset from wildcard pattern (SLURM chunking)
// =============================================================================
TChain* build_chain_subset(const TString& pattern, Long64_t fileStart, Long64_t fileCount)
{
  // Expand pattern into list of files using temporary chain
  TChain tmp("events");
  Long64_t matched = tmp.Add(pattern);
  std::cout << "Pattern matched " << matched << " files" << std::endl;

  TObjArray* fileList = tmp.GetListOfFiles();
  if (!fileList || fileList->GetEntries() <= 0)
  {
    std::cerr << "ERROR: pattern expansion returned 0 files." << std::endl;
    return nullptr;
  }

  Long64_t nAll = fileList->GetEntries();
  Long64_t start = (fileStart < 0) ? 0 : fileStart;
  Long64_t end   = (fileCount < 0) ? nAll : TMath::Min(start + fileCount, nAll);

  if (start >= end)
  {
    std::cerr << "ERROR: empty file range: [" << start << "," << end << ") out of nAll=" << nAll << std::endl;
    return nullptr;
  }

  TChain* ch = new TChain("events");
  for (Long64_t i = start; i < end; ++i)
  {
    TChainElement* el = (TChainElement*)fileList->At(i);
    if (!el) continue;
    ch->Add(el->GetTitle());
  }

  std::cout << "Using files [" << start << "," << end-1 << "] => " << (end-start) << " files in this job" << std::endl;
  return ch;
}

// =============================================================================
// Yield hist selection + phi kinematics generation
// =============================================================================

TH2* pick_yield_hist(TFile* Phi_yield, Double_t centrality)
{
  TH2 *h_y_pT = nullptr;

  if (!Phi_yield)
  {
    std::cout << "ERROR: No yield file found!" << std::endl;
    return nullptr;
  }

  if (centrality < 10.0)
    h_y_pT = (TH2*)Phi_yield->Get("h2PhiYpT010");
  else if (centrality >= 10.0 && centrality < 40.0)
    h_y_pT = (TH2*)Phi_yield->Get("h2PhiYpT1040");
  else
    h_y_pT = (TH2*)Phi_yield->Get("h2PhiYpT40100");

  return h_y_pT;
}

//
// Lambda directed flow parameterization:
//   v1^Lambda(sNN, centrality, lambda_pT, lambda_y)
// Based on your cent–pT–y function.
//
// Signature:
//   get_v1_lambda(sNN, centrality, lambda_pT, lambda_y)
//
// Returns:
//   v1^Lambda with |v1| limited to [-1, +1].

// getters.cpp — REPLACE your current get_v1_lambda with this version

Double_t get_v1_lambda(Double_t sNN,
  Double_t centrality,
  Double_t lambda_pT,
  Double_t lambda_y)
{

  // Precompute simple integer powers for efficiency
  Double_t lambda_pT2 = lambda_pT * lambda_pT;
  Double_t lambda_pT3 = lambda_pT2 * lambda_pT;
  Double_t lambda_y3  = lambda_y * lambda_y * lambda_y;

  Double_t v1 =
  (28.8635 / TMath::Power(sNN, 2.89092)) *
  (
  // First term: linear in lambda_y (rapidity-odd)
  (-0.0233 * centrality + 0.5413 * TMath::Power(centrality, 1.0/3.0)) *
  (163.536 / 18.0188 - 163.536 / (lambda_pT + 18.0188)) *
  lambda_y
  +
  // Second term: cubic in lambda_y (rapidity-odd)
  (-0.0056 * centrality + 0.377 * TMath::Power(centrality, 1.0/3.0)) *
  (0.6653 * lambda_pT
  - 0.6172 * lambda_pT2
  + 0.1154 * lambda_pT3) *
  lambda_y3
);

// Safety: enforce |v1| <= 1
if (v1 > 1.0)  v1 = 1.0;
if (v1 < -1.0) v1 = -1.0;

return v1;
}

// Helper: compute v1_phi from Lambda v1 (simple scaling, no pT rescaling)
//
// f_scale is a tunable factor (e.g. 0.7 or 1.0) controlling how large
// v1^phi is relative to v1^Lambda.

Double_t get_v1_phi_from_lambda(Double_t pT_phi,
  Double_t y_phi,
  Double_t centrality,
  Double_t sNN,
  Double_t f_scale)
{
  // Use the same pT for Lambda and phi (Option B)
  Double_t pT_Lambda_eff = pT_phi;

  // Lambda v1 at (pT_Lambda_eff, y_phi, centrality, sNN)
  Double_t v1_Lambda = get_v1_lambda(sNN, centrality, pT_Lambda_eff, y_phi);

  // Scale to get phi v1
  Double_t v1_phi = f_scale * v1_Lambda;

  // Safety: enforce |v1| <= 1
  if (v1_phi > 1.0)  v1_phi = 1.0;
  if (v1_phi < -1.0) v1_phi = -1.0;

  return v1_phi;
}

// Wrapper: return v2_Lambda(pT,y,cen,sNN)
Double_t get_v2_lambda(Double_t lambda_pT,
  Double_t lambda_y,
  Double_t centrality,
  Double_t sNN)
{
// Just call your existing function
Double_t v2 = get_V2(sNN, centrality, lambda_pT, lambda_y);

// Safety: limit |v2| <= 1 (optional but helps avoid extreme values)
if (v2 > 1.0)  v2 = 1.0;
if (v2 < -1.0) v2 = -1.0;

return v2;
}

// Helper: NCQ-scaling to get v2_phi from Lambda v2
Double_t get_v2_phi_from_lambda(Double_t pT_phi,
           Double_t y_phi,
           Double_t centrality,
           Double_t sNN)
{
// Map phi pT to "effective" Lambda pT using NCQ scaling
Double_t pT_Lambda_eff = pT_phi * 3.0 / 2.0;

// Lambda v2 at (pT_Lambda_eff, y_phi, centrality, sNN)
Double_t v2_Lambda = get_v2_lambda(pT_Lambda_eff, y_phi, centrality, sNN);

// phi v2 via NCQ scaling:
Double_t v2_phi = (2.0 / 3.0) * v2_Lambda;

// Safety: limit |v2| <= 1
if (v2_phi > 1.0)  v2_phi = 1.0;
if (v2_phi < -1.0) v2_phi = -1.0;

return v2_phi;
}