#include "set_phi_parameterization.hpp"

//======================================================================
// Main entry: set_phi_parameterization
//======================================================================

void set_phi_parameterization(TFile *Phi_yield,
                              Double_t fBVal,
                              UParticle &UPhi,
                              TRandom3 &rand)
{
    // ----------------------
    // 1. Centrality from impact parameter
    // ----------------------
    const Double_t centrality = get_centrality(fBVal);

    // Beam energy
    const Double_t sNN = 9.2; // GeV, NICA/MPD Bi+Bi setting

    // ----------------------
    // 2. Select yield histogram for (y, pT)
    // ----------------------
    if (!Phi_yield)
    {
        std::cerr << "ERROR in set_phi_parameterization: Phi_yield file is null" << std::endl;
        return;
    }

    TH2 *h_y_pT = nullptr;
    if (centrality < 10.0)
        h_y_pT = (TH2 *)Phi_yield->Get("h2PhiYpT010");
    else if (centrality < 40.0)
        h_y_pT = (TH2 *)Phi_yield->Get("h2PhiYpT1040");
    else
        h_y_pT = (TH2 *)Phi_yield->Get("h2PhiYpT40100");

    if (!h_y_pT)
    {
        std::cerr << "ERROR in set_phi_parameterization: cannot find h2PhiYpT*** in Phi_yield"
                  << " (centrality = " << centrality << ")" << std::endl;
        return;
    }

    // ----------------------
    // 3. Sample (y, pT) from yield histogram
    // ----------------------
    Double_t phi_y  = 0.0;
    Double_t phi_pT = 0.0;
    h_y_pT->GetRandom2(phi_y, phi_pT, &rand);

    // ----------------------
    // 4. Kinematics from (y, pT)
    // ----------------------
    const Double_t phi_mass = 1.019461; // GeV/c^2

    const Double_t mT = std::sqrt(phi_pT * phi_pT + phi_mass * phi_mass);
    const Double_t pz = mT * std::sinh(phi_y);
    const Double_t E  = mT * std::cosh(phi_y);

    // ----------------------
    // 5. Compute v1 and v2 at sampled (pT, y, centrality)
    // ----------------------
    const Double_t f_scale_v1_phi = 1.0;
    Double_t v1 = get_v1_phi_from_lambda(phi_pT, phi_y, centrality, sNN, f_scale_v1_phi);
    Double_t v2 = get_v2_phi_from_lambda(phi_pT, phi_y, centrality, sNN);

    // Extra safety
    if (v1 > 1.0)  v1 = 1.0;
    if (v1 < -1.0) v1 = -1.0;

    if (v2 > 1.0)  v2 = 1.0;
    if (v2 < -1.0) v2 = -1.0;

    // ----------------------
    // 6. Linear azimuthal PDF:
    //    f(phi) = 1 + 2*v1*cos(phi) + 2*v2*cos(2phi)
    //
    // Keep the conventional Fourier form.
    // Use TF1 only to evaluate the function.
    // Sample phi by rejection method.
    // ----------------------
    static TF1 f_phi_shape("f_phi_shape",
                           "1.0 + 2.0*[0]*TMath::Cos(x) + 2.0*[1]*TMath::Cos(2.0*x)",
                           0.0, 2.0 * TMath::Pi());

    f_phi_shape.SetParameter(0, v1);
    f_phi_shape.SetParameter(1, v2);
    f_phi_shape.SetNpx(1000);

    // Check if PDF becomes negative
    const Double_t fmin = f_phi_shape.GetMinimum(0.0, 2.0 * TMath::Pi());

    if (fmin < 0.0)
    {
        std::cerr << "WARNING in set_phi_parameterization: azimuthal PDF is negative"
                  << " (v1 = " << v1 << ", v2 = " << v2
                  << ", fmin = " << fmin << "). "
                  << "Sampling will use f(phi)=0 in negative region." << std::endl;
    }

    // Safe upper bound for rejection sampling
    Double_t fmax = 1.0 + 2.0 * std::abs(v1) + 2.0 * std::abs(v2);
    if (fmax <= 0.0)
        fmax = 1.0;

    Double_t phi_az = 0.0;
    Bool_t accepted = kFALSE;

    const Int_t maxTrials = 100000;

    for (Int_t iTry = 0; iTry < maxTrials; ++iTry)
    {
        const Double_t phi_try = rand.Uniform(0.0, 2.0 * TMath::Pi());
        Double_t fval = f_phi_shape.Eval(phi_try);

        if (fval < 0.0)
            fval = 0.0;

        const Double_t y_try = rand.Uniform(0.0, fmax);

        if (y_try < fval)
        {
            phi_az = phi_try;
            accepted = kTRUE;
            break;
        }
    }

    if (!accepted)
    {
        std::cerr << "WARNING in set_phi_parameterization: rejection sampling failed,"
                  << " fallback to uniform phi" << std::endl;
        phi_az = rand.Uniform(0.0, 2.0 * TMath::Pi());
    }
    // static int nDebug = 0;
// if (nDebug < 20)
// {
//     std::cout << "DEBUG phi: pT=" << phi_pT
//               << " y=" << phi_y
//               << " v1=" << v1
//               << " v2=" << v2
//               << " phi_az=" << phi_az
//               << " cos2phi=" << std::cos(2.0 * phi_az)
//               << std::endl;
//     ++nDebug;
// }

    // ----------------------
    // 7. Build final four-vector
    // ----------------------
    const Double_t px = phi_pT * std::cos(phi_az);
    const Double_t py = phi_pT * std::sin(phi_az);

    TLorentzVector vec;
    vec.SetPxPyPzE(px, py, pz, E);

    UPhi.SetMomentum(vec);
}