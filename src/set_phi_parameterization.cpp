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

    // Beam energy (center-of-mass)
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
    else if (centrality >= 10.0 && centrality < 40.0)
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
    // 4. Convert rapidity y to pseudorapidity eta and compute energy
    // ----------------------
    const Double_t phi_mass = 1.019461; // GeV/c^2
    Double_t mT   = std::sqrt(phi_pT * phi_pT + phi_mass * phi_mass); // Transverse mass
    Double_t pz   = mT * std::sinh(phi_y);                            // Longitudinal momentum
    Double_t p_abs = std::sqrt(pz * pz + phi_pT * phi_pT);            // Total momentum

    // Pseudorapidity eta (guard against division by zero)
    Double_t phi_eta = 0.0;
    if (p_abs > 0.0 && (p_abs - pz) > 0.0)
        phi_eta = 0.5 * std::log((p_abs + pz) / (p_abs - pz));

    // Total energy E_phi
    Double_t E_phi =
        std::sqrt(phi_pT * phi_pT * std::cosh(phi_eta) * std::cosh(phi_eta) +
                  phi_mass * phi_mass);

    // ----------------------
    // 5. Compute v1 and v2 at sampled (pT, y, centrality)
    // ----------------------
    // v1(phi): from Lambda via simple scaling
    const Double_t f_scale_v1_phi = 1.0; // can later tune, e.g. 0.7
    Double_t v1 = get_v1_phi_from_lambda(phi_pT, phi_y, centrality, sNN, f_scale_v1_phi);

    // v2(phi): from Lambda via NCQ scaling
    Double_t v2 = get_v2_phi_from_lambda(phi_pT, phi_y, centrality, sNN);

    // Note: v1 and v2 are clamped inside the getter functions to |v_n| <= 1

    // ----------------------
    // 6. Build azimuthal distribution dN/dphi and sample phi_az
    //    dN/dphi ∝ 1 + 2*v1*cos(phi) + 2*v2*cos(2*phi), shifted positive and normalized
    // ----------------------

    // Shape function: g(x) = 1 + 2*v1*cos(x) + 2*v2*cos(2*x)
    static TF1 f_phi_shape("f_phi_shape",
                           "1+2*[0]*TMath::Cos(x)+2*[1]*TMath::Cos(2*x)",
                           0.0, 2.0 * TMath::Pi());
    f_phi_shape.SetParameter(0, v1);
    f_phi_shape.SetParameter(1, v2);
    f_phi_shape.SetNpx(1000);

    // Numerical minimum of g(x) on [0, 2π]
    Double_t x_min = f_phi_shape.GetXmin();
    Double_t x_max = f_phi_shape.GetXmax();
    Double_t f_min = f_phi_shape.GetMinimum(x_min, x_max);

    // Positive function: g_pos(x) = g(x) - f_min
    static TF1 f_phi_pos("f_phi_pos",
                         "1+2*[0]*TMath::Cos(x)+2*[1]*TMath::Cos(2*x)-[2]",
                         0.0, 2.0 * TMath::Pi());
    f_phi_pos.SetParameter(0, v1);
    f_phi_pos.SetParameter(1, v2);
    f_phi_pos.SetParameter(2, f_min);
    f_phi_pos.SetNpx(1000);

    // Normalization: N = ∫ g_pos(x) dx over [0, 2π]
    Double_t norm = f_phi_pos.Integral(0.0, 2.0 * TMath::Pi());
    if (norm <= 0.0)
        norm = 1.0; // safety fallback

    // Final normalized PDF: f_norm(x) = g_pos(x) / norm
    static TF1 f_phi_norm("f_phi_norm",
                          "(1+2*[0]*TMath::Cos(x)+2*[1]*TMath::Cos(2*x)-[2])/[3]",
                          0.0, 2.0 * TMath::Pi());
    f_phi_norm.SetParameter(0, v1);
    f_phi_norm.SetParameter(1, v2);
    f_phi_norm.SetParameter(2, f_min);
    f_phi_norm.SetParameter(3, norm);
    f_phi_norm.SetNpx(1000);

    // Sample azimuthal angle
    Double_t phi_az = f_phi_norm.GetRandom(&rand);

    // ----------------------
    // 7. Set final four-vector for UPhi
    // ----------------------
    TLorentzVector vec;
    vec.SetPtEtaPhiE(phi_pT, phi_eta, phi_az, E_phi);

    UPhi.SetMomentum(vec);
}