#include "set_lambda_parameterization.hpp"

void set_lambda_parameterization(TFile *Lambda_yield,
                                 Double_t fBVal,
                                 UParticle &ULambda,
                                 TRandom3 &rand)
{ // Valeriy's function for properly lambda generation

    Double_t centrality = get_centrality(fBVal); // centrality for parameterization (b->centrality for Xe+Xe below)
    Double_t sNN = 2.87;                         // Energy of the collision in center-of-mass system

    TH2F *h_pt_y;
    if (centrality < 10)
        h_pt_y = (TH2F *)Lambda_yield->Get("h2PartYpT010");
    else if (centrality > 10 && centrality < 40)
        h_pt_y = (TH2F *)Lambda_yield->Get("h2PartYpT");
    else
        h_pt_y = (TH2F *)Lambda_yield->Get("h2PartYpT40100");

    Double_t lambda_pT; // pT of Lambda from pT-y TH2F
    Double_t lambda_y;  // rapidity of Lambda from pT-y TH2F

    h_pt_y->GetRandom2(lambda_y, lambda_pT, &rand);

    Double_t v1 = (28.8635 / (TMath::Power(sNN, 2.89092))) * ((-0.0233 * centrality + 0.5413 * TMath::Power(centrality, 1. / 3)) * (163.536 / 18.0188 - 163.536 / (lambda_pT + 18.0188)) * lambda_y + (-0.0056 * centrality + 0.377 * TMath::Power(centrality, 1. / 3)) * (0.6653 * lambda_pT - 0.6172 * TMath::Power(lambda_pT, 2) + 0.1154 * TMath::Power(lambda_pT, 3)) * TMath::Power(lambda_y, 3));                                                                                                                                       // v1 cent-pT-y func
    Double_t v2 = 1.05 * (0.8132 - 11.4 / TMath::Power(sNN, 2.1)) * ((-0.01105 * centrality + 0.000162 * TMath::Power(centrality, 2)) / (-0.01105 * 25 + 0.000162 * TMath::Power(25, 2))) * (-1) * ((0.32172 * TMath::Power(lambda_y, 1.805) - 0.1578) * lambda_pT + (0.05838 * TMath::Power(lambda_y, 3.2172) - 0.0179) * TMath::Exp(lambda_pT * (1.4562 * TMath::Power(lambda_y, 3.28814) + 0.22912)) * TMath::Sin(lambda_pT * (-1.004 * TMath::Power(lambda_y, 4.6398) + 3.88097) + (-3.6075 * TMath::Power(lambda_y, 4.6582) + 3.35966))); // v2 snn-cent-pT-y func

    if (v1 > 1)
        v1 = 1;
    else if (v1 < -1)
        v1 = -1;

    // v1 = 0.0;   //generation WITHOUT flows
    // v2 = 0.0;

    static TF1 f("f", "[0]*(1+2*[1]*TMath::Cos(x)+2*[2]*TMath::Cos(2*x))+[3]", 0, 2 * TMath::Pi());
    Double_t a1 = 1.0 + 2.0 * v1 + 2.0 * v2;
    Double_t a2 = 1.0 - 2.0 * v1 + 2.0 * v2;
    // Double_t a3 = 1 - v1 * v1 / (4 * v2) - 2 * v2; // maybe the reason is why f.SetParameter() become NaN

    Double_t a = (a1 < a2) ? a1 : a2;
    if (v2 != 0.0)
    {
        Double_t cos_min = -v1 / (4.0 * v2);
        if (cos_min >= -1.0 && cos_min <= 1.0)
        {
            Double_t a3 = 1.0 - v1 * v1 / (4.0 * v2) - 2.0 * v2;
            if (a3 < a)
                a = a3;
        }
    }
    if (a > 0.0)
        a = 0.0;                                         // Only shift if the minimum is actually below 0

    f.SetParameter(0, 1 / (2 * TMath::Pi() * (1 - a)));  // norm
    f.SetParameter(1, v1);                               // v1
    f.SetParameter(2, v2);                               // v2
    f.SetParameter(3, -a / (2 * TMath::Pi() * (1 - a))); // shift to have probability
    // f.SetNpx(10000);                                     // to get a better result when using TF1::GetRandom
    // f.SetNpx(100); // fixed: to speed up the execution time by a factor of 100

    Double_t phi = f.GetRandom(&rand);

    // Double_t fEnergyLambda = ULambda.GetMomentum().E();
    // Calculate total energy (E) properly
    TLorentzVector vec;
    // Calculate transverse mass (m_T)
    Double_t lambda_mass = 1.115683; // GeV/c² (Lambda mass)
    Double_t mT = sqrt(lambda_pT * lambda_pT + lambda_mass * lambda_mass);
    // Convert rapidity (y) to pseudorapidity (η)
    Double_t pz = mT * sinh(lambda_y); // longitudinal momentum
    Double_t lambda_eta = (pz != 0.0) ? TMath::ATanH(pz / sqrt(pz * pz + lambda_pT * lambda_pT)) : 0.0;
    Double_t fEnergyLambda = sqrt(lambda_pT * lambda_pT * cosh(lambda_eta) * cosh(lambda_eta) + lambda_mass * lambda_mass);
    vec.SetPtEtaPhiE(lambda_pT, lambda_eta, phi, fEnergyLambda);
    ULambda.SetMomentum(vec); // ULambda
}
