#include "post_phi_flow.h"

Double_t get_centrality(Double_t fBVal)
{
  if (fBVal < 3.44)
    return 5;
  else if (fBVal < 4.88)
    return 15;
  else if (fBVal < 5.84)
    return 25;
  else if (fBVal < 6.64)
    return 35;
  else if (fBVal < 7.44)
    return 45;
  else if (fBVal < 8.08)
    return 55;
  else if (fBVal < 8.72)
    return 65;
  else if (fBVal < 9.36)
    return 75;
  else if (fBVal < 9.84)
    return 85;
  else
    return 95;
}

static int centClass(Double_t cent)
{
  if (cent < 10.0)
    return 0; // 0-10
  if (cent < 40.0)
    return 1; // 10-40
  return 2;   // 40-100
}
static const char *centName(int ic)
{
  if (ic == 0)
    return "c010";
  if (ic == 1)
    return "c1040";
  return "c40100";
}

static inline double wrapToPi(double x)
{
  // wrap to (-pi, pi]
  x = std::fmod(x, 2.0 * TMath::Pi());
  if (x <= -TMath::Pi())
    x += 2.0 * TMath::Pi();
  if (x > TMath::Pi())
    x -= 2.0 * TMath::Pi();
  return x;
}

int post_phi_flow(const char *inPattern,
                  const char *outFile,
                  int enhancedFactor)
{
  // -----------------------------
  // Input
  // -----------------------------
  TChain ch("events");
  Long64_t nFiles = ch.Add(inPattern);
  std::cout << "[post] Pattern matched " << nFiles << " files\n";
  if (ch.GetEntries() <= 0)
  {
    std::cerr << "[post][ERROR] chain has 0 entries. Check pattern.\n";
    return 1;
  }

  UEvent *ev = nullptr;
  ch.SetBranchAddress("event", &ev);

  // -----------------------------
  // Output + QA histos
  // -----------------------------
  TFile *fout = TFile::Open(outFile, "RECREATE");
  if (!fout || fout->IsZombie())
  {
    std::cerr << "[post][ERROR] cannot create output: " << outFile << "\n";
    return 2;
  }

  TH1D *hB = new TH1D("hB", "impact parameter; b (fm); events", 200, 0, 12);
  TH1D *hCent = new TH1D("hCent", "centrality; centrality (%); events", 20, 0, 100);

  // Kinematics QA for each species
  TH2D *hYpT_phi = new TH2D("hYpT_phi", "phi kinematics; y; p_{T} (GeV/c)", 120, -2, 2, 120, 0, 3);
  TH2D *hYpT_proton = new TH2D("hYpT_proton", "proton kinematics; y; p_{T} (GeV/c)", 120, -2, 2, 120, 0, 3);
  TH2D *hYpT_pion = new TH2D("hYpT_pion", "pion kinematics; y; p_{T} (GeV/c)", 120, -2, 2, 120, 0, 3);

  // QA for multiplicity(impact parameter) based on protons, pions and phi
  // NOTE: we store <N_species> vs b (event-averaged multiplicity in each b-bin)
  TProfile *hMult_proton = new TProfile("hMult_proton", "proton multiplicity vs b; b (fm); N_{p}", 200, 0, 12);
  TProfile *hMult_pion = new TProfile("hMult_pion", "pion multiplicity vs b; b (fm); N_{#pi^{+}}+N_{#pi^{-}}", 200, 0, 12);
  TProfile *hMult_phi = new TProfile("hMult_phi", "phi multiplicity vs b; b (fm); N_{#phi}", 200, 0, 12);

  // Per centrality class:
  //  - 1D profiles vs y and pT
  // (same logic as for phi, but added for proton/pion, and phi renamed to *_phi_*)
  TProfile *pV1_y_phi[3];
  TProfile *pV2_y_phi[3];
  TProfile *pV1_pT_phi[3];
  TProfile *pV2_pT_phi[3];

  TProfile *pV1_y_proton[3];
  TProfile *pV2_y_proton[3];
  TProfile *pV1_pT_proton[3];
  TProfile *pV2_pT_proton[3];

  TProfile *pV1_y_pion[3];
  TProfile *pV2_y_pion[3];
  TProfile *pV1_pT_pion[3];
  TProfile *pV2_pT_pion[3];

  for (int ic = 0; ic < 3; ++ic)
  {
    TString tag = centName(ic);

    // phi (renamed)
    // pV1_y_phi[ic]   = new TProfile("pV1_y_phi_"+tag,   "v_{1}^{#phi}(y); y; <cos(#Delta)>",   50, -2.0, 2.0);
    // pV2_y_phi[ic]   = new TProfile("pV2_y_phi_"+tag,   "v_{2}^{#phi}(y); y; <cos2#Delta>",    50, -2.0, 2.0);
    pV1_y_phi[ic] = new TProfile("pV1_y_phi_" + tag, "v_{1}^{#phi}(y); y; <cos(#varphi-#Psi_{RP})>", 50, -1.0, 1.0);
    pV2_y_phi[ic] = new TProfile("pV2_y_phi_" + tag, "v_{2}^{#phi}(y); y; <cos2(#varphi-#Psi_{RP})>", 50, -1.0, 1.0);

    pV1_pT_phi[ic] = new TProfile("pV1_pT_phi_" + tag, "v_{1}^{#phi}(p_{T}); p_{T}; <cos(#varphi-#Psi_{RP})>", 50, 0.4, 2.);
    pV2_pT_phi[ic] = new TProfile("pV2_pT_phi_" + tag, "v_{2}^{#phi}(p_{T}); p_{T}; <cos2(#varphi-#Psi_{RP})>", 50, 0.4, 2.);

    // proton
    // pV1_y_proton[ic]   = new TProfile("pV1_y_proton_"+tag,   "v_{1}^{p}(y); y; <cos(#Delta)>",   50, -2.0, 2.0);
    // pV2_y_proton[ic]   = new TProfile("pV2_y_proton_"+tag,   "v_{2}^{p}(y); y; <cos2#Delta>",    50, -2.0, 2.0);
    pV1_y_proton[ic] = new TProfile("pV1_y_proton_" + tag, "v_{1}^{p}(y); y; <cos(#varphi-#Psi_{RP})>", 50, -1.0, 1.0);
    pV2_y_proton[ic] = new TProfile("pV2_y_proton_" + tag, "v_{2}^{p}(y); y; <cos2(#varphi-#Psi_{RP})>", 50, -1.0, 1.0);
    pV1_pT_proton[ic] = new TProfile("pV1_pT_proton_" + tag, "v_{1}^{p}(p_{T}); p_{T}; <cos(#varphi-#Psi_{RP})>", 50, 0.4, 2.);
    pV2_pT_proton[ic] = new TProfile("pV2_pT_proton_" + tag, "v_{2}^{p}(p_{T}); p_{T}; <cos2(#varphi-#Psi_{RP})>", 50, 0.4, 2.);

    // pion (pi+ + pi- together)
    pV1_y_pion[ic] = new TProfile("pV1_y_pion_" + tag, "v_{1}^{#pi}(y); y; <cos(#varphi-#Psi_{RP})>", 50, -1.0, 1.0);
    pV2_y_pion[ic] = new TProfile("pV2_y_pion_" + tag, "v_{2}^{#pi}(y); y; <cos2(#varphi-#Psi_{RP})>", 50, -1.0, 1.0);

    pV1_pT_pion[ic] = new TProfile("pV1_pT_pion_" + tag, "v_{1}^{#pi}(p_{T}); p_{T}; <cos(#varphi-#Psi_{RP})>", 50, 0.4, 2.);
    pV2_pT_pion[ic] = new TProfile("pV2_pT_pion_" + tag, "v_{2}^{#pi}(p_{T}); p_{T}; <cos2(#varphi-#Psi_{RP})>", 50, 0.4, 2.);
  }

  // -----------------------------
  // Event loop
  // -----------------------------
  const Long64_t nEv = ch.GetEntries();
  std::cout << "[post] Events = " << nEv << "\n";
  std::cout << "[post] enhancedFactor = " << enhancedFactor << "\n";

  Long64_t nPhiUsed = 0;
  Long64_t nProtUsed = 0;
  Long64_t nPionUsed = 0;

  // counts of ACCEPTED phi types (global)
  Long64_t nPhiOrigUsed = 0;
  Long64_t nPhiEnhUsed = 0;

  for (Long64_t iEv = 0; iEv < nEv; ++iEv)
  {
    if (iEv % 200 == 0)
      std::cout << "[post] " << iEv << " / " << nEv << "\n";

    ch.GetEntry(iEv);
    if (!ev)
      continue;

    const double b = ev->GetB();
    const double cent = get_centrality(b);
    const int ic = centClass(cent);

    // TRUE reaction plane (UniGen stores it as event phi)
    const double psiRP = ev->GetPhi(); // radians

    hB->Fill(b);
    hCent->Fill(cent);

    // event multiplicities (filled into <N>(b) profiles)
    int nProtonEv = 0;
    int nPionEv = 0;
    int nPhiEv = 0;

    const int nPa = ev->GetNpa();
    int enhLeft = enhancedFactor;
    for (int i = 0; i < nPa; ++i)
    {
      UParticle *p = ev->GetParticle(i);
      if (!p)
        continue;

      const int pdg = p->GetPdg();

      const TLorentzVector mom = p->GetMomentum();
      const double pT = mom.Pt();
      const double y = mom.Rapidity();
      const double phi = mom.Phi();

      // Δ = φ - ΨRP
      const double d = wrapToPi(phi - psiRP);
      const double c1 = std::cos(d);
      const double c2 = std::cos(2.0 * d);

      // -----------------------------
      // phi meson
      // -----------------------------
      if (pdg == 333)
      {
        // enhanced selection (tag: Mate == -9)
        const bool isEnhanced = (p->GetMate() == -9);

        // Keep ALL originals
        if (isEnhanced)
        {
          // Keep only a controlled number of enhanced per event
          if (enhLeft <= 0) continue;
          nPhiEnhUsed++;
          enhLeft--;
          hYpT_phi->Fill(y, pT);

          pV1_y_phi[ic]->Fill(y, c1);
          pV2_y_phi[ic]->Fill(y, c2);
  
          pV1_pT_phi[ic]->Fill(pT, c1);
          pV2_pT_phi[ic]->Fill(pT, c2);
        }
        else
        {
          nPhiOrigUsed++;
          hYpT_phi->Fill(y, pT);

          pV1_y_phi[ic]->Fill(y, c1);
          pV2_y_phi[ic]->Fill(y, c2);
  
          pV1_pT_phi[ic]->Fill(pT, c1);
          pV2_pT_phi[ic]->Fill(pT, c2);
        }

        // hYpT_phi->Fill(y, pT);

        // pV1_y_phi[ic]->Fill(y, c1);
        // pV2_y_phi[ic]->Fill(y, c2);

        // pV1_pT_phi[ic]->Fill(pT, c1);
        // pV2_pT_phi[ic]->Fill(pT, c2);

        nPhiUsed++;
        nPhiEv++;
        continue;
      }

      // -----------------------------
      // proton / antiproton
      // -----------------------------
      if (std::abs(pdg) == 2212)
      {
        hYpT_proton->Fill(y, pT);

        pV1_y_proton[ic]->Fill(y, c1);
        pV2_y_proton[ic]->Fill(y, c2);

        pV1_pT_proton[ic]->Fill(pT, c1);
        pV2_pT_proton[ic]->Fill(pT, c2);

        nProtUsed++;
        nProtonEv++;
        continue;
      }

      // -----------------------------
      // pion (pi+ + pi-)
      // -----------------------------
      if (std::abs(pdg) == 211)
      {
        hYpT_pion->Fill(y, pT);

        pV1_y_pion[ic]->Fill(y, c1);
        pV2_y_pion[ic]->Fill(y, c2);

        pV1_pT_pion[ic]->Fill(pT, c1);
        pV2_pT_pion[ic]->Fill(pT, c2);

        nPionUsed++;
        nPionEv++;
        continue;
      }
    }

    // fill multiplicity QA profiles per event
    hMult_proton->Fill(b, nProtonEv);
    hMult_pion->Fill(b, nPionEv);
    hMult_phi->Fill(b, nPhiEv);
  }

  // std::cout << "[post] phi used    = " << nPhiUsed << "\n";
  // std::cout << "[post] proton used = " << nProtUsed << "\n";
  // std::cout << "[post] pion used   = " << nPionUsed << "\n";

  std::cout << "[post] phi used total      = " << nPhiUsed << "\n";
  std::cout << "[post] phi orig used     = " << nPhiOrigUsed << "\n";
  std::cout << "[post] phi enhanced used = " << nPhiEnhUsed << "\n";
  std::cout << "[post] proton used         = " << nProtUsed << "\n";
  std::cout << "[post] pion used           = " << nPionUsed << "\n";

  // -----------------------------
  // Write
  // -----------------------------
  fout->cd();
  hB->Write();
  hCent->Write();

  hYpT_phi->Write();
  hYpT_proton->Write();
  hYpT_pion->Write();

  hMult_proton->Write();
  hMult_pion->Write();
  hMult_phi->Write();

  for (int ic = 0; ic < 3; ++ic)
  {
    pV1_y_phi[ic]->Write();
    pV2_y_phi[ic]->Write();
    pV1_pT_phi[ic]->Write();
    pV2_pT_phi[ic]->Write();

    pV1_y_proton[ic]->Write();
    pV2_y_proton[ic]->Write();
    pV1_pT_proton[ic]->Write();
    pV2_pT_proton[ic]->Write();

    pV1_y_pion[ic]->Write();
    pV2_y_pion[ic]->Write();
    pV1_pT_pion[ic]->Write();
    pV2_pT_pion[ic]->Write();
  }

  fout->Close();
  std::cout << "[post] Done. Wrote " << outFile << "\n";
  return 0;
}
