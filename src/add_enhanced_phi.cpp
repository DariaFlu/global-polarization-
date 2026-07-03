// add_enhanced_phi.C
//
// Goal:
//    selected from wildcard pattern (fileStart, fileCount)
//
// IMPORTANT FIX (serious remark):
//  - Now RNG is created once in main and passed by reference everywhere.
//
// Run example (local):
//   root -l -b -q 'add_enhanced_phi.C+("out_000000-000099.root", 0, 100)'
//
// Run example (SLURM array):
//   FILES_PER_TASK=100
//   START=$((SLURM_ARRAY_TASK_ID * FILES_PER_TASK))
//   root -l -b -q "add_enhanced_phi.C+(\"out_${START}.root\", ${START}, ${FILES_PER_TASK})"
//
// Notes:
//  - This macro still uses yield QA file (qa_xexe_phi.root) for (y,pT) sampling.
//  - Here just splited the input files and write independent output ROOTs.

#include "add_enhanced_phi.hpp"

// =============================================================================
// MAIN EXECUTION FUNCTION (SLURM-ready)
// =============================================================================
int add_enhanced_phi(const TString &outputFile,
                     Long64_t fileStart,
                     Long64_t fileCount,
                     const TString &pattern,
                     const TString &qaFile,
                     Int_t enhanceStat,
                     UInt_t seed)
{

  // --- Seed logic ---
  // If seed==0, derive from SLURM_ARRAY_TASK_ID (so each array task differs)
  if (seed == 0)
  {
    const char *a = gSystem->Getenv("SLURM_ARRAY_TASK_ID");
    if (a)
      seed = 12345u + (UInt_t)atoi(a);
    else
      seed = 12345u;
  }

  TRandom3 rand(seed); // Random number generator with fixed seed
  std::cout << "RNG seed = " << seed << std::endl;

  // --- Yield file ---
  TFile *Phi_yield = TFile::Open(qaFile, "READ");
  if (!Phi_yield || Phi_yield->IsZombie())
  {
    std::cerr << "ERROR: cannot open " << qaFile << std::endl;
    return 1;
  }

  // Input chain: subset of files
  TChain *inChain = build_chain_subset(pattern, fileStart, fileCount);
  if (!inChain)
  {
    Phi_yield->Close();
    return 2;
  }

  // --- Create or open output ROOT file ---
  // Try to open in UPDATE mode first (append to existing file)
  TFile *outFile = TFile::Open(outputFile, "RECREATE");
  if (!outFile || outFile->IsZombie())
  {
    // If file cannot be opened or is corrupted, recreate it
    delete outFile;
    outFile = TFile::Open(outputFile, "RECREATE");
    std::cout << "Creating new output file: " << outputFile << std::endl;
  }

  // --- Copy URun object from first input file (run-level information) ---
  // Get the name of the first file from TChain
  TObjArray *fileElements = inChain->GetListOfFiles();
  if (fileElements && fileElements->GetEntries() > 0)
  {
    TChainElement *chEl = (TChainElement *)fileElements->At(0);
    TString firstFileName = chEl->GetTitle();

    TFile *fFirst = TFile::Open(firstFileName, "READ");
    if (fFirst && !fFirst->IsZombie())
    {
      URun *inRun = nullptr;
      fFirst->GetObject("run", inRun);
      if (inRun)
      {
        // Write run object into the output file (overwrite if exists)
        outFile->cd();
        inRun->Write("run", TObject::kOverwrite);
        std::cout << "Copied URun from " << firstFileName << std::endl;
      }
      else
      {
        std::cerr << "WARNING: URun object not found in "
                  << firstFileName << std::endl;
      }
      fFirst->Close();
    }
  }

  // --- Event-level objects ---
  UEvent *inEvent = nullptr;       // Input event (from UniGen format)
  UEvent *outEvent = new UEvent(); // Output event that will contain original + enhanced particles
  inChain->SetBranchAddress("event", &inEvent);

  outFile->cd();

  // Get or create TTree
  // NOTE: Tree name cannot contain spaces if you want easy TChain later
  TTree *outTree = new TTree("events", "UniGen events (enhanced phi)");
  outTree->Branch("event", "UEvent", &outEvent, 32000, 99);

  // Dummy array for UParticle constructor children field
  Int_t child_null[2] = {0, 0}; // Children indices placeholder for UParticle

  // Temporary particle object used as container for Phi
  // (keep as in your original style)
  UParticle phi(1, 1, 1, 1, 1, 1, 1, child_null,
                1., 1., 1., 1., 1., 1., 1., 1., 1.);

  // --- List of variables ---
  // Int_t enhanceStat = 2;  // now passed as argument

  // --- List of histos ---
  TH1D *hPhiVarphi = new TH1D("hPhiVarphi", "#phi #varphi;#varphi [rad];Counts", 50, 0, 2. * TMath::Pi());
  TProfile *hPhiCos1 = new TProfile("hPhiCos1", "<cos(#Delta#phi)>", 50, 0, 2. * TMath::Pi(), -1., 1.);
  TProfile *hPhiCos2 = new TProfile("hPhiCos2", "<cos2(#Delta#phi)>", 50, 0, 2. * TMath::Pi(), -1., 1.);

  // --- Main event loop over input chain ---
  Long64_t nEvents = inChain->GetEntries();
  std::cout << "Events: " << nEvents << std::endl;

  for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++)
  {
    if (iEvent % 200 == 0)
    {
      std::cout << "Event No " << iEvent << " / " << nEvents << std::endl;
    }

    // Load event from input chain
    inChain->GetEntry(iEvent);

    outEvent->Clear(); // Remove particles from previous iteration

    // Copy event header information (global properties)
    outEvent->SetEventNr(inEvent->GetEventNr());
    outEvent->SetB(inEvent->GetB());
    outEvent->SetPhi(inEvent->GetPhi());
    outEvent->SetNes(inEvent->GetNes());
    outEvent->SetStepNr(inEvent->GetStepNr());
    outEvent->SetStepT(inEvent->GetStepT());

    // --- First pass: count phi in the event ---
    Int_t phiCounter = 0;
    for (Int_t i = 0; i < inEvent->GetNpa(); i++)
    {
      UParticle *part = inEvent->GetParticle(i);
      if (part && part->GetPdg() == 333)
        phiCounter++;
    }

    // --- Second pass: loop over all particles and process phi ---
    for (Int_t i = 0; i < inEvent->GetNpa(); i++)
    {
      UParticle *part = inEvent->GetParticle(i);
      if (!part)
        continue;

      // If this is NOT a phi, simply copy it to output
      if (part->GetPdg() != 333)
      {
        outEvent->AddParticle(*part);
        continue;
      }

      // If this is a phi (PDG 333) - reparameterize it
      phi = *part;
      set_phi_parameterization(Phi_yield, inEvent->GetB(), phi, rand);
      outEvent->AddParticle(phi);

      // ---- histos (kept) ----
      const Double_t phiAbs = get_positive_phi(phi.GetMomentum().Phi());
      const Double_t dphi = wrap_delta_0_2pi(phiAbs - inEvent->GetPhi());
      hPhiVarphi->Fill(phiAbs);
      hPhiCos1->Fill(dphi, std::cos(dphi));
      hPhiCos2->Fill(dphi, std::cos(2.0 * dphi));
    }

    // If there were no phi's at all, we artificially create ONE with random position/momentum.
    // --- FIX: do it once outside the loop and do NOT leak memory ---
    if (phiCounter == 0)
    {
      // Random position (x, y, z, t) for the new phi
      TLorentzVector newPhiPos(
          get_random_value_rand(rand, 0, 12),
          get_random_value_rand(rand, 0, 12),
          get_random_value_rand(rand, 0, 20),
          0);

      // Dummy momentum (will be overwritten in set_phi_parameterization)
      TLorentzVector newPhiMom(1., 1., 1., 1.);

      // Index for newly added particle in event
      Int_t idx = outEvent->GetNpa();

      // Create artificial phi with PDG=333
      UParticle newPhi(idx, 333, 1, 1, 1, -15, -1, child_null, newPhiMom, newPhiPos, 1.0);

      set_phi_parameterization(Phi_yield, inEvent->GetB(), newPhi, rand);

      outEvent->AddParticle(newPhi);

      // ---- histos (kept) ----
      const Double_t phiAbs = get_positive_phi(newPhi.GetMomentum().Phi());
      const Double_t dphi = wrap_delta_0_2pi(phiAbs - inEvent->GetPhi());
      hPhiVarphi->Fill(phiAbs);
      hPhiCos1->Fill(dphi, std::cos(dphi));
      hPhiCos2->Fill(dphi, std::cos(2.0 * dphi));
    }

    // Total number of particles in output event (after copy + possible artificial phi)
    Int_t nParticles = outEvent->GetNpa();

    // --- Artificial enhancement of phi statistics ---
    Double_t fEnhanceStat = enhanceStat;
    while (fEnhanceStat > 0)
    {                          // Repeat until enhancement counter is exhausted
      Int_t enhancedFlag = -9; // Mark these as "enhanced" phis

      // Dummy random momentum (not really used; overwritten in parameterization)
      TLorentzVector mom_rand(1., 1., 1., 1.);

      // Randomized position for enhanced phi inside some small volume
      TLorentzVector pos_rand(
          get_random_value_rand(rand, 0, 0.03),
          get_random_value_rand(rand, 0, 0.03),
          get_random_value_rand(rand, 0, 0.03),
          get_random_value_rand(rand, 0, 0.03));

      // Index for newly added particles in event
      Int_t indexOfEnhancedParticle = nParticles++;

      // Enhanced phi
      // UniGen -> https://git.jinr.ru/nica/mpdroot/-/blob/dev/simulation/generators/unigenFormat/UParticle.h
      // Enhanced phi with fMate == -9
      UParticle enhancedPhi(indexOfEnhancedParticle, 333, 1, 1, 1, enhancedFlag, // Use Mate to tag enhanced particles
                            -1, child_null,
                            1., 1., 1., 1., 1., 1., 1., 1., 1.); // Dummy parameters

      enhancedPhi.SetPosition(pos_rand);

      // Sample Phi kinematics from parameterization/yield
      set_phi_parameterization(Phi_yield, inEvent->GetB(), enhancedPhi, rand);

      // Fill QA histograms for enhanced sample
      const Double_t phiAbs = get_positive_phi(enhancedPhi.GetMomentum().Phi());
      const Double_t dphi = wrap_delta_0_2pi(phiAbs - inEvent->GetPhi());
      hPhiVarphi->Fill(phiAbs);
      hPhiCos1->Fill(dphi, std::cos(dphi));
      hPhiCos2->Fill(dphi, std::cos(2.0 * dphi));

      // Add enhanced Phi and its products to output event and vectors
      outEvent->AddParticle(enhancedPhi);

      // Decrease enhancement counter
      fEnhanceStat--;
    }

    // Store current event into the output tree
    outTree->Fill();
  }

  outFile->cd();
  outTree->Write("", TObject::kOverwrite);

  delete outEvent;
  delete inChain;

  hPhiVarphi->Write();
  hPhiCos1->Write();
  hPhiCos2->Write();

  outFile->Close();
  Phi_yield->Close();

  std::cout << "Wrote output: " << outputFile << std::endl;
  return 0;
}
