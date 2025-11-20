
#include "simulate_lambda_decays.hpp"


// Main function to simulate Λ → p + π- decays and write them into an output tree.
// - inputFile  : UniGen input file(s) with "events" TTree
// - outputFile : ROOT file where "decays" tree is stored/updated
// - confInFile : ROOT file with parameterization / yields for Lambda
// - enhanceStat: factor to artificially enhance Lambda statistics
void simulate_lambda_decays(TString inputFile, TString outputFile, TString confInFile, Int_t enhanceStat) {

    // Input chain and output file/tree
    TChain *inChain = nullptr;
    TFile *outFile = nullptr;
    TTree *outTree = nullptr;

    // Variables stored in the output tree
    double cos_theta_p = 0;  // Proton angle relative to some axis (for anisotropy / polarization analysis)
    Int_t eventID;           // Event ID (to be filled if needed)

    // Vectors that will be written to the TTree as branches
    std::vector<UParticle>* ULambda = nullptr;  // Container for Lambda hyperons
    std::vector<UParticle>* UProton = nullptr;  // Container for decay protons
    std::vector<UParticle>* UPion = nullptr;    // Container for decay pions
    std::vector<ROOT::Math::XYZVector>* vecPolarization = nullptr; // Polarization vectors per particle/event
    std::vector<ROOT::Math::XYZVector>* vecPol = nullptr;          // Alternative polarization container (naming overlap)

    // --- Open input as TChain of "events" ---
    inChain = new TChain("events");
    inChain->Add(inputFile);

    // --- Create or open output ROOT file ---
    // Try to open in UPDATE mode first (append to existing file)
    outFile = TFile::Open(outputFile, "UPDATE");
    if (!outFile || outFile->IsZombie()) {
        // If file cannot be opened or is corrupted, recreate it
        delete outFile;
        outFile = TFile::Open(outputFile, "RECREATE");
        std::cout << "Creating new output file: " << outputFile << std::endl;
    }

    // --- Copy URun object from first input file (run-level information) ---
    // Get the name of the first file from TChain
    TObjArray* fileElements = inChain->GetListOfFiles();
    if (fileElements && fileElements->GetEntries() > 0) {
        TChainElement* chEl = (TChainElement*)fileElements->At(0);
        TString firstFileName = chEl->GetTitle();

        TFile* fFirst = TFile::Open(firstFileName, "READ");
        if (fFirst && !fFirst->IsZombie()) {
            URun* inRun = nullptr;
            fFirst->GetObject("run", inRun);
            if (inRun) {
                // Write run object into the output file (overwrite if exists)
                outFile->cd();
                inRun->Write("run", TObject::kOverwrite);
                std::cout << "Copied URun from " << firstFileName << std::endl;
            } else {
                std::cerr << "WARNING: URun object not found in " 
                          << firstFileName << std::endl;
            }
            fFirst->Close();
        }
    }

    // --- Event-level objects ---
    UEvent  *inEvent = nullptr;      // Input event (from UniGen format)
    UEvent  *outEvent = new UEvent(); // Output event that will contain original + generated decay particles
    inChain->SetBranchAddress("event", &inEvent);

    // Try to get existing "decays" tree (if appending)
    outTree = (TTree*)outFile->Get("decays");
    bool newTree = false;
    if (!outTree) {
        // No existing tree: create a new one
        outTree = new TTree("decays", "Lambda decay products");
        newTree = true;
        std::cout << "Creating new output tree" << std::endl;

        // --- Define branches for a newly created tree ---
        outTree->Branch("event", &outEvent);
        // outTree->Branch("run", &outRun); // If needed, run-level object could be added

        outTree->Branch("cos_theta_p", &cos_theta_p);
        outTree->Branch("eventID", &eventID, "eventID/I");
        outTree->Branch("ULambda", &ULambda);
        outTree->Branch("UProton", &UProton);
        outTree->Branch("UPion", &UPion);
        outTree->Branch("LambdaPolarization", &vecPol);
        outTree->Branch("Polarization", &vecPolarization); // TODO: check naming consistency with vecPol
    } else {
        // Existing tree found: append new entries
        std::cout << "Appending to existing tree with " << outTree->GetEntries() << " entries" << std::endl;

        // Connect branches to already existing tree structure
        outTree->SetBranchAddress("event", &outEvent);
        // outTree->SetBranchAddress("run", &outRun);

        outTree->SetBranchAddress("cos_theta_p", &cos_theta_p);
        outTree->SetBranchAddress("eventID", &eventID);
        outTree->SetBranchAddress("ULambda", &ULambda);
        outTree->SetBranchAddress("UProton", &UProton);
        outTree->SetBranchAddress("UPion", &UPion);
        outTree->SetBranchAddress("LambdaPolarization", &vecPol);
        outTree->SetBranchAddress("Polarization", &vecPolarization);
    }

    // --- Open configuration / yield file for Lambda parameterization ---
    TFile* Lambda_yield = TFile::Open(confInFile,"READ");

    // --- Histograms for QA / debugging ---
    TH1D *hLambdaPt  = new TH1D("hLambdaPt", "Lambda pT;pT [GeV/c];Counts", 100, 0, 2);
    TH1D *hCosTheta  = new TH1D("hCosTheta", "Proton cos(#theta);cos(#theta);Counts", 100, -1, 1);
    TH1D *hLambdaPhi = new TH1D("hLambdaPhi", "Lambda (#Delta#varphi);#Delta#varphi;Counts", 100, 0, 2.*TMath::Pi());
    TH2D *hAngVsPt   = new TH2D("hAngVsPt", "cos(#theta) vs Lambda pT;pT [GeV/c];cos(#theta)", 50, 0, 2, 50, -1, 1);

    // Dummy array for UParticle constructor children field
    Int_t dummy[2];

    // Temporary particle objects used as containers for Lambda and its decay products
    UParticle lambda(1, 1, 1, 1, 1, 1, 1, dummy,
                     1., 1., 1., 1., 1., 1., 1., 1., 1.);
    UParticle proton(1, 1, 1, 1, 1, 1, 1, dummy,
                     1., 1., 1., 1., 1., 1., 1., 1., 1.);
    UParticle pion(1, 1, 1, 1, 1, 1, 1, dummy,
                   1., 1., 1., 1., 1., 1., 1., 1., 1.);

    // --- Physics constants and parameters ---
    TRandom3* rand = new TRandom3(0);  // Random number generator with fixed seed
    const double mLambda = 1.115683;   // Λ mass [GeV/c^2]
    const double mProton = 0.938272;   // Proton mass [GeV/c^2]
    const double mPion   = 0.139570;   // Pion mass [GeV/c^2]
    Double_t fSigmaPolVal = 1.5;       // Width for polarization smearing (model parameter)
    const double anisotropy = 0.732;   // Decay anisotropy parameter (αΛ)
    Int_t child_null[2] = {0, 0};      // Children indices placeholder for UParticle

    // --- Main event loop over input chain ---
    Long64_t nEvents = inChain->GetEntries();
    std::cout << "Events: " << nEvents << std::endl;

    for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {

        // Load event from input chain
        inChain->GetEntry(iEvent);
        outEvent->Clear(); // Remove particles from previous iteration
        Int_t lambdaCounter = 0;

        // Copy event header information (global properties)
        outEvent->SetEventNr(inEvent->GetEventNr());
        outEvent->SetB(inEvent->GetB());
        outEvent->SetPhi(inEvent->GetPhi());
        outEvent->SetNes(inEvent->GetNes());
        outEvent->SetStepNr(inEvent->GetStepNr());
        outEvent->SetStepT(inEvent->GetStepT());

        // Clear per-event containers (must be allocated beforehand)
        vecPol->clear();
        vecPolarization->clear();
        ULambda->clear();
        UProton->clear();
        UPion->clear();

        // --- First pass: count Lambdas in the event ---
        for (Int_t i = 0; i < inEvent->GetNpa(); i++) {
            UParticle* part = inEvent->GetParticle(i);
            if (part->GetPdg() == 3122) lambdaCounter++;
        }

        // --- Second pass: loop over all particles and process Lambdas ---
        for (Int_t i = 0; i < inEvent->GetNpa(); ++i) {

            UParticle* part = inEvent->GetParticle(i);

            // If we reached the last particle and no Λ's were found,
            // we artificially create one with random position/momentum.
            if (i == inEvent->GetNpa()-1 && lambdaCounter == 0) {

                // Random position (x, y, z, t) for the new Lambda
                TLorentzVector newLambdaPos( 
                    get_random_value(0, 12), 
                    get_random_value(0, 12), 
                    get_random_value(0, 20), 
                    0
                );

                // Dummy momentum (E,px,py,pz) = (1,1,1,1)
                TLorentzVector newLambdaMom( 1., 1., 1., 1. );

                // Create artificial Lambda with PDG=3122
                UParticle* newpart = new UParticle(i, 3122, 1, 1, 1, -15, -1,
                                                   child_null, newLambdaMom, newLambdaPos, 0 );
                lambdaCounter++;
                // Keep original particle in event
                outEvent->AddParticle(*part);
                ROOT::Math::XYZVector nullPol(0, 0, 0);
                vecPolarization->push_back(nullPol);

                // Treat the newly created Lambda as current "part"
                part = newpart;
            }

            // If this is NOT a Lambda, simply copy it to output and push null polarization
            if (part->GetPdg() != 3122) {
                outEvent->AddParticle(*part);
                ROOT::Math::XYZVector nullPol(0, 0, 0);
                vecPolarization->push_back(nullPol);
                continue; 
            } 
            // If this is a Lambda (PDG 3122), process its decay
            else lambdaCounter++;

            // Copy current Lambda to local container
            lambda = *part;

            // Set kinematic / yield parameterization based on impact parameter, etc.
            set_lambda_parameterization(Lambda_yield, inEvent->GetB(), lambda); 
            Double_t fPolY = get_mean_polarization(2.87, get_centrality(inEvent->GetB()));

            // --- Lambda 4-momentum in lab frame ---
            TLorentzVector lambda_lab(lambda.Px(), lambda.Py(), lambda.Pz(), lambda.E());
            hLambdaPt->Fill(lambda_lab.Pt());

            // --- Boost to Lambda rest frame ---
            TVector3 beta = -lambda_lab.BoostVector();
            TLorentzVector lambda_rest = lambda_lab;
            lambda_rest.Boost(beta);

            // Two-body decay momentum in the rest frame (p*)
            double pStar = sqrt((mLambda*mLambda - (mProton + mPion)*(mProton + mPion)) *
                                (mLambda*mLambda - (mProton - mPion)*(mProton - mPion))) /
                                (2*mLambda);

            // Random azimuthal angle for decay
            double phi = rand->Uniform(0, 2*TMath::Pi());

            // Get polarization vector for current Lambda
            ROOT::Math::XYZVector pol = get_pol_lambda(lambda, fPolY/100., fSigmaPolVal);
            vecPol->push_back(pol);
            vecPolarization->push_back(pol);

            // Generate cos(theta) according to anisotropic distribution
            double cos_theta = get_costh(0.732, pol.R());
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

            // Numerical safety for |cosθ| slightly above 1 due to rounding
            if (TMath::Abs(cos_theta) >= 1.0) cos_theta = TMath::Sign(1.0, cos_theta);
            else sin_theta = TMath::Sqrt((1. - cos_theta) * (1. + cos_theta));

            // Define unit vector in proton direction (before rotation with respect to polarization axis)
            TVector3 unit(
                sin_theta * cos(phi),
                sin_theta * sin(phi),
                cos_theta
            );

            // Rotate the direction so that z-axis aligns with polarization vector
            TVector3 unit_tv3(unit.X(), unit.Y(), unit.Z());
            TVector3 pol_tv3(pol.X(), pol.Y(), pol.Z());
            unit_tv3.RotateUz(pol_tv3.Unit());
            unit = TVector3(unit_tv3.X(), unit_tv3.Y(), unit_tv3.Z());

            // Now recalculate angles after rotation
            phi = unit.Phi();
            cos_theta = TMath::Cos(unit.Theta());
            sin_theta = TMath::Sin(unit.Theta());

            if (TMath::Abs(cos_theta) >= 1.0) cos_theta = TMath::Sign(1.0, cos_theta);
            else sin_theta = TMath::Sqrt((1. - cos_theta) * (1. + cos_theta));

            // --- Proton momentum in Lambda rest frame ---
            TVector3 p_proton_rest = pStar * TVector3(
                sin_theta * cos(phi),
                sin_theta * sin(phi),
                cos_theta
            );

            // Pion momentum is opposite in two-body decay
            TVector3 p_pion_rest = -p_proton_rest;

            // --- 4-vectors of decay products in rest frame ---
            TLorentzVector proton_rest(p_proton_rest, sqrt(pStar*pStar + mProton*mProton));
            TLorentzVector pion_rest  (p_pion_rest,   sqrt(pStar*pStar + mPion  *mPion  ));

            // Positions of decay products are set to the parent Lambda position in lab frame
            TLorentzVector proton_lab_pos = lambda.GetPosition();
            TLorentzVector pion_lab_pos   = lambda.GetPosition();

            Double_t fWeight = 0;         // Optional event weight (unused here)
            Int_t enhancedFlag = 0;       // Flag for enhanced statistics (0 = original Lambda)

            // Create UParticle objects for decay products
            proton = UParticle(i, 2212, 0, i, i, enhancedFlag, -1,
                               child_null, proton_rest, proton_lab_pos, fWeight);
            pion   = UParticle(i, -211, 0, i, i, enhancedFlag, -1,
                               child_null, pion_rest,   pion_lab_pos,   fWeight);

            // Add original Lambda to output event (not replaced)
            outEvent->AddParticle(lambda);

            // Store Lambda and its decay products into auxiliary vectors
            ULambda->push_back(lambda);
            UProton->push_back(proton);
            UPion->push_back(pion);

            // Fill QA histograms
            hCosTheta->Fill(cos_theta);
            hAngVsPt->Fill(lambda_lab.Pt(), cos_theta);
            hLambdaPhi->Fill(get_positive_phi(lambda.GetMomentum().Phi()));
        }

        // Total number of particles in original event
        Int_t nParticles = inEvent->GetNpa();

        // --- Artificial enhancement of Lambda statistics ---
        Double_t fEnhanceStat = enhanceStat;
        while (fEnhanceStat > 0) { // Repeat until enhancement counter is exhausted
            Int_t enhancedFlag = -9; // Mark these as "enhanced" Lambdas

            // Dummy random momentum (not really used; overwritten in parameterization)
            TLorentzVector mom_rand( 1., 1., 1., 1. );

            // Randomized position for enhanced Lambda inside some small volume
            TLorentzVector pos_rand(
                get_random_value(0, 0.03),
                get_random_value(0, 0.03),
                get_random_value(0, 0.03),
                get_random_value(0, 0.03)
            );

            // Temporary enhanced Lambda
            UParticle enhancedLambda(1, 1, 1, 1, 1, 1, 1, dummy,
                                     1., 1., 1., 1., 1., 1., 1., 1., 1.);

            enhancedLambda.SetPosition(pos_rand);
            enhancedLambda.SetMate(enhancedFlag); // Use Mate to tag enhanced particles

            // Sample Lambda kinematics from parameterization/yield
            set_lambda_parameterization(Lambda_yield, inEvent->GetB(), enhancedLambda); 
            Double_t fPolY = get_mean_polarization(2.87, get_centrality(inEvent->GetB()));

            // 4-momentum in lab frame
            TLorentzVector lambda_lab = enhancedLambda.GetMomentum();
            hLambdaPt->Fill(lambda_lab.Pt());

            TVector3 beta = -lambda_lab.BoostVector();
            // p* for decay (computed again; could be factored out)
            double pStar = sqrt((mLambda*mLambda - (mProton + mPion)*(mProton + mPion)) *
                                (mLambda*mLambda - (mProton - mPion)*(mProton - mPion))) /
                                (2*mLambda);

            double phi = rand->Uniform(0, 2*TMath::Pi());

            // Polarization for enhanced Lambda
            ROOT::Math::XYZVector pol = get_pol_lambda(enhancedLambda, fPolY/100., fSigmaPolVal);
            vecPol->push_back(pol);
            vecPolarization->push_back(pol);

            double cos_theta = get_costh(0.732, pol.R());
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

            if (TMath::Abs(cos_theta) >= 1.0) cos_theta = TMath::Sign(1.0, cos_theta);
            else sin_theta = TMath::Sqrt((1. - cos_theta) * (1. + cos_theta));

            TVector3 unit(
                sin_theta * cos(phi),
                sin_theta * sin(phi),
                cos_theta
            );

            TVector3 unit_tv3(unit.X(), unit.Y(), unit.Z());
            TVector3 pol_tv3(pol.X(), pol.Y(), pol.Z());
            unit_tv3.RotateUz(pol_tv3.Unit());
            unit = TVector3(unit_tv3.X(), unit_tv3.Y(), unit_tv3.Z());

            phi = unit.Phi();
            cos_theta = TMath::Cos(unit.Theta());
            sin_theta = TMath::Sin(unit.Theta());

            if (TMath::Abs(cos_theta) >= 1.0) cos_theta = TMath::Sign(1.0, cos_theta);
            else sin_theta = TMath::Sqrt((1. - cos_theta) * (1. + cos_theta));

            // Proton momentum in Lambda rest frame
            TVector3 p_proton_rest = pStar * TVector3(
                sin_theta * cos(phi),
                sin_theta * sin(phi),
                cos_theta
            );

            // Pion momentum is opposite
            TVector3 p_pion_rest = -p_proton_rest;

            // 4-vectors of decay products in rest frame
            TLorentzVector proton_rest(p_proton_rest, sqrt(pStar*pStar + mProton*mProton));
            TLorentzVector pion_rest  (p_pion_rest,   sqrt(pStar*pStar + mPion  *mPion  ));

            // Index for newly added particles in event
            Int_t indexOfEnhancedParticle = nParticles++;

            // Enhanced decay products (tagged with enhancedFlag)
            UParticle protonEnhanced(
                indexOfEnhancedParticle, 2212, 0, indexOfEnhancedParticle, indexOfEnhancedParticle, enhancedFlag, -1,
                child_null, proton_rest, enhancedLambda.GetPosition(), 0
            );
            UParticle pionEnhanced(
                indexOfEnhancedParticle, -211, 0, indexOfEnhancedParticle, indexOfEnhancedParticle, enhancedFlag, -1,
                child_null, pion_rest, enhancedLambda.GetPosition(), 0
            );

            // Fill QA histograms for enhanced sample
            hCosTheta->Fill(cos_theta);
            hAngVsPt->Fill(lambda_lab.Pt(), cos_theta);
            hLambdaPhi->Fill(get_positive_phi(enhancedLambda.GetMomentum().Phi()));

            // Add enhanced Lambda and its products to output event and vectors
            outEvent->AddParticle(enhancedLambda);

            ULambda->push_back(enhancedLambda);
            UProton->push_back(protonEnhanced);
            UPion->push_back(pionEnhanced);

            // Decrease enhancement counter
            fEnhanceStat--;
        }

        // Store current event into the output tree
        outTree->Fill();
    }

    // --- Write tree to file (overwrite "decays" if exists) ---
    outFile->cd();
    outTree->Write("",TObject::kOverwrite);

    // Optionally, histograms could also be written to the file
    // hLambdaPt->Write("",TObject::kOverwrite);
    // hCosTheta->Write("",TObject::kOverwrite);
    // hLambdaPhi->Write("",TObject::kOverwrite);
    // hAngVsPt->Write("",TObject::kOverwrite);

    // --- Simple QA canvas with four plots ---
    TCanvas *c1 = new TCanvas("c1", "Lambda_Decays", 1200, 800);
    c1->Divide(2,2);
    c1->cd(1); hLambdaPt->Draw();
    c1->cd(2); hCosTheta->Draw();
    c1->cd(3); hAngVsPt->Draw("colz");
    c1->cd(4); hLambdaPhi->Draw();

    c1->SaveAs("lambda_decays_plots.png");

    // --- Cleanup ---
    delete c1;
    delete outEvent;
    delete rand;
    delete inChain;

    outFile->Close();
    Lambda_yield->Close();
}

