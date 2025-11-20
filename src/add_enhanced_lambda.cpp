
#include "add_enhanced_lambda.hpp"

// Main function to enhance Λ production and write into an output tree.
// - inputFile  : UniGen input file(s) with "events" TTree
// - outputFile : ROOT file where "decays" tree is stored/updated
// - confInFile : ROOT file with parameterization / yields for Lambda
// - enhanceStat: factor to artificially enhance Lambda statistics
void add_enhanced_lambda(TString inputFile, TString outputFile, TString confInFile, Int_t enhanceStat) {

    // Input chain and output file/tree
    TChain *inChain = nullptr;
    TFile *outFile = nullptr;
    TTree *outTree = nullptr;
    
    // Variables stored in the output tree
    double cos_theta_p = 0;  // Proton angle relative to some axis (for anisotropy / polarization analysis)
    Int_t eventID;           // Event ID (to be filled if needed)
    // Vector that will be written to the TTree as branches
    std::vector<ROOT::Math::XYZVector>* vecPolarization = nullptr;

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
        
        // Setup branches for new tree
        outTree->Branch("event", &outEvent);

        outTree->Branch("cos_theta_p", &cos_theta_p);
        outTree->Branch("eventID", &eventID, "eventID/I");
        outTree->Branch("Polarization", &vecPolarization);//switch with vecPol


    } else {
        // Existing tree found: append new entries
        std::cout << "Appending to existing tree with " << outTree->GetEntries() << " entries" << std::endl;

        // Connect branches for existing tree
        outTree->SetBranchAddress("event", &outEvent);

        outTree->SetBranchAddress("cos_theta_p", &cos_theta_p);
        outTree->SetBranchAddress("eventID", &eventID);
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
    // Temporary particle object used as containers for Lambda 
    UParticle lambda(1, 1, 1, 1, 1, 1, 1, dummy,
                    1., 1., 1., 1., 1., 1., 1., 1., 1.);

    // --- Physics constants and parameters ---
    TRandom3* rand = new TRandom3(0);  // Random number generator with fixed seed
    const double mLambda = 1.115683;   // Λ mass [GeV/c^2]
    const double mProton = 0.938272;   // Proton mass [GeV/c^2]
    const double mPion   = 0.139570;   // Pion mass [GeV/c^2]
    Double_t fSigmaPolVal = 1.5;       // Width for polarization smearing (model parameter)
    const double anisotropy = 0.732;   // Decay anisotropy parameter (αΛ)
    Int_t child_null[2] = {0, 0};      // Children indices placeholder for UParticle

    // Process events
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

        // Clear per-event container
        vecPolarization->clear();

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
            if(i == inEvent->GetNpa()-1 && lambdaCounter == 0){
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

            // Get polarization vector for current Lambda
            ROOT::Math::XYZVector pol = get_pol_lambda(lambda, fPolY/100., fSigmaPolVal);
            vecPolarization->push_back(pol);

            // Add original Lambda to output event (not replaced)
            outEvent->AddParticle(lambda);

            // Fill histograms
            hLambdaPhi->Fill(get_positive_phi(lambda.GetMomentum().Phi()));

            // --- Artificial enhancement of Lambda statistics ---
            Double_t fEnhanceStat = enhanceStat;
            while(fEnhanceStat > 1){ // Repeat until enhancement counter is exhausted
                Int_t enhancedFlag = -9; // Mark these as "enhanced" Lambdas

                // Dummy random momentum (not really used; overwritten in parameterization)
                TLorentzVector mom_rand( 1., 1., 1., 1. );
                
                // Randomized position for enhanced Lambda inside some small volume
                TLorentzVector pos_rand(
                    get_random_value(lambda.X(), 0.03),
                    get_random_value(lambda.Y(), 0.03),
                    get_random_value(lambda.Z(), 0.03),
                    get_random_value(lambda.T(), 0.03)
                );

                // Temporary enhanced Lambda
                UParticle enhancedLambda(lambda);
                enhancedLambda.SetPosition(pos_rand);
                enhancedLambda.SetMate(enhancedFlag);

                // Sample Lambda kinematics from parameterization/yield
                set_lambda_parameterization(Lambda_yield, inEvent->GetB(), enhancedLambda); 
                fPolY = get_mean_polarization(2.87, get_centrality(inEvent->GetB()));
                
                // 4-momentum in lab frame
                TLorentzVector lambda_lab = enhancedLambda.GetMomentum();
                hLambdaPt->Fill(lambda_lab.Pt());
                
                // Polarization for enhanced Lambda
                ROOT::Math::XYZVector pol = get_pol_lambda(enhancedLambda, fPolY/100., fSigmaPolVal);
                vecPolarization->push_back(pol);

                // Fill QA histograms for enhanced sample
                hLambdaPhi->Fill(get_positive_phi(enhancedLambda.GetMomentum().Phi()));

                // Add enhanced Lambda and its products to output event and vectors
                outEvent->AddParticle(enhancedLambda);

                // Decrease enhancement counter
                fEnhanceStat--;
            }
        }
        // Store current event into the output tree
        outTree->Fill();
    }
    // --- Write tree to file (overwrite "decays" if exists) ---
    outFile->cd();
    outTree->Write("",TObject::kOverwrite);

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
