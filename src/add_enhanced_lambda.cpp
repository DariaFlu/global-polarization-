
#include "add_enhanced_lambda.hpp"


void add_enhanced_lambda(TString inputFile, TString outputFile, TString confInFile, Int_t enhanceStat) {

    TChain *inChain = nullptr;
    TFile *outFile = nullptr;
    TTree *outTree = nullptr;
    
    // Important variables that were accidentally removed
    double cos_theta_p = 0;  // Proton angle relative to some axis (for anisotropy)
    Int_t eventID;
    std::vector<ROOT::Math::XYZVector>* vecPolarization = nullptr;


    // Open input file as TChain
    inChain = new TChain("events");
    inChain->Add(inputFile);

    // Create output file and tree
    // TFile *outFile = new TFile(outputFile, "RECREATE");

    outFile = TFile::Open(outputFile, "UPDATE");
    if (!outFile || outFile->IsZombie()) {
        delete outFile;
        outFile = TFile::Open(outputFile, "RECREATE");
        std::cout << "Creating new output file: " << outputFile << std::endl;

    }

    // --- Copy URun object from first input file ---
    // Извлекаем имя первого файла из TChain
    TObjArray* fileElements = inChain->GetListOfFiles();
    if (fileElements && fileElements->GetEntries() > 0) {
        TChainElement* chEl = (TChainElement*)fileElements->At(0);
        TString firstFileName = chEl->GetTitle();

        TFile* fFirst = TFile::Open(firstFileName, "READ");
        if (fFirst && !fFirst->IsZombie()) {
            URun* inRun = nullptr;
            fFirst->GetObject("run", inRun);
            if (inRun) {
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


    // Event variables
    UEvent  *inEvent = nullptr;
    UEvent  *outEvent = new UEvent();
    inChain->SetBranchAddress("event", &inEvent);

    outTree = (TTree*)outFile->Get("decays");
    bool newTree = false;
    if (!outTree) {

        outTree = new TTree("decays", "Lambda decay products");
        newTree = true;
        std::cout << "Creating new output tree" << std::endl;
        
        // Setup branches for new tree
        outTree->Branch("event", &outEvent);
        // outTree->Branch("run", &outRun);

        outTree->Branch("cos_theta_p", &cos_theta_p);
        outTree->Branch("eventID", &eventID, "eventID/I");
        outTree->Branch("Polarization", &vecPolarization);//switch with vecPol


    } else {

        std::cout << "Appending to existing tree with " << outTree->GetEntries() << " entries" << std::endl;

        // Connect branches for existing tree
        outTree->SetBranchAddress("event", &outEvent);
        // outTree->SetBranchAddress("run", &outRun);

        outTree->SetBranchAddress("cos_theta_p", &cos_theta_p);
        outTree->SetBranchAddress("eventID", &eventID);
        outTree->SetBranchAddress("Polarization", &vecPolarization);

    }


    TFile* Lambda_yield = TFile::Open(confInFile,"READ");

    // Histograms
    TH1D *hLambdaPt  = new TH1D("hLambdaPt", "Lambda pT;pT [GeV/c];Counts", 100, 0, 2);
    TH1D *hCosTheta  = new TH1D("hCosTheta", "Proton cos(#theta);cos(#theta);Counts", 100, -1, 1);
    TH1D *hLambdaPhi = new TH1D("hLambdaPhi", "Lambda (#Delta#varphi);#Delta#varphi;Counts", 100, 0, 2.*TMath::Pi());
    TH2D *hAngVsPt   = new TH2D("hAngVsPt", "cos(#theta) vs Lambda pT;pT [GeV/c];cos(#theta)", 50, 0, 2, 50, -1, 1);

    // Particle instances
    Int_t dummy[2];
    UParticle lambda(1, 1, 1, 1, 1, 1, 1, dummy, 1., 1., 1., 1., 1., 1., 1., 1., 1.);

    // Physics parameters
    TRandom3* rand = new TRandom3(0);  
    const double mLambda = 1.115683;    
    Double_t fSigmaPolVal = 1.5;
    const double anisotropy = 0.732;
    Int_t child_null[2] = {0, 0};


    // Process events
    Long64_t nEvents = inChain->GetEntries();

    std::cout << "Events: " << nEvents << std::endl;
    for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {

        inChain->GetEntry(iEvent);
        outEvent->Clear(); // Clear previous event
        Int_t lambdaCounter = 0;

        // Copy event header information
        outEvent->SetEventNr(inEvent->GetEventNr());
        outEvent->SetB(inEvent->GetB());
        outEvent->SetPhi(inEvent->GetPhi());
        outEvent->SetNes(inEvent->GetNes());
        outEvent->SetStepNr(inEvent->GetStepNr());
        outEvent->SetStepT(inEvent->GetStepT());

        vecPolarization->clear();
        //first we COUNT lambdas
        for (Int_t i = 0; i < inEvent->GetNpa(); i++) {
            UParticle* part = inEvent->GetParticle(i);
            if (part->GetPdg() == 3122) lambdaCounter++;
        }

        for (Int_t i = 0; i < inEvent->GetNpa(); ++i) {

            UParticle* part = inEvent->GetParticle(i);
            if(i == inEvent->GetNpa()-1 && lambdaCounter == 0){

                TLorentzVector newLambdaPos( 
                                            get_random_value(0, 12), 
                                            get_random_value(0, 12), 
                                            get_random_value(0, 20), 
                                            0
                                            );
                TLorentzVector newLambdaMom( 1., 1., 1., 1. );

                UParticle* newpart = new UParticle(i, 3122, 1, 1, 1, -15, -1,
                                    child_null, newLambdaMom, newLambdaPos, 0 );
                lambdaCounter++;
                outEvent->AddParticle(*part);
                ROOT::Math::XYZVector nullPol(0, 0, 0);
                vecPolarization->push_back(nullPol);
                part = newpart;
            }

            if (part->GetPdg() != 3122) {
                outEvent->AddParticle(*part);
                ROOT::Math::XYZVector nullPol(0, 0, 0);
                vecPolarization->push_back(nullPol);
                continue; 
            }// Select Lambdas (PDG code 3122)
            else lambdaCounter++;

            
            lambda = *part;
            set_lambda_parameterization(Lambda_yield, inEvent->GetB(), lambda); 
            Double_t fPolY = get_mean_polarization(2.87, get_centrality(inEvent->GetB()));

            // Get Lambda 4-momentum
            TLorentzVector lambda_lab(lambda.Px(), lambda.Py(), lambda.Pz(), lambda.E());
            hLambdaPt->Fill(lambda_lab.Pt());

            ROOT::Math::XYZVector pol = get_pol_lambda(lambda, fPolY/100., fSigmaPolVal);
            vecPolarization->push_back(pol);

            outEvent->AddParticle(lambda);

            //std::cout << "Lambda counter in event #" << iEvent << " b4 enhance : " << ULambda->size() << " lambdas\n"; 
            // Fill histograms
            hLambdaPhi->Fill(get_positive_phi(lambda.GetMomentum().Phi()));

            Double_t fEnhanceStat = enhanceStat;
            while(fEnhanceStat > 1){ //enhancing of lambdas
                Int_t enhancedFlag = -9;

                TLorentzVector mom_rand( 1., 1., 1., 1. );
                
                TLorentzVector pos_rand(
                    get_random_value(lambda.X(), 0.03),
                    get_random_value(lambda.Y(), 0.03),
                    get_random_value(lambda.Z(), 0.03),
                    get_random_value(lambda.T(), 0.03)
                );

                UParticle enhancedLambda(lambda);
                enhancedLambda.SetPosition(pos_rand);
                enhancedLambda.SetMate(enhancedFlag);

                set_lambda_parameterization(Lambda_yield, inEvent->GetB(), enhancedLambda); 
                fPolY = get_mean_polarization(2.87, get_centrality(inEvent->GetB()));
                
                TLorentzVector lambda_lab = enhancedLambda.GetMomentum();
                hLambdaPt->Fill(lambda_lab.Pt());
                
                ROOT::Math::XYZVector pol = get_pol_lambda(enhancedLambda, fPolY/100., fSigmaPolVal);
                vecPolarization->push_back(pol);

                hLambdaPhi->Fill(get_positive_phi(enhancedLambda.GetMomentum().Phi()));

                outEvent->AddParticle(enhancedLambda);
                //std::cout << "ooooo Lambda counter in event #" << iEvent << " inside enhance : " << ULambda->size() << " lambdas\n"; 
                fEnhanceStat--;
            }
        }
        //std::cout << "Lambda counter in event #" << iEvent << " after add & enhance : " << ULambda->size() << " lambdas\n"; 

        outTree->Fill();
    }

    outFile->cd();
    outTree->Write("",TObject::kOverwrite);

    // hLambdaPt->Write("",TObject::kOverwrite);
    // hCosTheta->Write("",TObject::kOverwrite);
    // hLambdaPhi->Write("",TObject::kOverwrite);
    // hAngVsPt->Write("",TObject::kOverwrite);

    // Plot results
    TCanvas *c1 = new TCanvas("c1", "Lambda_Decays", 1200, 800);
    c1->Divide(2,2);
    c1->cd(1); hLambdaPt->Draw();
    c1->cd(2); hCosTheta->Draw();
    c1->cd(3); hAngVsPt->Draw("colz");
    c1->cd(4); hLambdaPhi->Draw();

    c1->SaveAs("lambda_decays_plots.png");

    // Cleanup
    delete c1;
    delete outEvent;
    delete rand;
    delete inChain;

    outFile->Close();
    Lambda_yield->Close();
}
