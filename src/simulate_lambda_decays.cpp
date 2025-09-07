
#include "simulate_lambda_decays.hpp"


void simulate_lambda_decays(TString inputFile, TString outputFile, TString confInFile, Int_t enhanceStat) {

    TChain *inChain = nullptr;
    TFile *outFile = nullptr;
    TTree *outTree = nullptr;
    
    // Important variables that were accidentally removed
    double cos_theta_p = 0;  // Proton angle relative to some axis (for anisotropy)
    Int_t eventID;
    std::vector<UParticle>* ULambda = nullptr;
    std::vector<UParticle>* UProton = nullptr;
    std::vector<UParticle>* UPion = nullptr;
    std::vector<ROOT::Math::XYZVector>* vecPolarization = nullptr;
    std::vector<ROOT::Math::XYZVector>* vecPol = nullptr;


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
        outTree->Branch("ULambda", &ULambda);
        outTree->Branch("UProton", &UProton);
        outTree->Branch("UPion", &UPion);
        outTree->Branch("LambdaPolarization", &vecPol);
        outTree->Branch("Polarization", &vecPolarization);//switch with vecPol


    } else {

        std::cout << "Appending to existing tree with " << outTree->GetEntries() << " entries" << std::endl;

        // Connect branches for existing tree
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


    TFile* Lambda_yield = TFile::Open(confInFile,"READ");

    // Histograms
    TH1D *hLambdaPt  = new TH1D("hLambdaPt", "Lambda pT;pT [GeV/c];Counts", 100, 0, 2);
    TH1D *hCosTheta  = new TH1D("hCosTheta", "Proton cos(#theta);cos(#theta);Counts", 100, -1, 1);
    TH1D *hLambdaPhi = new TH1D("hLambdaPhi", "Lambda (#Delta#varphi);#Delta#varphi;Counts", 100, 0, 2.*TMath::Pi());
    TH2D *hAngVsPt   = new TH2D("hAngVsPt", "cos(#theta) vs Lambda pT;pT [GeV/c];cos(#theta)", 50, 0, 2, 50, -1, 1);

    // Particle instances
    Int_t dummy[2];
    UParticle lambda(1, 1, 1, 1, 1, 1, 1, dummy, 1., 1., 1., 1., 1., 1., 1., 1., 1.);
    UParticle proton(1, 1, 1, 1, 1, 1, 1, dummy, 1., 1., 1., 1., 1., 1., 1., 1., 1.);
    UParticle pion(1, 1, 1, 1, 1, 1, 1, dummy, 1., 1., 1., 1., 1., 1., 1., 1., 1.);

    // Physics parameters
    TRandom3* rand = new TRandom3(0);  
    const double mLambda = 1.115683;
    const double mProton = 0.938272;
    const double mPion = 0.139570;
    Double_t fSigmaPolVal = 1.5;
    const double anisotropy = 0.732;
    const double polPercent = 0.6;
    Double_t fBMin = 3.44;
    Double_t fBMax = 7.44;
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

        vecPol->clear();
        vecPolarization->clear();
        ULambda->clear();
        UProton->clear();
        UPion->clear();
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
            
            // Boost to Lambda rest frame
            TVector3 beta = -lambda_lab.BoostVector();
            TLorentzVector lambda_rest = lambda_lab;
            lambda_rest.Boost(beta);
            
            // Generate decay in rest frame with anisotropy
            double pStar = sqrt((mLambda*mLambda - (mProton + mPion)*(mProton + mPion)) *
                                (mLambda*mLambda - (mProton - mPion)*(mProton - mPion))) /
                                (2*mLambda);
            
            double phi = rand->Uniform(0, 2*TMath::Pi());
            ROOT::Math::XYZVector pol = get_pol_lambda(lambda, fPolY/100., fSigmaPolVal);
            vecPol->push_back(pol);
            vecPolarization->push_back(pol);

            double cos_theta = get_costh(0.732, pol.R());
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

            if (TMath::Abs(cos_theta) >= 1.0) cos_theta = TMath::Sign(1.0, cos_theta);
            else sin_theta = TMath::Sqrt((1. - cos_theta) * (1. + cos_theta));

            TVector3 unit = TVector3(
                sin_theta * cos(phi),
                sin_theta * sin(phi),
                cos_theta
            );
            TVector3 unit_tv3(unit.X(), unit.Y(), unit.Z());
            TVector3 pol_tv3(pol.X(), pol.Y(), pol.Z());
            unit_tv3.RotateUz(pol_tv3.Unit());
            unit = TVector3(unit_tv3.X(), unit_tv3.Y(), unit_tv3.Z());

            // unit.RotateUz(pol.Unit()); // rotate Z to norm (with extra phi rotation which is random)

            phi = unit.Phi();
            cos_theta = TMath::Cos(unit.Theta());
            sin_theta = TMath::Sin(unit.Theta());

            if (TMath::Abs(cos_theta) >= 1.0) cos_theta = TMath::Sign(1.0, cos_theta);
            else sin_theta = TMath::Sqrt((1. - cos_theta) * (1. + cos_theta));

            //Generate proton momentum in lambda rest frame
            TVector3 p_proton_rest = pStar * TVector3(
                sin_theta * cos(phi),
                sin_theta * sin(phi),
                cos_theta
            );
            
            // Pion momentum (opposite to proton)
            TVector3 p_pion_rest = -p_proton_rest;

            // Create 4-vectors in rest frame
            TLorentzVector proton_rest(p_proton_rest, sqrt(pStar*pStar + mProton*mProton));
            TLorentzVector pion_rest(p_pion_rest, sqrt(pStar*pStar + mPion*mPion));

            TLorentzVector proton_lab_pos = lambda.GetPosition();
            TLorentzVector pion_lab_pos   = lambda.GetPosition();
            
            Double_t fWeight = 0;
            Int_t enhancedFlag = 0;

            proton = UParticle(i, 2212, 0, i, i, enhancedFlag, -1, child_null, proton_rest, proton_lab_pos, fWeight);
            pion   = UParticle(i, -211, 0, i, i, enhancedFlag, -1, child_null, pion_rest, pion_lab_pos, fWeight);

            outEvent->AddParticle(lambda);

            ULambda->push_back(lambda);
            UProton->push_back(proton);
            UPion->push_back(pion);

            //std::cout << "Lambda counter in event #" << iEvent << " b4 enhance : " << ULambda->size() << " lambdas\n"; 

            // Fill histograms
            hCosTheta->Fill(cos_theta);
            hAngVsPt->Fill(lambda_lab.Pt(), cos_theta);
            hLambdaPhi->Fill(get_positive_phi(lambda.GetMomentum().Phi()));


            Double_t fEnhanceStat = enhanceStat;
            while(fEnhanceStat > 1){ //enhancing of lambdas
                enhancedFlag = -9;

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
                
                TVector3 beta = -lambda_lab.BoostVector();
                double pStar = sqrt((mLambda*mLambda - (mProton + mPion)*(mProton + mPion)) *
                                    (mLambda*mLambda - (mProton - mPion)*(mProton - mPion))) /
                                    (2*mLambda);

                double phi = rand->Uniform(0, 2*TMath::Pi());

                ROOT::Math::XYZVector pol = get_pol_lambda(enhancedLambda, fPolY/100., fSigmaPolVal);
                vecPol->push_back(pol);
                vecPolarization->push_back(pol);

                double cos_theta = get_costh(0.732, pol.R());
                double sin_theta = sqrt(1.0 - cos_theta*cos_theta);


                if (TMath::Abs(cos_theta) >= 1.0) cos_theta = TMath::Sign(1.0, cos_theta);
                else sin_theta = TMath::Sqrt((1. - cos_theta) * (1. + cos_theta));

                TVector3 unit = TVector3(
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


                TVector3 p_proton_rest = pStar * TVector3(
                    sin_theta * cos(phi),
                    sin_theta * sin(phi),
                    cos_theta
                );
                
                TVector3 p_pion_rest = -p_proton_rest;

                TLorentzVector proton_rest(p_proton_rest, sqrt(pStar*pStar + mProton*mProton));
                TLorentzVector pion_rest(p_pion_rest, sqrt(pStar*pStar + mPion*mPion));

                UParticle protonEnhanced(
                    i, 2212, 0, i, i, enhancedFlag, -1,
                    child_null, proton_rest, enhancedLambda.GetPosition(), 0
                );
                UParticle pionEnhanced(
                    i, -211, 0, i, i, enhancedFlag, -1,
                    child_null, pion_rest, enhancedLambda.GetPosition(), 0
                );

                hCosTheta->Fill(cos_theta);
                hAngVsPt->Fill(lambda_lab.Pt(), cos_theta);
                hLambdaPhi->Fill(get_positive_phi(enhancedLambda.GetMomentum().Phi()));

                outEvent->AddParticle(enhancedLambda);

                ULambda->push_back(enhancedLambda);
                UProton->push_back(protonEnhanced);
                UPion->push_back(pionEnhanced);
    
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
