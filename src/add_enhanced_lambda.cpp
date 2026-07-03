#include "add_enhanced_lambda.hpp"

// Main function to enhance Λ production and write into an output tree.
// - outputFile : ROOT file where "decays" tree is stored/updated
// - qaFile     : ROOT file with parameterization / yields for Lambda
// - enhanceStat: factor to artificially enhance Lambda statistics
void add_enhanced_lambda(const TString &outputFile,
                         Long64_t fileStart,
                         Long64_t fileCount,
                         const TString &pattern,
                         const TString &qaFile,
                         Int_t enhanceStat,
                         UInt_t seed)
{
    // --- Seed logic ---
    if (seed == 0)
    {
        const char *a = gSystem->Getenv("SLURM_ARRAY_TASK_ID");
        if (a)
            seed = 12345u + (UInt_t)atoi(a);
        else
            seed = 12345u;
    }
    std::cout << "RNG seed = " << seed << std::endl;

    // Input chain and output file/tree
    TChain *inChain = nullptr;
    TTree *outTree = nullptr;

    // Variables stored in the output tree
    double cos_theta_p = 0;
    Int_t eventID;
    std::vector<ROOT::Math::XYZVector> *vecPolarization = nullptr;

    inChain = new TChain("events");
    TString strPattern(pattern);

    if (strPattern.Contains("%"))
    {
        for (long long i = fileStart; i < fileStart + fileCount; ++i)
        {
            TString filename = TString::Format(pattern, i, i);
            inChain->Add(filename);
            std::cout << "Added to chain: " << filename << std::endl;
        }
    }
    else
    {
        inChain->Add(pattern);
        std::cout << "Added wildcard pattern to chain: " << pattern << std::endl;
    }

    if (inChain->GetListOfFiles()->GetEntries() == 0)
    {
        std::cerr << "CRITICAL ERROR: No files were added to the TChain. Check your pattern: " << pattern << std::endl;
        return;
    }

    TFile *outFile = TFile::Open(outputFile, "UPDATE");
    if (!outFile || outFile->IsZombie())
    {
        delete outFile;
        outFile = TFile::Open(outputFile, "RECREATE");
        std::cout << "Creating new output file: " << outputFile << std::endl;
    }

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

    UEvent *inEvent = nullptr;
    UEvent *outEvent = new UEvent();
    inChain->SetBranchAddress("event", &inEvent);

    outFile->cd();

    outTree = (TTree *)outFile->Get("decays");

    bool newTree = false;
    if (!outTree)
    {
        outTree = new TTree("decays", "Lambda decay products");
        newTree = true;
        std::cout << "Creating new output tree" << std::endl;

        outTree->Branch("event", &outEvent);
        outTree->Branch("cos_theta_p", &cos_theta_p);
        outTree->Branch("eventID", &eventID, "eventID/I");
        outTree->Branch("Polarization", &vecPolarization);
    }
    else
    {
        std::cout << "Appending to existing tree with " << outTree->GetEntries() << " entries" << std::endl;

        outTree->SetBranchAddress("event", &outEvent);
        outTree->SetBranchAddress("cos_theta_p", &cos_theta_p);
        outTree->SetBranchAddress("eventID", &eventID);
        outTree->SetBranchAddress("Polarization", &vecPolarization);
    }

    TFile *Lambda_yield = TFile::Open(qaFile, "READ");

    TH1D *hLambdaPt = new TH1D("hLambdaPt", "Lambda pT;pT [GeV/c];Counts", 100, 0, 2);
    TH1D *hCosTheta = new TH1D("hCosTheta", "Proton cos(#theta);cos(#theta);Counts", 100, -1, 1);
    TH1D *hLambdaPhi = new TH1D("hLambdaPhi", "Lambda (#Delta#varphi);#Delta#varphi;Counts", 100, 0, 2. * TMath::Pi());
    TH2D *hAngVsPt = new TH2D("hAngVsPt", "cos(#theta) vs Lambda pT;pT [GeV/c];cos(#theta)", 50, 0, 2, 50, -1, 1);

    Int_t dummy[2];
    UParticle lambda(1, 1, 1, 1, 1, 1, 1, dummy,
                     1., 1., 1., 1., 1., 1., 1., 1., 1.);

    TRandom3 rand(seed);

    const double mLambda = 1.115683;
    const double mProton = 0.938272;
    const double mPion = 0.139570;
    Double_t fSigmaPolVal = 1.5;
    const double anisotropy = 0.732;
    Int_t child_null[2] = {0, 0};

    Long64_t nEvents = inChain->GetEntries();
    std::cout << "Events: " << nEvents << std::endl;

    for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++)
    {
        if (iEvent % 200 == 0)
        {
            std::cout << "Event No " << iEvent << " / " << nEvents << std::endl;
        }

        inChain->GetEntry(iEvent);
        outEvent->Clear();

        outEvent->SetEventNr(inEvent->GetEventNr());
        outEvent->SetB(inEvent->GetB());
        outEvent->SetPhi(inEvent->GetPhi());
        outEvent->SetNes(inEvent->GetNes());
        outEvent->SetStepNr(inEvent->GetStepNr());
        outEvent->SetStepT(inEvent->GetStepT());

        vecPolarization->clear();

        // --- First pass: count Lambdas in the event ---
        Int_t lambdaCounter = 0;
        for (Int_t i = 0; i < inEvent->GetNpa(); i++)
        {
            UParticle *part = inEvent->GetParticle(i);
            if (part->GetPdg() == 3122)
                lambdaCounter++;
        }

        // --- Second pass: loop over all particles and process Lambdas ---
        for (Int_t i = 0; i < inEvent->GetNpa(); ++i)
        {
            UParticle *part = inEvent->GetParticle(i);
            if (!part)
                continue;

            // If this is NOT a lambda, simply copy it to output
            if (part->GetPdg() != 3122)
            {
                outEvent->AddParticle(*part);
                ROOT::Math::XYZVector nullPol(0, 0, 0);
                vecPolarization->push_back(nullPol);
                continue;
            }
            lambda = *part;
            set_lambda_parameterization(Lambda_yield, inEvent->GetB(), lambda, rand);
            Double_t fPolY = get_mean_polarization(2.87, get_centrality(inEvent->GetB()));

            TLorentzVector lambda_lab(lambda.Px(), lambda.Py(), lambda.Pz(), lambda.E());
            hLambdaPt->Fill(lambda_lab.Pt());

            ROOT::Math::XYZVector pol = get_pol_lambda(lambda, fPolY / 100., fSigmaPolVal);
            vecPolarization->push_back(pol);

            outEvent->AddParticle(lambda);
            hLambdaPhi->Fill(get_positive_phi(lambda.GetMomentum().Phi()));
        }

        // --- Artificial Lambda creation logic ---
        if (lambdaCounter == 0)
        {
            TLorentzVector newLambdaPos(
                get_random_value_rand(rand, 0, 12),
                get_random_value_rand(rand, 0, 12),
                get_random_value_rand(rand, 0, 20),
                0);

            TLorentzVector newLambdaMom(1., 1., 1., 1.);

            UParticle newLambda(outEvent->GetNpa(), 3122, 1, 1, 1, -15, -1,
                                child_null, newLambdaMom, newLambdaPos, 0);

            set_lambda_parameterization(Lambda_yield, inEvent->GetB(), newLambda, rand);

            outEvent->AddParticle(newLambda);

            Double_t fPolY = get_mean_polarization(2.87, get_centrality(inEvent->GetB()));
            ROOT::Math::XYZVector pol = get_pol_lambda(newLambda, fPolY / 100., fSigmaPolVal);
            vecPolarization->push_back(pol);

            TLorentzVector lambda_lab = newLambda.GetMomentum();
            hLambdaPt->Fill(lambda_lab.Pt());
            hLambdaPhi->Fill(get_positive_phi(lambda_lab.Phi()));
        }

        Int_t nParticles = inEvent->GetNpa();

        // --- Artificial enhancement of Lambda statistics ---
        Double_t fEnhanceStat = enhanceStat;
        while (fEnhanceStat > 0)
        {
            Int_t enhancedFlag = -9;

            TLorentzVector mom_rand(1., 1., 1., 1.);

            TLorentzVector pos_rand(
                get_random_value_rand(rand, 0, 0.03),
                get_random_value_rand(rand, 0, 0.03),
                get_random_value_rand(rand, 0, 0.03),
                get_random_value_rand(rand, 0, 0.03));
            Int_t indexOfEnhancedParticle = nParticles++;

            UParticle enhancedLambda(indexOfEnhancedParticle, 3122, 1, 1, 1, enhancedFlag,
                                     -1, child_null,
                                     1., 1., 1., 1., 1., 1., 1., 1., 1.);
            enhancedLambda.SetPosition(pos_rand);

            set_lambda_parameterization(Lambda_yield, inEvent->GetB(), enhancedLambda, rand);
            Double_t fPolY = get_mean_polarization(2.87, get_centrality(inEvent->GetB()));

            TLorentzVector lambda_lab = enhancedLambda.GetMomentum();
            hLambdaPt->Fill(lambda_lab.Pt());

            ROOT::Math::XYZVector pol = get_pol_lambda(enhancedLambda, fPolY / 100., fSigmaPolVal);
            vecPolarization->push_back(pol);

            hLambdaPhi->Fill(get_positive_phi(enhancedLambda.GetMomentum().Phi()));

            outEvent->AddParticle(enhancedLambda);

            fEnhanceStat--;
        }
        outTree->Fill();
    }

    outFile->cd();
    outTree->Write("", TObject::kOverwrite);

    TCanvas *c1 = new TCanvas("c1", "Lambda_Decays", 1200, 800);
    c1->Divide(2, 2);
    c1->cd(1);
    hLambdaPt->Draw();
    c1->cd(2);
    hCosTheta->Draw();
    c1->cd(3);
    hAngVsPt->Draw("colz");
    c1->cd(4);
    hLambdaPhi->Draw();

    c1->SaveAs("lambda_decays_plots.png");
    c1->SaveAs("lambda_decays_plots.C");

    delete c1;
    delete outEvent;
    delete inChain;

    outFile->Close();
    Lambda_yield->Close();
}