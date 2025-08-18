

#ifdef __CLING__
//#pragma link C++ class LambdaAnalysisData+;
//#pragma link C++ class UUEvent+;
//#pragma link C++ class std::vector<UParticle>+;
//#pragma link C++ class std::vector<TVector3>+;
#endif

#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TProfile.h>
#include <TLegend.h>
//#include <TRotationZ.h>
#include <TROOT.h>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TChain.h>
#include "UEvent.h"
#include "UParticle.h"

#include "TInterpreter.h"


void simulate_lambda_decays(TString inputFile, TString outputFile, TString confInFile, Int_t enhanceStat = 1);
void calc_global_polarization(TString InFileName, TString OutFileName, Int_t enhancedFlag = -100);
void set_lambda_parameterization(TFile* Lambda_yield, Double_t fBVal, UParticle &ULambda); 

Double_t get_mean_polarization(Double_t sNN, Double_t centrality);// Value in % 
Double_t get_random_value(Double_t fMean, Double_t fSigma);
Double_t get_positive_phi(const Double_t& phi);
Double_t get_centrality  (Double_t fBVal);
Double_t get_costh(Double_t alpha, Double_t pol = 0.6);
Double_t get_V2   (Double_t sNN, Double_t centrality, Double_t lambda_pT, Double_t lambda_y);
Int_t    get_number_of_bin(Double_t fValue, Double_t fMinValue, Double_t fMaxValue, Int_t NBins); 
ROOT::Math::XYZVector get_pol_lambda(UParticle& lambda, Double_t _fpoly, Double_t _fSigmaPol = 0.3);


class TString;
class TClonesArray;
class UParticle;

void calc_pol_vs_Nenh(TString InFileName, TString OutFileName, std::vector<Int_t> &enhancedFlag);


void simulate_lambda_decays(TString inputFile, TString outputFile, TString confInFile, Int_t enhanceStat) {

    // gInterpreter->GenerateDictionary("vector<UParticle>", "vector;UParticle.h");
    // gInterpreter->GenerateDictionary("vector<TVector3>", "vector;TVector3.h");

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

    // std::vector<ROOT::Math::XYZVector>* vecPol_XYZ = nullptr;


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

    // TTree *outTree = new TTree("decays", "Lambda decay products");

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
        outTree->Branch("cos_theta_p", &cos_theta_p);
        outTree->Branch("eventID", &eventID, "eventID/I");
        outTree->Branch("ULambda", &ULambda);
        outTree->Branch("UProton", &UProton);
        outTree->Branch("UPion", &UPion);
        outTree->Branch("Polarization", &vecPol);
        outTree->Branch("LambdaPolarization", &vecPolarization);


    } else {

        // vecPol.clear();
        // ULambda.clear();
        // UProton.clear();
        // UPion.clear();
        
        std::cout << "Appending to existing tree with " << outTree->GetEntries() << " entries" << std::endl;

        // Connect branches for existing tree
        outTree->SetBranchAddress("event", &outEvent);
        outTree->SetBranchAddress("cos_theta_p", &cos_theta_p);
        outTree->SetBranchAddress("eventID", &eventID);
        outTree->SetBranchAddress("ULambda", &ULambda);
        outTree->SetBranchAddress("UProton", &UProton);
        outTree->SetBranchAddress("UPion", &UPion);
        outTree->SetBranchAddress("Polarization", &vecPol);
        outTree->SetBranchAddress("LambdaPolarization", &vecPolarization);

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

        //std::cout << "=====\nLambda counter in event #" << iEvent << " : " << lambdaCounter << " lambdas\n";

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
            // outEvent->AddParticle(proton);
            // outEvent->AddParticle(pion);

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
                // unit.RotateUz(pol.Unit());
    
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
                // outEvent->AddParticle(protonEnhanced);
                // outEvent->AddParticle(pionEnhanced);

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


void calc_global_polarization(TString InFileName, TString OutFileName, Int_t enhancedFlag){
    gStyle->SetOptStat(0);  // Disable stats globally
    TF1* fitPhiDistro = new TF1("fitPhiDistro", "[0]*(1+2*[1]*TMath::Sin(x)+2*[2]*TMath::Cos(x))", 0, 2*TMath::Pi()); //Fitting function
    fitPhiDistro->SetParameter(1, 4.*(TMath::Pi()*0.732)/(8.*100)); //initial guess for fit
    // TF1* fitPhiDistro = new TF1("fitPhiDistro", 
    //     "[0]*(1 + [1]*sin(x) + [2]*cos(x) + [3]*sin(2*x) + [4]*cos(2*x))", 
    //     0, 2*TMath::Pi());
    const double anisotropy = 0.732;  // Strength of anisotropy (-1 to 1)

    TF1* fitPhiCollectiveFlowDistro = new TF1("fitPhiCollectiveFlowDistro", "[0]*(1+2*[1]*TMath::Cos(2*x))", 0, 2*TMath::Pi()); //Fitting function
    // const double anisotropy = 0.732;  // Strength of anisotropy (-1 to 1)

    //--Parameters for binning
    Int_t NpTBins = 10;//number of pT bins
    Int_t NYBins = 10; //number of rapidity bins
    Double_t fpTMax_bin = 1.6; // max pT [GeV/c]
    Double_t fpTMin_bin = 0.2; // min pT [GeV/c]
    Double_t fYMax_bin =  1.5; // max rapidity
    Double_t fYMin_bin = -1.5; // min rapidity

    //--Histograms--//
    TH1D* hProtonLambdaFrame_phi = new TH1D("hProtonLambdaFrame_phi", " ; #Delta#phi, rad;Counts", 50,  0, 2*TMath::Pi()); //proton phi* distribution (lambda frame)
    TH1D* hLambdaLabFrame_pT = new TH1D("hLambdaLabFrame_pT", " lambda lab frame ; p_{T}, GeV/c; Counts", 50,  0., 2.0); //proton pT distribution (lab frame)
    TH1D* hlambdaLabFrame_Y  = new TH1D("hlambdaLabFrame_Y", " lambda lab frame ; Y; Counts", 50,  -1.5, 1.5); //lambda Y distribution (lab frame)
    TH1D* hlambdaLabFrame_phiDistr  = new TH1D("hlambdaLabFrame_phiDistr", " lambda lab frame ; #Delta#phi; Counts", 100, 0, 2*TMath::Pi()); //lambda Y distribution (lab frame)
    TProfile* fProfV2pT = new TProfile("fProfV2pT", "v_{2} vs pT", 20, 0, 2, -10, 10);
    TProfile* fProfV1Y = new TProfile("fProfV1Y", "v_{1} vs Y", 20, -1, 1, -2, 2);

    TH1D* hPsiEP = new TH1D("hPsiEP","hPsiEP ; #Psi_{EP}; Counts",100, -TMath::Pi(), TMath::Pi());
    TH2D* hQvector = new TH2D("hQvector", "Q-vector Components", 100, -100, 100, 100, -100, 100);
    TH1D* hResolution = new TH1D("hResolution","hResolution", 100,-1, 1);

    TH1D* hProtonLabFrame_phiDistr  = new TH1D("hProtonLabFrame_phiDistr", " proton lab frame ; #Delta#phi; Counts", 100, 0, 2*TMath::Pi()); 
    TH1D* hPionLabFrame_phiDistr  = new TH1D("hPionLabFrame_phiDistr", " pion lab frame ; #Delta#phi; Counts", 100, 0, 2*TMath::Pi()); 
    TH1D* hProtonLambdaFrame_phiDistr  = new TH1D("hProtonLambdaFrame_phiDistr", " proton lambda frame  ; #Delta#phi; Counts", 100, 0, 2*TMath::Pi()); 
    TH1D* hPionLambdaFrame_phiDistr  = new TH1D("hPionLambdaFrame_phiDistr", " pion lambda frame  ; #Delta#phi; Counts", 100, 0, 2*TMath::Pi()); 

//---------_TEST
//--LAB
    TH2D* hProtonMomentumLab_XY = new TH2D("hProtonMomentumLab_XY","hProtonMomentumLab_XY ; p_{X}; p_{Y}", 100, -3, 3, 100, -3, 3);
    TH2D* hProtonMomentumLab_YZ = new TH2D("hProtonMomentumLab_YZ","hProtonMomentumLab_YZ ; p_{Y}; p_{Z}", 100, -3, 3, 100, -3, 3);
    TH1D* hProtonMomentumLab_Mag  = new TH1D("hProtonMomentumLab_Mag", "hProtonMomentumLab_Mag; Momentum.Rho(); Counts", 100, -3, 3); 

    TH2D* hPionMomentumLab_XY = new TH2D("hPionMomentumLab_XY","hPionMomentumLab_XY ; p_{X}; p_{Y}",100, -3, 3, 100, -3, 3);
    TH2D* hPionMomentumLab_YZ = new TH2D("hPionMomentumLab_YZ","hPionMomentumLab_YZ ; p_{Y}; p_{Z}",100, -3, 3, 100, -3, 3);
    TH1D* hPionMomentumLab_Mag  = new TH1D("hPionMomentumLab_Mag", "hPionMomentumLab_Mag; Momentum.Rho(); Counts", 100, -3, 3); 

//--REST
    TH2D* hProtonMomentumRest_XY = new TH2D("hProtonMomentumRest_XY","hProtonMomentumRest_XY ; p_{X}; p_{Y}",100, -3, 3, 100, -3, 3);
    TH2D* hProtonMomentumRest_YZ = new TH2D("hProtonMomentumRest_YZ","hProtonMomentumRest_YZ ; p_{Y}; p_{Z}",100, -3, 3, 100, -3, 3);
    TH1D* hProtonMomentumLRest_Mag  = new TH1D("hProtonMomentumLRest_Mag", "hProtonMomentumLRest_Mag; Momentum.Rho(); Counts", 100, -3, 3); 

    TH2D* hPionMomentumRest_XY = new TH2D("hPionMomentumRest_XY","hPionMomentumRest_XY ; p_{X}; p_{Y}",100, -3, 3, 100, -3, 3);
    TH2D* hPionMomentumRest_YZ = new TH2D("hPionMomentumRest_YZ","hPionMomentumRest_YZ ; p_{Y}; p_{Z}",100, -3, 3, 100, -3, 3);
    TH1D* hPionMomentumRest_Mag  = new TH1D("hPionMomentumRest_Mag", "hPionMomentumRest_Mag; Momentum.Rho(); Counts", 100, -3, 3);  


    TH2D* hProtonPionMomentumLab_XX = new TH2D("hProtonPionMomentumLab_XX","hProtonPionMomentumLab_XX ; proton p_{X}; pion p_{Y}",100, -3, 3, 100, -3, 3);
    TH2D* hProtonPionMomentumLab_YY = new TH2D("hProtonPionMomentumLab_YY","hProtonPionMomentumLab_YY ; proton p_{Y}; pion p_{Y}",100, -3, 3, 100, -3, 3);
    TH2D* hProtonPionMomentumLab_ZZ = new TH2D("hProtonPionMomentumLab_ZZ","hProtonPionMomentumLab_ZZ ; proton p_{Z}; pion p_{Z}",100, -3, 3, 100, -3, 3);
    TH2D* hProtonPionMomentumLab_Mag = new TH2D("hProtonPionMomentumLab_Mag","hProtonPionMomentumLab_Mag ; proton |p|; pion |p|",100, -3, 3, 100, -3, 3);


//---------_TEST

    TH2D* hYPt = new TH2D("hYPt", "hYPt", 100, -2, 2, 100, 0, 3);
    // Set axis titles
    fProfV2pT->GetYaxis()->SetTitle("v_{2}");
    fProfV2pT->GetXaxis()->SetTitle("p_{T} (GeV/c)");

    TH1D* hProtonLambdaFrame_phi_pTBins[NpTBins]; //pT binnig
    TH1D* hProtonLambdaFrame_phi_YBins [NYBins ]; //rapidity binnig

    TH1D* hlambdaLabFrame_phiDistrBin[NpTBins]; 
    //loop for TH1D object instance 
    for(Int_t iHisto = 0; iHisto < NpTBins; iHisto++) hProtonLambdaFrame_phi_pTBins[iHisto] = new TH1D(TString::Format("hProtonLambdaFrame_phi_pT_%i",iHisto), " ; #Delta#phi, rad;Counts", 50,  0, 2*TMath::Pi());
    for(Int_t iHisto = 0; iHisto < NYBins;  iHisto++) hProtonLambdaFrame_phi_YBins [iHisto] = new TH1D(TString::Format("hProtonLambdaFrame_phi_Y_%i",iHisto), " ; #Delta#phi, rad;Counts", 50,  0, 2*TMath::Pi());
    
    Int_t NpTBins_v2 = 4;
    for(Int_t iHisto = 0; iHisto < NpTBins_v2; iHisto++) hlambdaLabFrame_phiDistrBin [iHisto] = new TH1D(TString::Format("hlambdaLabFrame_phiDistrBin_%i",iHisto), " ; #Delta#phi, rad;Counts", 50,  0, 2*TMath::Pi());

    
    // Int_t dummy[2];
    // UParticle *lambda = new UParticle(1, 1, 1, 1, 1, 1, 1, dummy, 1., 1., 1., 1., 1., 1., 1., 1., 1.); //lambda particle instance
    // UParticle *proton = new UParticle(1, 1, 1, 1, 1, 1, 1, dummy, 1., 1., 1., 1., 1., 1., 1., 1., 1.); //proton particle instance
    // UParticle *pion   = new UParticle(1, 1, 1, 1, 1, 1, 1, dummy, 1., 1., 1., 1., 1., 1., 1., 1., 1.); //pion particle instance

    // TFile* inFile  = TFile::Open(InFileName, "READ"); //Open input file with polarized proton
    // TTree *inTree = (TTree*)inFile->Get("decays");
    
    TChain *inChain = new TChain("decays");
    //inChain->Add(InFileName);
    TString pathIn = "/scratch3/dflusova/afterburner/out/"; //directory with input files
    TString fileIn = "result_urqmd_xexe_2.87gev_mf_6195240"; //file for storing produced protons 
    for (int i = 1; i < 2000; i++){
        // if (!(TFile::Open(TString::Format("%s%s_%i.mcini.root",pathIn.Data(), fileIn.Data(), i), "READ")->IsZombie())){
        TString fname = TString::Format("%s%s_%i.mcini.root",pathIn.Data(), fileIn.Data(), i);
        TFile* f = TFile::Open(fname, "READ");
        // TTree* t = (TTree*) f->Get("decays");
        if (!f || f->IsZombie()){
            std::cout << "ne\n";
        }
        else{
            inChain->Add(fname);
            f->Close();
        }
        
    }

    TFile *outFile = TFile::Open(OutFileName, "RECREATE");//Resulting out file
    
    // Setup branches

    std::vector<UParticle> *ULambda = nullptr;
    std::vector<UParticle> *UProton = nullptr; 
    std::vector<UParticle> *UPion = nullptr; 
    std::vector<ROOT::Math::XYZVector> *vecPol = nullptr;

    UEvent *inEvent = new UEvent ();

    inChain->SetBranchAddress("event", &inEvent);

    inChain->SetBranchAddress("Polarization", &vecPol);
    inChain->SetBranchAddress("UProton", &UProton);
    inChain->SetBranchAddress("ULambda", &ULambda);
    inChain->SetBranchAddress("UPion", &UPion);

    // Double_t fBMin = 3.44;
    // Double_t fBMax = 6.64;
    Double_t fCenMin = 10;//%
    Double_t fCenMax = 40;//%
    Double_t fpTMin = 0.4;//GeV/C
    Double_t fpTMax = 2.0;//GeV/C
    Double_t fYmin = -0.75;
    Double_t fYmax = 0.75;
    Double_t n = 2.0; //harmonic
    //Loop over events
    // Long64_t nEvents = inTree->GetEntries();
    Long64_t nEvents = inChain->GetEntries();

    std::cout << "entries size: " << nEvents << std::endl;
    // TH1D* hTest = new TH1D("hTest","hTest", 50, 0, 2*TMath::Pi());
    for(Long64_t iEvent = 0; iEvent < nEvents; iEvent++){
        // inTree->GetEntry(iEvent);
        inChain->GetEntry(iEvent);

        // std::cout<<"Event No"<<iEvent<<std::endl;
        // std::cout << "vector size: " << vecPol->size() << std::endl;
        // std::cout<<"vecPol->at(iTrack).X()"<<std::endl;
        // std::cout<<vecPol->at(0).X()<<std::endl;

        Int_t lambdaCounter = 0;
        Double_t Qx = 0, Qy = 0;
        Double_t QxA = 0, QyA = 0, QxB = 0, QyB = 0;
        Double_t fPsi = 0;
        for(size_t iTrack = 0; iTrack < vecPol->size(); iTrack++){
            
            if(get_centrality(inEvent->GetB()) < fCenMin || get_centrality(inEvent->GetB()) > fCenMax) continue;
            // if(!ULambda->at(iTrack)) continue;
            if (ULambda->size() < 1) continue;
            // if (ULambda->at(iTrack).GetMomentum().Rapidity() < fYmin  || ULambda->at(iTrack).GetMomentum().Rapidity() > fYmax ) continue;
            // if (ULambda->at(iTrack).GetMomentum().Pt()       < fpTMin || ULambda->at(iTrack).GetMomentum().Pt()       > fpTMax) continue;
            
            // Qx += TMath::Cos(n * lambda.GetMomentum().Phi());
            // Qy += TMath::Sin(n * lambda.GetMomentum().Phi());
            Qx += TMath::Cos(n * ULambda->at(iTrack).GetMomentum().Phi());
            Qy += TMath::Sin(n * ULambda->at(iTrack).GetMomentum().Phi());
            hQvector->Fill(Qx,Qy);
            if(ULambda->at(iTrack).GetMomentum().Rapidity() > 0 ){
                QxA+=TMath::Cos(n * ULambda->at(iTrack).GetMomentum().Phi());
                QyA+=TMath::Sin(n * ULambda->at(iTrack).GetMomentum().Phi());
            }
            else{
                QxB+=TMath::Cos(n * ULambda->at(iTrack).GetMomentum().Phi());
                QyB+=TMath::Sin(n * ULambda->at(iTrack).GetMomentum().Phi());
            }

            Double_t phiStar = get_positive_phi(UProton->at(iTrack).GetMomentum().Phi()); //Get phi*
            hProtonLambdaFrame_phi->Fill(phiStar); //phi* distribution 

            //if bin isn't uderflow or overflow then put phi*-distribution in histogram[number of bin]
            if(get_number_of_bin(ULambda->at(iTrack).GetMomentum().Pt(),       fpTMin_bin, fpTMax_bin, NpTBins) > -1) hProtonLambdaFrame_phi_pTBins[(get_number_of_bin(ULambda->at(iTrack).GetMomentum().Pt(),       fpTMin_bin, fpTMax_bin, NpTBins))]->Fill(phiStar); 
            if(get_number_of_bin(ULambda->at(iTrack).GetMomentum().Rapidity(), fYMin_bin,  fYMax_bin,  NYBins ) > -1) hProtonLambdaFrame_phi_YBins [(get_number_of_bin(ULambda->at(iTrack).GetMomentum().Rapidity(), fYMin_bin,  fYMax_bin,  NYBins ))]->Fill(phiStar); 


            //Filling pT and Y histograms
            hLambdaLabFrame_pT->Fill(ULambda->at(iTrack).GetMomentum().Pt());
            hlambdaLabFrame_Y ->Fill(ULambda->at(iTrack).GetMomentum().Rapidity());
            hlambdaLabFrame_phiDistr->Fill(get_positive_phi(ULambda->at(iTrack).GetMomentum().Phi()));
            //std::cout << "ymin:" << fYmin << " ymax: " << fYmax << " Y: " << lambda.GetMomentum().Rapidity() << std::endl;
            // if (lambda.GetMomentum().Rapidity() > fYmin  && lambda.GetMomentum().Rapidity() < fYmax ){
            //     std::cout << "Pt " << ULambda->at(iTrack).GetMomentum().Pt() << " CosPhi " << TMath::Cos(2*(ULambda->at(iTrack).GetMomentum().Phi())) << std::endl;
            //      fProfV2pT->Fill(ULambda->at(iTrack).GetMomentum().Pt(),      TMath::Cos(2*(ULambda->at(iTrack).GetMomentum().Phi())));
            // }
            // if (lambda.GetMomentum().Pt()       > fpTMin && lambda.GetMomentum().Pt()       < fpTMax) fProfV1Y ->Fill(ULambda->at(iTrack).GetMomentum().Rapidity(),TMath::Cos(get_positive_phi(fPsi-ULambda->at(iTrack).GetMomentum().Phi())));
            if (ULambda->at(iTrack).GetMomentum().Rapidity() > fYmin  && ULambda->at(iTrack).GetMomentum().Rapidity() < fYmax ) fProfV2pT->Fill(ULambda->at(iTrack).GetMomentum().Pt(),      TMath::Cos(2*(ULambda->at(iTrack).GetMomentum().Phi())));
            if (ULambda->at(iTrack).GetMomentum().Pt()       > fpTMin && ULambda->at(iTrack).GetMomentum().Pt()        < fpTMax) fProfV1Y ->Fill(ULambda->at(iTrack).GetMomentum().Rapidity(),TMath::Cos(get_positive_phi(ULambda->at(iTrack).GetMomentum().Phi())));
           
            hYPt->Fill(ULambda->at(iTrack).GetMomentum().Rapidity(), ULambda->at(iTrack).GetMomentum().Pt());
            lambdaCounter++;


//---TEST
TLorentzVector vLambda_lab  = ULambda->at(iTrack).GetMomentum();
TLorentzVector vProton_rest = UProton->at(iTrack).GetMomentum();
TLorentzVector vPion_rest   = UPion  ->at(iTrack).GetMomentum();

TLorentzVector vProton_lab = vProton_rest;
TLorentzVector vPion_lab   = vPion_rest;

TVector3 beta = -vLambda_lab.BoostVector();

vProton_lab.Boost(-beta);
vPion_lab.Boost(-beta);

hProtonLabFrame_phiDistr->Fill(get_positive_phi(vProton_lab.Phi()));
hPionLabFrame_phiDistr->Fill(get_positive_phi(vPion_lab.Phi()));
hProtonLambdaFrame_phiDistr->Fill(get_positive_phi(vProton_rest.Phi()));
hPionLambdaFrame_phiDistr  ->Fill(get_positive_phi(vPion_rest.Phi()));



//--LAB
hProtonMomentumLab_XY->Fill(vProton_lab.Px(),vProton_lab.Py());
hProtonMomentumLab_YZ->Fill(vProton_lab.Py(),vProton_lab.Pz());
hProtonMomentumLab_Mag->Fill(vProton_lab.Rho());

hPionMomentumLab_XY->Fill(vPion_lab.Px(),vPion_lab.Py());
hPionMomentumLab_YZ->Fill(vPion_lab.Py(),vPion_lab.Pz());
hPionMomentumLab_Mag->Fill(vPion_lab.Rho());


//--REST
hProtonMomentumRest_XY->Fill(vProton_rest.Px(),vProton_rest.Py());
hProtonMomentumRest_YZ->Fill(vProton_rest.Py(),vProton_rest.Pz());
hProtonMomentumLRest_Mag->Fill(vProton_rest.Rho());

hPionMomentumRest_XY->Fill(vPion_rest.Px(),vPion_rest.Py());
hPionMomentumRest_YZ->Fill(vPion_rest.Py(),vPion_rest.Pz());
hPionMomentumRest_Mag->Fill(vPion_rest.Rho());



hProtonPionMomentumLab_XX->Fill(vProton_lab.Px(),vPion_lab.Px());
hProtonPionMomentumLab_YY->Fill(vProton_lab.Py(),vPion_lab.Py());
hProtonPionMomentumLab_ZZ->Fill(vProton_lab.Pz(),vPion_lab.Pz());
hProtonPionMomentumLab_Mag->Fill(vProton_lab.Rho(),vPion_lab.Rho());


//---!TEST        

        }
        //std::cout<<"number of lambda within event "<<lambdaCounter<<std::endl;
        Double_t Psi_n = (1./n) * TMath::ATan2(Qy, Qx);
        Double_t Psi_A = (1./n) * TMath::ATan2(QyA, QxA);
        Double_t Psi_B = (1./n) * TMath::ATan2(QyB, QxB);
        //std::cout << "res: " << Psi_A << " " <<  Psi_B << std::endl;
        Double_t resolution = TMath::Sqrt(TMath::Cos(n * (Psi_A - Psi_B)));

        // std::cout<<"Psi = "<<Psi_n<<std::endl;
        hPsiEP->Fill(Psi_n);
        hResolution->Fill(resolution);
    }


    //--pT binning--//
    TH1D* hphi_pTBins= new TH1D("hphi_pTBins", " p_{T} binning; p_{T}, GeV/c; P_{#lambda} %", NpTBins,  fpTMin_bin, fpTMax_bin); // Histogram with pT bins
    for(Int_t iHisto = 0; iHisto < NpTBins; iHisto++) { //Loop over pT bins
        hProtonLambdaFrame_phi_pTBins[iHisto]->Fit(fitPhiDistro ,"L");
        hProtonLambdaFrame_phi_pTBins[iHisto]->Write();
        hphi_pTBins->SetBinContent(iHisto+1, (fitPhiDistro->GetParameter(1)*8)/(TMath::Pi()*anisotropy)*100);
        hphi_pTBins->SetBinError  (iHisto+1, (fitPhiDistro->GetParError(1)*8)/(TMath::Pi()*anisotropy)*100);
    }

    hphi_pTBins->SetMarkerStyle(23); // Filled circle (ROOT default)
    hphi_pTBins->SetMarkerSize(1.2); // Slightly larger markers
    hphi_pTBins->Write();

    //--Y binning--//
    //same logic like in previous one
    TH1D* hphi_YBins= new TH1D("hphi_YBins", " rapidity binning; Y; P_{#lambda} %", NYBins,  fYMin_bin, fYMax_bin);
    for(Int_t iHisto = 0; iHisto < NYBins; iHisto++) {
        hProtonLambdaFrame_phi_YBins[iHisto]->Fit(fitPhiDistro, "L");
        hProtonLambdaFrame_phi_YBins[iHisto]->Write();
        hphi_YBins->SetBinContent(iHisto+1, (fitPhiDistro->GetParameter(1)*8)/(TMath::Pi()*anisotropy)*100);
        hphi_YBins->SetBinError  (iHisto+1, (fitPhiDistro->GetParError(1)*8)/ (TMath::Pi()*anisotropy)*100);
    }

    hphi_YBins->SetMarkerStyle(23); // Filled circle (ROOT default)
    hphi_YBins->SetMarkerSize(1.2); // Slightly larger markers
    hphi_YBins->Write();

    hProtonLambdaFrame_phi->Fit(fitPhiDistro);//Fitting phi* full distribution  (without pT or Y binning)
    TLegend* lProtonLambdaFrame_phi = new TLegend(0.50, 0.58, 0.9, 0.89);
    lProtonLambdaFrame_phi->SetHeader("p_{0}(1 + 2p_{1}sin(#phi) + 2p_{2}cos(#phi))", "C");
    lProtonLambdaFrame_phi->SetTextSize(0.05); 
    lProtonLambdaFrame_phi->SetBorderSize(0);   
    lProtonLambdaFrame_phi->AddEntry("", Form("p_{1} = %.3f #pm %.3f", fitPhiDistro->GetParameter(1), fitPhiDistro->GetParError(1)), "");
    lProtonLambdaFrame_phi->AddEntry("", Form("P_{#Lambda} = %.2f #pm %.2f %%", (fitPhiDistro->GetParameter(1)*8)/(TMath::Pi()*anisotropy)*100, (fitPhiDistro->GetParError(1)*8)/(TMath::Pi()*anisotropy)*100), "");
    lProtonLambdaFrame_phi->AddEntry("", Form("P_{#Lambda} = (8p_{1}) / (%.3f#pi)", anisotropy), "");

    // Plot results
    TCanvas *cProton = new TCanvas("cProton", "proton", 1200, 800);
    cProton->Divide(2,2);
    cProton->cd(1); hProtonLambdaFrame_phi->Draw(); lProtonLambdaFrame_phi->Draw("same");
    cProton->cd(2); hLambdaLabFrame_pT->Draw();
    cProton->cd(3); hphi_pTBins->Draw("PE");
    cProton->cd(4); hphi_YBins->Draw("PE");

    // cProton->SaveAs("proton_polarization_plots.png");
    cProton->SaveAs((enhancedFlag > -1 ) ? TString::Format("picture/proton_polarization_plots_%i.png",enhancedFlag) : "picture/proton_polarization_plots.png");
    cProton->Write();
    hLambdaLabFrame_pT->Write();
    hlambdaLabFrame_Y->Write();
    hProtonLambdaFrame_phi->Write();
    hlambdaLabFrame_phiDistr->Fit(fitPhiCollectiveFlowDistro);
    hlambdaLabFrame_phiDistr->Write();
    TCanvas *cCollectiveFlow = new TCanvas("cCollectiveFlow", "cCollectiveFlow", 1200, 800);
    hlambdaLabFrame_phiDistr->Draw();
    // Get fit parameters
    double v2 = fitPhiCollectiveFlowDistro->GetParameter(1);
    double v2_error = fitPhiCollectiveFlowDistro->GetParError(1);
    double constant = fitPhiCollectiveFlowDistro->GetParameter(0);
    double constant_error = fitPhiCollectiveFlowDistro->GetParError(0);

    TLegend *llambdaLabFrame_phiDistr = new TLegend(0.7, 0.7, 0.89, 0.89);
    llambdaLabFrame_phiDistr->AddEntry(fitPhiCollectiveFlowDistro, Form("v_{2} = %.3f #pm %.3f", v2, v2_error), "l");
    llambdaLabFrame_phiDistr->SetBorderSize(0);   
    llambdaLabFrame_phiDistr->Draw("same");

    cCollectiveFlow->SaveAs((enhancedFlag > -1 ) ? TString::Format("picture/lambda_collective_flow_%i.png",enhancedFlag) : "picture/lambda_collective_flow.png");
    fProfV2pT->Write();
    fProfV1Y->Write();
    hPsiEP->Write();
    hQvector->Write();
    hResolution->Write();
    hYPt->Write();
    // outFile->Close();
    // inFile->Close();
    // std::cout<<"Close() "<<std::endl;

    // std::cout<<"delete 1"<<std::endl;



    //--LAB
hProtonMomentumLab_XY->Write();
hProtonMomentumLab_YZ->Write();
hProtonMomentumLab_Mag->Write();

hPionMomentumLab_XY->Write();
hPionMomentumLab_YZ->Write();
hPionMomentumLab_Mag->Write();


//--REST
hProtonMomentumRest_XY->Write();
hProtonMomentumRest_YZ->Write();
hProtonMomentumLRest_Mag->Write();

hPionMomentumRest_XY->Write();
hPionMomentumRest_YZ->Write();
hPionMomentumRest_Mag->Write();


hProtonPionMomentumLab_XX->Write();
hProtonPionMomentumLab_YY->Write();
hProtonPionMomentumLab_ZZ->Write();
hProtonPionMomentumLab_Mag->Write();
//---!TEST        

    hProtonLabFrame_phiDistr->Write();
    hPionLabFrame_phiDistr->Write();
    hProtonLambdaFrame_phiDistr->Write();
    hPionLambdaFrame_phiDistr->Write();

    delete hLambdaLabFrame_pT; delete hlambdaLabFrame_Y; delete hlambdaLabFrame_phiDistr; delete hProtonLambdaFrame_phi; delete fProfV1Y; delete hPsiEP; delete hQvector; delete fProfV2pT;//delete hProtonLambdaFrame_phi_YBins; delete hProtonLambdaFrame_phi_pTBins;
    // delete hResolution;
    // std::cout<<"delete 2"<<std::endl;

    for(TH1D* histo : hProtonLambdaFrame_phi_pTBins) delete histo;
    for(TH1D* histo : hProtonLambdaFrame_phi_YBins ) delete histo;
    // for(TH1D* histo : hlambdaLabFrame_phiDistrBin  ) delete histo;
    // std::cout<<"delete 3"<<std::endl;

    delete cCollectiveFlow; delete cProton;

    // std::cout<<"delete 4"<<std::endl;


}

void calc_pol_vs_Nenh(TString InFileName, TString OutFileName, std::vector<Int_t> &enhancedFlag){

    // TFile *inFile = TFile::Open(InFileName);

    TChain *inChain = new TChain("decays");
    //inChain->Add(InFileName);
    TString pathIn = "/scratch3/dflusova/afterburner/out/"; //directory with input files
    TString fileIn = "result_urqmd_xexe_2.87gev_mf_6195240"; //file for storing produced protons 
    for (int i = 1; i < 2000; i++){
        TString fname = TString::Format("%s%s_%i.mcini.root",pathIn.Data(), fileIn.Data(), i);
        TFile* f = TFile::Open(fname, "READ");
        // TTree* t = (TTree*) f->Get("decays");
        if (!f || f->IsZombie()){
            std::cout << "ne\n";
        }
        else{
            inChain->Add(fname);
            f->Close();
        }    
    }

    std::vector<UParticle> *ULambda = nullptr;
    std::vector<UParticle> *UProton = nullptr; 
    std::vector<UParticle> *UPion = nullptr; 
    std::vector<ROOT::Math::XYZVector > *vecPol = nullptr;

    UEvent *inEvent = new UEvent ();

    inChain->SetBranchAddress("event", &inEvent);

    inChain->SetBranchAddress("Polarization", &vecPol);
    inChain->SetBranchAddress("UProton", &UProton);
    inChain->SetBranchAddress("ULambda", &ULambda);
    inChain->SetBranchAddress("UPion", &UPion);

    Long64_t nEvents = inChain->GetEntries();

    std::cout << "entries size: " << nEvents << std::endl;
    Int_t nBinsEnh = *(enhancedFlag.end() - 1) - enhancedFlag.at(0);
    TH2D* hPx_Nenh   = new TH2D("hPx_Nenh","hPx_Nenh", nBinsEnh, enhancedFlag.at(0), *(enhancedFlag.end() - 1), 50, -1., 1.);
    TH2D* hPy_Nenh   = new TH2D("hPy_Nenh","hPy_Nenh", nBinsEnh, enhancedFlag.at(0), *(enhancedFlag.end() - 1), 50, -1., 1.);
    TH2D* hPmag_Nenh = new TH2D("hPmag_Nenh","hPmag_Nenh", nBinsEnh, enhancedFlag.at(0), *(enhancedFlag.end() - 1), 50, -1., 1.);

    TH2D* hPx_Nlamb   = new TH2D("hPx_Nlamb","hPx_Nlamb", nBinsEnh, enhancedFlag.at(0), *(enhancedFlag.end() - 1), 50, -1., 1.);
    TH2D* hPy_Nlamb   = new TH2D("hPy_Nlambd","hPy_Nlambd", nBinsEnh, enhancedFlag.at(0), *(enhancedFlag.end() - 1), 50, -1., 1.);
    TH2D* hPmag_lamb = new TH2D("hPmag_lambd","hPmag_lambd", nBinsEnh, enhancedFlag.at(0), *(enhancedFlag.end() - 1), 50, -1., 1.);

    TProfile* hPx_pT = new TProfile("hPx_pT", "hPx_pT", 50, 0., 3.0, -50, 50);
    TProfile* hPy_pT = new TProfile("hPy_pT", "hPy_pT", 50, 0., 3.0, -50, 50);
    TProfile* hPz_pT = new TProfile("hPz_pT", "hPz_pT", 50, 0., 3.0, -50, 50);


    TProfile* hPx_Y = new TProfile("hPx_y", "hPx_y", 50, -1., 1., -50, 50);
    TProfile* hPy_Y = new TProfile("hPy_y", "hPy_y", 50, -1., 1., -50, 50);
    TProfile* hPz_Y = new TProfile("hPz_Y", "hPz_Y", 50, -1., 1., -50, 50);

    Double_t fpTMin = 0.4;//GeV/C
    Double_t fpTMax = 2.0;//GeV/C
    Double_t fYmin = -0.75;
    Double_t fYmax = 0.75;

    for(Long64_t iEvent = 0; iEvent < nEvents; iEvent++){
        // inTree->GetEntry(iEvent);
        inChain->GetEntry(iEvent);
        int lambdas = 0;
        for (size_t pols_i = 0; pols_i < vecPol->size(); pols_i++){
            //std::cout << "pz: " << vecPol->at(pols_i).Z() * 100 << std::endl;
            Double_t pT = ULambda->at(pols_i).GetMomentum().Pt();
            Double_t Y = ULambda->at(pols_i).GetMomentum().Rapidity();
            if (pT < fpTMin || pT > fpTMax || Y < fYmin || Y > fYmax) continue;
            hPx_pT->Fill(ULambda->at(pols_i).GetMomentum().Pt(), vecPol->at(pols_i).X() * 100);
            hPy_pT->Fill(ULambda->at(pols_i).GetMomentum().Pt(), vecPol->at(pols_i).Y() * 100);
            hPz_pT->Fill(ULambda->at(pols_i).GetMomentum().Pt(), vecPol->at(pols_i).Z() * 100);
            hPx_Y->Fill(ULambda->at(pols_i).GetMomentum().Rapidity(), vecPol->at(pols_i).X() * 100);
            hPy_Y->Fill(ULambda->at(pols_i).GetMomentum().Rapidity(), vecPol->at(pols_i).Y() * 100);
            hPz_Y->Fill(ULambda->at(pols_i).GetMomentum().Rapidity(), vecPol->at(pols_i).Z() * 100);
        }
    }


    // for(size_t iEnh = 0; iEnh < enhancedFlag.size(); iEnh++){
    //     Int_t enhFlag = enhancedFlag[iEnh];
    //     std::cout << "eng flag: " << enhFlag << std::endl;

    //     for(Long64_t iEvent = 0; iEvent < nEvents; iEvent++){
    //         // inTree->GetEntry(iEvent);
    //         inChain->GetEntry(iEvent);
    //         int lambdas = 0;
    //         for (size_t pols_i = 0; pols_i < vecPol->size(); pols_i++){
    //             hPx_Nlamb->Fill(enhFlag, vecPol->at(pols_i).X());
    //             hPy_Nlamb->Fill(enhFlag, vecPol->at(pols_i).Y());
    //             hPmag_lamb->Fill(enhFlag, TMath::Sqrt(vecPol->at(pols_i).X()*vecPol->at(pols_i).X()+vecPol->at(pols_i).Y()*vecPol->at(pols_i).Y()));
    //             if (ULambda->at(pols_i).GetMate() == -9){
    //                 lambdas++;
    //                 hPx_Nenh->Fill(enhFlag, vecPol->at(pols_i).X());
    //                 hPy_Nenh->Fill(enhFlag, vecPol->at(pols_i).Y());
    //                 hPmag_Nenh->Fill(enhFlag, TMath::Sqrt(vecPol->at(pols_i).X()*vecPol->at(pols_i).X()+vecPol->at(pols_i).Y()*vecPol->at(pols_i).Y()));
    //                 if (lambdas > enhFlag) break;
    //             }
    //         }
    //     }
    // }


    TFile *outFile = TFile::Open(OutFileName, "RECREATE");//Resulting out file
    hPx_Nenh->Write();
    hPy_Nenh->Write();
    hPmag_Nenh->Write();

    hPx_Nlamb->Write();
    hPy_Nlamb->Write();
    hPmag_lamb->Write();

    hPx_pT->Write();
    hPy_pT->Write();
    hPz_pT->Write();
    hPx_Y->Write();
    hPy_Y->Write();
    hPz_Y->Write();

    outFile->Close();
}


Double_t get_positive_phi(const Double_t& phi){
    Double_t phi_rot = phi;
    if (phi < 0) phi_rot += 2.*TMath::Pi();
    return phi_rot;
}


Double_t get_costh(Double_t alpha, Double_t pol) {

    // TF1* randomCosTheta = new TF1("randomCosTheta", "1+(TMath::Pi()*[0]*[1])/8*x", -1.0, 1.0);
    // TF1* randomCosTheta = new TF1("randomCosTheta", "1+[0]*[1]*x", -1.0, 1.0);
    TF1* randomCosTheta = new TF1("randomCosTheta", "(1+[0]*x)", -1.0, 1.0);
    randomCosTheta->SetNpx(10000);

    // Set the parameters
    // randomCosTheta->SetParameters(alpha);

    // randomCosTheta->SetParameters(alpha, pol);
    randomCosTheta->SetParameters(alpha);

    // Generate a random value from the distribution
    Double_t value = randomCosTheta->GetRandom();
    // std::cout<<"random costh = "<<value<<std::endl;
    // Clean up
    delete randomCosTheta;
    
    return value;

}

Double_t get_random_value(Double_t fMean, Double_t fSigma)
{
    TRandom3* rand = new TRandom3(0);  // Seed with 0 for reproducibility
    Double_t randomVal = rand->Gaus(fMean, fSigma);
    delete rand;
    return randomVal;
}


Int_t get_number_of_bin(Double_t fValue, Double_t fMinValue, Double_t fMaxValue, Int_t NBins) { //Function return number of bin
    if (NBins <= 0) return -1;
    if (fValue < fMinValue) return -1;  // Underflow
    if (fValue >= fMaxValue) return -1; // Overflow
    
    Double_t binWidth = (fMaxValue - fMinValue)/NBins;
    //std::cout << "bin_value" <<  (fValue - fMinValue)/binWidth << std::endl;
    Int_t bin = static_cast<Int_t>((fValue - fMinValue)/binWidth);
    
    // Protect against edge case where fValue == fMaxValue
    return (bin < NBins) ? bin : NBins-1;
}

    void set_lambda_parameterization(TFile* Lambda_yield, Double_t fBVal, UParticle &ULambda){ //Valeriy's function for properly lambda generation

    TRandom3* rand = new TRandom3(0);  // Seed with 0 for reproducibility

	Double_t centrality = get_centrality(fBVal); //centrality for parameterization (b->centrality for Xe+Xe below)
	Double_t sNN = 2.87; // Energy of the collision in center-of-mass system

	TH2F* h_pt_y;
	if(centrality<10) h_pt_y = (TH2F*) Lambda_yield->Get("h2PartYpT010");
	else if(centrality>10 && centrality <40) h_pt_y = (TH2F*) Lambda_yield->Get("h2PartYpT");
	else h_pt_y = (TH2F*) Lambda_yield->Get("h2PartYpT40100");


    Double_t lambda_pT; // pT of Lambda from pT-y TH2F
	Double_t lambda_y;  // rapidity of Lambda from pT-y TH2F

	h_pt_y->GetRandom2(lambda_y,lambda_pT,rand);

	Double_t v1 = (28.8635/(TMath::Power(sNN,2.89092))) 
                * ((-0.0233*centrality+0.5413* TMath::Power(centrality,1./3) ) 
                * (163.536/18.0188 - 163.536/(lambda_pT +18.0188) )* lambda_y + ( -0.0056*centrality+0.377*TMath::Power(centrality,1./3) ) 
                * (0.6653* lambda_pT - 0.6172 * TMath::Power(lambda_pT,2) +0.1154 * TMath::Power(lambda_pT,3) ) * TMath::Power(lambda_y,3) ); // v1 cent-pT-y func
    // Double_t v2 = 0.07; //uniform distribution
    Double_t v2 = get_V2(sNN, centrality, lambda_pT, lambda_y);
    // std::cout<<"================ "<<v2<<std::endl;
    // if(v2 == 0) v2 =0.07;
    if( v1  > 1 )  v1 =1;
    else if(v1<-1) v1=-1;
    // if( v2  > 1 ) v2 =1;
    // std::cout<<" v2 value =  "<<v2<<std::endl;
    // generate phi according to v1 and v2
    // v1 = 0.0;   //generation WITHOUT flows
    // v2 = 0.0;

    static TF1 f("f", "[0]*(1+2*[1]*TMath::Cos(x)+2*[2]*TMath::Cos(2*x))+[3]", 0,2*TMath::Pi());
    Double_t a1=1+2*v1+2*v2;
    Double_t a2=1-2*v1+2*v2;
    Double_t a3=1-v1*v1/(4*v2)-2*v2;
    Double_t a=0; //
    if (a1<a) a=a1;  // find analytic minimun to shift
    if (a2<a) a=a2;  // find analytic minimun to shift
    if (a3<a) a=a3;  // find analytic minimun to shift
    f.SetParameter(0,1/(2*TMath::Pi()*(1-a))); // norm
    f.SetParameter(1,v1);  // v1
    f.SetParameter(2,v2);  // v2
    f.SetParameter(3,-a/(2*TMath::Pi()*(1-a))); // shift to have probability
    f.SetNpx(10000);  // to get a better result when using TF1::GetRandom
    Double_t phi=f.GetRandom(rand);

    // Double_t fEnergyLambda = ULambda.GetMomentum().E();
    // Calculate total energy (E) properly
    TLorentzVector vec;
    // Calculate transverse mass (m_T)
    Double_t lambda_mass = 1.115683; // GeV/c² (Lambda mass)
    Double_t mT = sqrt(lambda_pT * lambda_pT + lambda_mass * lambda_mass);
    // Convert rapidity (y) to pseudorapidity (η)
    Double_t pz = mT * sinh(lambda_y); // longitudinal momentum
    Double_t lambda_eta = (pz != 0.0) ? TMath::ATanH(pz / sqrt(pz*pz + lambda_pT*lambda_pT)) : 0.0;
    Double_t fEnergyLambda = sqrt(lambda_pT*lambda_pT * cosh(lambda_eta)*cosh(lambda_eta) + lambda_mass*lambda_mass);
    vec.SetPtEtaPhiE(lambda_pT, lambda_eta, phi, fEnergyLambda);
    ULambda.SetMomentum(vec); //ULambda

    delete rand;


}


Double_t get_V2(Double_t sNN, Double_t centrality, Double_t lambda_pT, Double_t lambda_y) {

	// Double_t v2 = 1.05*(0.8132-11.4/TMath::Power(sNN,2.1)) 
    //                   *((-0.01105*centrality + 0.000162*TMath::Power(centrality,2) )/(-0.01105*25 + 0.000162*TMath::Power(25,2) )) 
    //                   *(-1) * ( (0.32172*TMath::Power(lambda_y,1.805)-0.1578)*lambda_pT + (0.05838*TMath::Power(lambda_y,3.2172) - 0.0179)
    //                   *TMath::Exp(lambda_pT*(1.4562*TMath::Power(lambda_y,3.28814) + 0.22912))
    //                   *TMath::Sin(lambda_pT*(-1.004*TMath::Power(lambda_y,4.6398) + 3.88097)+(-3.6075*TMath::Power(lambda_y,4.6582) + 3.35966))); //v2 snn-cent-pT-y func
    // if (sNN <= 0 || centrality < 0 || lambda_pT < 0 || lambda_y < 0) return 0.;  // or throw an error
    lambda_y = TMath::Abs(lambda_y);
    
    // =============================================
    // 1. Energy-dependent term (depends on sNN)
    // =============================================
    Double_t energy_term = 1.05 * (0.8132 - 11.4 / TMath::Power(sNN, 2.1));

    // =============================================
    // 2. Centrality-dependent term (normalized at centrality=25)
    // =============================================
    Double_t centrality_numerator = 
        -0.01105 * centrality + 0.000162 * TMath::Power(centrality, 2);
    Double_t centrality_denominator = 
        -0.01105 * 25 + 0.000162 * TMath::Power(25, 2);
    Double_t centrality_term = centrality_numerator / centrality_denominator;

    // =============================================
    // 3. Sign flip (-1 factor)
    // =============================================
    Double_t sign_flip = -1;

    // =============================================
    // 4. Rapidity (lambda_y) and pT (lambda_pT) dependent term
    // =============================================
    // 4a. Linear term in lambda_pT (depends on lambda_y)
    Double_t linear_pT_term = 
        (0.32172 * TMath::Power(lambda_y, 1.805) - 0.1578) * lambda_pT;

    // 4b. Nonlinear term (exponential + sinusoidal dependence)
    Double_t exp_argument = 
        lambda_pT * (1.4562 * TMath::Power(lambda_y, 3.28814) + 0.22912);
    Double_t sin_phase = 
        lambda_pT * (-1.004 * TMath::Power(lambda_y, 4.6398) + 3.88097) 
        + (-3.6075 * TMath::Power(lambda_y, 4.6582) + 3.35966);
    
    Double_t nonlinear_term = 
        (0.05838 * TMath::Power(lambda_y, 3.2172) - 0.0179) 
        * TMath::Exp(exp_argument) 
        * TMath::Sin(sin_phase);

    // Combine all terms into final v2 expression
    Double_t v2 = 
        energy_term 
        * centrality_term 
        * sign_flip 
        * (linear_pT_term + nonlinear_term);

    // std::cout << "ERROR: NaN in nonlinear_term. Args: "
    //     << "exp_argument=" << exp_argument << ", sin_phase=" << sin_phase <<" , centrality_denominator = " << centrality_denominator << std::endl;

    // if (TMath::Abs(centrality_denominator) < 1e-10) {
    //     return 0.;  // Avoid division by zero
    // }

    // if (TMath::IsNaN(nonlinear_term)) {
    //     std::cout << "ERROR: NaN in nonlinear_term. Args: "
    //                 << "exp_argument=" << exp_argument << ", sin_phase=" << sin_phase << std::endl;
    // }
    // if (!TMath::Finite(sin_phase)) {
    //     return 0.;
    // }

    return v2;

}

Double_t get_mean_polarization(Double_t sNN, Double_t centrality){ return (2.8569/ TMath::Power(sNN,0.955513) ) * (2.4702 - 0.0461*centrality + 0.0042 * TMath::Power(centrality, 2)); }// Value in % 

Double_t get_centrality(Double_t fBVal){
    if(fBVal<3.44) {return 5;}
    else if(fBVal<4.88) {return 15;}
    else if(fBVal<5.84) {return 25;}
    else if(fBVal<6.64) {return 35;}
    else if(fBVal<7.44) {return 45;}
    else if(fBVal<8.08) {return 55;}
    else if(fBVal<8.72) {return 65;}
    else if(fBVal<9.36) {return 75;}
    else if(fBVal<9.84) {return 85;}
    else {return 95;}
}


ROOT::Math::XYZVector get_pol_lambda(UParticle& lambda, Double_t _fpoly, Double_t _fSigmaPol){ 
    
    Double_t fPhi = 0;

    Double_t fpolx  = 0; 
    Double_t fpolSx = 0.07;  
    Double_t fpoly  = _fpoly;
    Double_t fpolSy = 0.07; 
    Double_t fpolz  = 0;
    Double_t fpolSz = 0.07;

    TRandom3* rand = new TRandom3(0);  // Seed with 0 for reproducibility
    Double_t polX = rand->Gaus(fpolx,fpolSx);
    Double_t polY = rand->Gaus(fpoly,fpolSy);
    Double_t polZ = rand->Gaus(fpolz,fpolSz);

    // Double_t polX = gRandom->Gaus(fpolx,fpolSx);  // generate polarization direction
    // Double_t polY = gRandom->Gaus(fpoly,fpolSy);  // generate polarization direction /100
    // Double_t polZ = gRandom->Gaus(fpolz,fpolSz);  // generate polarization direction
    ROOT::Math::XYZVector polarizationVec = ROOT::Math::XYZVector(polX, polY, polZ);

    // std::cout<<"TVector3(polX, polY, polZ) "<<polX<<" "<<polY<<" "<<polZ<<std::endl;

    ROOT::Math::RotationZ rotateRP(fPhi);  // set rotation transformation
    polarizationVec = rotateRP*polarizationVec;  // rotate reaction plane

    // Float_t polmag = TMath::Sqrt(polarizationVec.mag());
    Float_t polmag = TMath::Sqrt(polarizationVec.mag2());

    // special case, overpolarized -> 100% polarized
    if (polmag>1.) {
        polarizationVec *= 1./polmag;  // scale to 1
        polmag = 1.;
    }
    // special case, no polarization -> generate random direction
    if (polmag == 0) {
    // random unitary vector, see https://mathworld.wolfram.com/SpherePointPicking.html
        Double_t x1,x2,R2,R;
        do {
            // x1 = 1-2*gRandom->Rndm();
            // x2 = 1-2*gRandom->Rndm();
            x1 = 1-2*rand->Rndm();
            x2 = 1-2*rand->Rndm();
            R2 = x1*x1+x2*x2;
        } while(R2 >= 1);
        R=2*TMath::Sqrt(1-R2);
        polarizationVec.SetXYZ(x1*R, x2*R, 1-2*R2);  // random unitary vector - background for signal
    }

    // Float_t xxx = gRandom->Rndm();
    Float_t xxx = rand->Rndm();
    if (xxx < 1./2.*(1.-polmag)) {  // probability of spin flip
        polarizationVec *= -1.;  // spin flip according to mean polarization
    }

    // std::cout<<"polmag =  "<<polmag<<std::endl;

    //std::cout<<"Polarization generated (X,Y,Z) "<<polarizationVec.X()<<" "<<polarizationVec.Y()<<" "<<polarizationVec.Z()<<std::endl;
    delete rand;
    return polarizationVec;
}

