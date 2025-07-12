#ifdef __CLING__
#pragma link C++ class LambdaAnalysisData+;
#pragma link C++ class UUEvent+;
#endif

#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include "UEvent.h"
#include "UParticle.h"

class TString;
class TClonesArray;
class UParticle;

class UUEvent : public UEvent { //Custom inherited class. Here are added new methods and polarization implementation
    private:
        Int_t fEventNr;
        Double_t fB;
        Double_t fPhi;
        Int_t fNes;
        Int_t fStepNr;
        Double_t fStepT;
        Int_t fNpa;
        TString fComment;

        //TClonesArray *fParticles;

        TClonesArray *fLambdas = new TClonesArray("UParticle");
        TClonesArray *fProtons = new TClonesArray("UParticle");
        TClonesArray *fPions = new TClonesArray("UParticle");

        TVector3 polarizationAxis;
        TRotation rotation;
    public :
        UUEvent() {
            fLambdas = new TClonesArray("UParticle");
            fProtons = new TClonesArray("UParticle");
            fPions   = new TClonesArray("UParticle");
        }
        // Copy constructor
         UUEvent(const UUEvent& right) {
            fLambdas = new TClonesArray("UParticle",100);
            fProtons = new TClonesArray("UParticle",100);
            fPions = new TClonesArray("UParticle",100);

            // Copy lambdas
            for (Int_t i = 0; i < right.GetNLambdas(); i++) {
                UParticle* p = right.GetLambda(i);
                if (p) new ((*fLambdas)[i]) UParticle(*p);
            }

            // Copy protons
            for (Int_t i = 0; i < right.GetNProtons(); i++) {
                UParticle* p = right.GetProton(i);
                if (p) new ((*fProtons)[i]) UParticle(*p);
            }

            // Copy pions
            for (Int_t i = 0; i < right.GetNPions(); i++) {
                UParticle* p = right.GetPion(i);
                if (p) new ((*fPions)[i]) UParticle(*p);
            }
            fEventNr = right.GetEventNr();
            fB = right.GetB();
            fPhi = right.GetPhi();
            fNes = right.GetNes();
            fStepNr = right.GetStepNr();
            fStepT = right.GetStepT();
            fNpa = right.GetNpa();
            fComment = "";
        }

        void Init(const UEvent& right){
            fEventNr = right.GetEventNr();
            fB = right.GetB();
            fPhi = right.GetPhi();
            fNes = right.GetNes();
            fStepNr = right.GetStepNr();
            fStepT = right.GetStepT();
            fNpa = right.GetNpa();
            fComment = "";
        }

        virtual ~UUEvent() {
            delete fLambdas;
            delete fProtons;
            delete fPions;
            // delete fParticles;
        }
        Int_t GetNLambdas() const { return fLambdas->GetEntriesFast(); }
        Int_t GetNProtons() const { return fProtons->GetEntriesFast(); }
        Int_t GetNPions() const { return fPions->GetEntriesFast(); }

        void AddLambda(UParticle* particle) {
            new ((*fLambdas)[GetNLambdas()]) UParticle(*particle);
        }
    
        void AddProton(UParticle* particle) {
            new ((*fProtons)[GetNProtons()]) UParticle(*particle);
        }
    
        void AddPion(UParticle* particle) {
            new ((*fPions)[GetNPions()]) UParticle(*particle);
        }

        UParticle *GetLambda(Int_t index) const {
            if (index < 0 || index >= GetNLambdas()){
                std::cout << "OOB\n";
                return nullptr;
            } 
            return static_cast<UParticle*>(fLambdas->At(index));
        }
    
        UParticle *GetProton(Int_t index) const {
            if (index < 0 || index >= GetNProtons()) return nullptr;
            return static_cast<UParticle*>(fProtons->At(index));
        }
    
        UParticle *GetPion(Int_t index) const {
            if (index < 0 || index >= GetNPions()) return nullptr;
            return static_cast<UParticle*>(fPions->At(index));
        }

        void SetPolarizationAxis(){

            TVector3 b(fB, 0, 0); // Impact parameter direction
            TVector3 z(0, 0, 1);   // Beam axis
            polarizationAxis = z.Cross(b.Unit()); // y = z × b̂
            rotation.SetZAxis(polarizationAxis).RotateZ(0);
        }

        void RotateProtons() {
            const Int_t nProtons = fProtons->GetEntriesFast();
            for (Int_t i = 0; i < nProtons; i++) {
                UParticle* proton = static_cast<UParticle*>(fProtons->At(i));
                if (!proton) continue;
    
                TLorentzVector p = proton->GetMomentum();
                p.Transform(rotation);
                proton->SetMomentum(p);
            }
        }

        void RotateProton(UParticle& proton){
            TLorentzVector p = proton.GetMomentum();
            p.Transform(rotation);
            proton.SetMomentum(p);
        }

        void Clear(){
            UEvent::Clear();
            fLambdas->Clear();
            fProtons->Clear();
            fPions->Clear();
        }

        ClassDef(UUEvent, 3); // Version number
    };


void get_particles(const char* filename, Int_t fPdg);
void simulate_lambda_decays(TString inputFile, TString outputFile, Int_t flag = 1, Int_t enhanceStat = 1);
void calc_global_polarization(TString InFileName, TString OutFileName);

Double_t GetRandomValue   (Double_t fMean, Double_t fSigma);
Int_t    returnNumberOfBin(Double_t fValue, Double_t fMinValue, Double_t fMaxValue, Int_t NBins); 

double generateCosTheta(double alpha) {
    TRandom3* rand = new TRandom3(0);
    double r = rand->Uniform(0, 1); // Random number [0,1]
    delete rand;

    // Solve quadratic equation for inverse CDF:
    // r = 0.5*(1 + x) + alpha/4*(x² - 1)
    // Rearranged to: (alpha/2)x² + x - (- 1 + alpha/2 + 2r) = 0
    double a = alpha/2.0;
    double b = 1.0;
    // double c = -(1.0 + alpha/2.0 - 2.0*r);
    double c = -(-1.0+alpha/2.0 + 2.0*r);

    double discriminant = b*b - 4.0*a*c;
    double x1 = (-b + sqrt(discriminant))/(2.0*a);
    double x2 = (-b - sqrt(discriminant))/(2.0*a);

    // Return the root within [-1,1]
    return (fabs(x1) <= 1.0) ? x1 : x2;
}

Double_t positivePhi(const Double_t& phi){
    Double_t phi_rot = phi;
    if (phi < 0) phi_rot += 2.*TMath::Pi();
    return phi_rot;
}


void simulate_lambda_decays(TString inputFile, TString outputFile, Int_t flag, Int_t enhanceStat) {
    // Open input file
    TFile *inFile = TFile::Open(inputFile, "READ");
    TTree *tree = (TTree*)inFile->Get("events");
    // Create output file and tree
    TFile *outFile;
    TTree *outTree;

    UEvent  *inEvent = nullptr;
    UEvent  *outEvent  = new UEvent ();  // Unigen UEvent format
    UUEvent *outUEvent = new UUEvent();  // Custom class which inherited from UEvent

    tree->SetBranchAddress("event", &inEvent);

    //Histograms
    TH1D *hLambdaPt;
    TH1D *hCosTheta;
    TH1D *hProtonPhi;
    TH2D *hAngVsPt;

    Int_t dummy[2];
    UParticle lambda(1, 1, 1, 1, 1, 1, 1, dummy, 1., 1., 1., 1., 1., 1., 1., 1., 1.); //lambda particle instance
    UParticle proton(1, 1, 1, 1, 1, 1, 1, dummy, 1., 1., 1., 1., 1., 1., 1., 1., 1.); //proton particle instance
    UParticle pion(1, 1, 1, 1, 1, 1, 1, dummy, 1., 1., 1., 1., 1., 1., 1., 1., 1.);   //pion particle instance

    double cos_theta_p = 0;  // Proton angle relative to some axis (for anisotropy)
    Double_t fB;

    if(flag == 1){    // Create output file and tree
        outFile = new TFile(outputFile, "RECREATE");
        outTree = new TTree("decays", "Lambda decay products");

        // Setup branches
        outTree->Branch("event", &outEvent);
        outTree->Branch("user_event", "UUEvent", &outUEvent);

        hLambdaPt = new TH1D("hLambdaPt", "Lambda pT;pT [GeV/c];Counts", 100, 0, 3);
        hCosTheta = new TH1D("hCosTheta", "Proton cos(#theta);cos(#theta);Counts", 100, -1, 1);
        hProtonPhi= new TH1D("hProtonPhi", "Proton (#Delta#varphi);#Delta#varphi;Counts", 100, 0, 2.*TMath::Pi());

        hAngVsPt  = new TH2D("hAngVsPt", "cos(#theta) vs Lambda pT;pT [GeV/c];cos(#theta)", 
                                   50, 0, 2, 50, -1, 1);
        outTree->Branch("cos_theta_p", &cos_theta_p);
    }
    else{ // Update output file and tree
        outFile = TFile::Open(outputFile, "UPDATE");
        outTree  = (TTree*) outFile->Get("decays");
        if (outTree){
            std::cout << "Tree found!\n";
            outTree->SetBranchAddress("event", &outEvent);
            outTree->SetBranchAddress("user_event", &outUEvent); 
        } 

        hLambdaPt = (TH1D*) outFile->Get("hLambdaPt");
        hCosTheta = (TH1D*) outFile->Get("hCosTheta");
        hProtonPhi= (TH1D*) outFile->Get("hProtonPhi");
        hAngVsPt  = (TH2D*) outFile->Get("hAngVsPt");

        // Connect branches
        outTree->SetBranchAddress("cos_theta_p", &cos_theta_p);
    }   

    TRandom3* rand = new TRandom3(0);  
    // Lambda mass and decay parameters
    const double mLambda = 1.115683;  // GeV/c²
    const double mProton = 0.938272;  // GeV/c²
    const double mPion   = 0.139570;  // GeV/c²
    
    const double anisotropy = 0.732;  // Strength of anisotropy (-1 to 1)
    // Process events
    Long64_t nEvents = tree->GetEntries();
    std::cout << "Events: " << nEvents << std::endl;
    int j=0;
    for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {
        tree->GetEntry(iEvent);
        outUEvent->Init(*inEvent);
        outEvent->Clear(); // Clear previous event
        outUEvent->Clear();

        for (Int_t i = 0; i < inEvent->GetNpa(); i++) {
            UParticle* part = inEvent->GetParticle(i);
            if (part->GetPdg() != 3122) continue; // Select Lambdas (PDG code 3122)
            lambda = *part;
            //Adding to custom UEvent class:
            outEvent->AddParticle(*part);
            outUEvent->UEvent::AddParticle(*part);
            outUEvent->AddLambda(part);
            
            // Get Lambda 4-momentum
            TLorentzVector lambda_lab(part->Px(), part->Py(), part->Pz(), part->E());
            hLambdaPt->Fill(lambda_lab.Pt());
            
            // Boost to Lambda rest frame
            TVector3 beta = -lambda_lab.BoostVector();
            TLorentzVector lambda_rest = lambda_lab;
            lambda_rest.Boost(beta);
            
            // Generate decay in rest frame with anisotropy
            double pStar = sqrt((mLambda*mLambda - (mProton + mPion)*(mProton + mPion)) *
                                (mLambda*mLambda - (mProton - mPion)*(mProton - mPion))) /
                                (2*mLambda);
            
            // Generate random angles phi
            double phi = rand->Uniform(0, 2*TMath::Pi());

            double alpha = 0.732;
            double cos_theta = generateCosTheta(alpha);
            double sin_theta = sqrt(1 - cos_theta*cos_theta);
            
            // Store proton angle
            cos_theta_p = cos_theta;
            
            // Proton momentum in Lambda rest frame
            TVector3 p_proton_rest(
                pStar * sin_theta * cos(phi),
                pStar * sin_theta * sin(phi),
                pStar * cos_theta
            );
            
            // Pion momentum (opposite to proton)
            TVector3 p_pion_rest = -p_proton_rest;

            // Create 4-vectors in rest frame
            TLorentzVector proton_rest(p_proton_rest, sqrt(pStar*pStar + mProton*mProton));
            TLorentzVector pion_rest(p_pion_rest, sqrt(pStar*pStar + mPion*mPion));

            TLorentzVector proton_lab_pos = part->GetPosition();
            TLorentzVector pion_lab_pos   = part->GetPosition();
            
            Int_t child_null[2];
            Double_t fWeight = 0;
            Int_t enhancedFlag = 0;

            proton = UParticle(i, 2212, 0, i, i, enhancedFlag, -1, child_null, proton_rest, proton_lab_pos, fWeight); //   UParticle(index, pdg, status, parent, parentDecay, mate, decay, child[2], mom,  pos,weight);
            pion   = UParticle(i, -211, 0, i, i, enhancedFlag, -1, child_null, pion_rest,   pion_lab_pos,   fWeight);

            //Proton rotation
            outUEvent->SetPolarizationAxis();
            outUEvent->RotateProton(proton);
            outUEvent->AddProton(&proton);
            outUEvent->AddPion(&pion);

            // Fill histograms
            hCosTheta->Fill(cos_theta);
            hAngVsPt->Fill(lambda_lab.Pt(), cos_theta);
            hProtonPhi->Fill(positivePhi(proton.GetMomentum().Phi()));
            //hProtonPhi->Fill(proton.GetMomentum().Phi());

            Double_t fEnhanceStat = enhanceStat;
            while(fEnhanceStat > 1){ //enhancing of lambdas
                enhancedFlag = -9;

                TLorentzVector mom_rand( //generating momentum
                    GetRandomValue(part->Px(), 0.03),
                    GetRandomValue(part->Py(), 0.03),
                    GetRandomValue(part->Pz(), 0.03),
                    GetRandomValue(part->E(), 0.03)
                );
                
                TLorentzVector pos_rand( //generating position
                    GetRandomValue(part->X(), 0.03),
                    GetRandomValue(part->Y(), 0.03),
                    GetRandomValue(part->Z(), 0.03),
                    GetRandomValue(part->T(), 0.03)
                );

                UParticle enhancedLambda(*part);
                enhancedLambda.SetPosition(pos_rand);
                enhancedLambda.SetMate(enhancedFlag);

                // Process decay. Same logic
                TLorentzVector lambda_lab = enhancedLambda.GetMomentum();
                hLambdaPt->Fill(lambda_lab.Pt());
                
                TVector3 beta = -lambda_lab.BoostVector();
                double pStar = sqrt((mLambda*mLambda - (mProton + mPion)*(mProton + mPion)) *
                                 (mLambda*mLambda - (mProton - mPion)*(mProton - mPion))) /
                                 (2*mLambda);
                
                double phi = rand->Uniform(0, 2*TMath::Pi());
                double cos_theta = generateCosTheta(anisotropy);
                double sin_theta = sqrt(1 - cos_theta*cos_theta);
                
                // Create decay products
                TVector3 p_proton_rest(
                    pStar * sin_theta * cos(phi),
                    pStar * sin_theta * sin(phi),
                    pStar * cos_theta
                );
                TVector3 p_pion_rest = -p_proton_rest;

                TLorentzVector proton_rest(p_proton_rest, sqrt(pStar*pStar + mProton*mProton));
                TLorentzVector pion_rest(p_pion_rest, sqrt(pStar*pStar + mPion*mPion));

                // Create final particles
                Int_t child_null[2] = {0, 0};
                UParticle protonEnhanced(
                    i, 2212, 0, i, i, enhancedFlag, -1,
                    child_null, proton_rest, enhancedLambda.GetPosition(), 0
                );
                UParticle pionEnhanced(
                    i, -211, 0, i, i, enhancedFlag, -1,
                    child_null, pion_rest, enhancedLambda.GetPosition(), 0
                );

                outUEvent->SetPolarizationAxis();
                outUEvent->RotateProton(protonEnhanced);

                // Add to event
                outUEvent->AddLambda(&enhancedLambda);
                outUEvent->AddProton(&protonEnhanced);
                outUEvent->AddPion(&pionEnhanced);

                // Fill histograms
                hCosTheta->Fill(cos_theta);
                hAngVsPt->Fill(lambda_lab.Pt(), cos_theta);
                hProtonPhi->Fill(positivePhi(protonEnhanced.GetMomentum().Phi()));

                fEnhanceStat--;
            }

        }
        outTree->Fill();
    }

    outFile->cd();
    outTree->Write("",TObject::kOverwrite);
    hLambdaPt->Write("",TObject::kOverwrite);
    hCosTheta->Write("",TObject::kOverwrite);
    hProtonPhi->Write("",TObject::kOverwrite);
    hAngVsPt->Write("",TObject::kOverwrite);

    // Plot results
    TCanvas *c1 = new TCanvas("c1", "Lambda_Decays", 1200, 800);
    c1->Divide(2,2);
    c1->cd(1); hLambdaPt->Draw();
    c1->cd(2); hCosTheta->Draw();
    c1->cd(3); hAngVsPt->Draw("colz");
    c1->cd(4); hProtonPhi->Draw();

    c1->SaveAs("lambda_decays_plots.png");

    // Cleanup
    delete c1;
    delete outUEvent;
    delete outEvent;
    delete rand;

    outFile->Close();
    inFile->Close();
}



void calc_global_polarization(TString InFileName, TString OutFileName){
    gStyle->SetOptStat(0);  // Disable stats globally
    TF1* fitPhiDistro = new TF1("fitPhiDistro", "[0]*(1+2*[1]*TMath::Sin(x)+2*[2]*TMath::Cos(x))", 0, 2*TMath::Pi()); //Fitting function
    const double anisotropy = 0.732;  // Strength of anisotropy (-1 to 1)

    //--Parameters for binning
    Int_t NpTBins = 10;//number of pT bins
    Int_t NYBins = 15; //number of rapidity bins
    Double_t fpTMax = 0.9; // max pT [GeV/c]
    Double_t fpTMin = 0.2; // min pT [GeV/c]
    Double_t fYMax =  1.2; // max rapidity
    Double_t fYMin = -1.2; // min rapidity

    //--Histograms--//
    TH1D* hProtonLambdaFrame_phi = new TH1D("hProtonLambdaFrame_phi", " ; #Delta#phi, rad;Counts", 50,  0, 2*TMath::Pi()); //proton phi* distribution (lambda frame)
    TH1D* hProtonLabFrame_pT = new TH1D("hProtonLabFrame_pT", " proton lab frame ; p_{T}, GeV/c; Counts", 50,  0, 2.); //proton pT distribution (lab frame)
    TH1D* hProtonLabFrame_Y  = new TH1D("hProtonLabFrame_Y", " proton lab frame ; Y; Counts", 50,  -2., 2.); //proton Y distribution (lab frame)
    TH1D* hProtonLambdaFrame_phi_pTBins[NpTBins]; //pT binnig
    TH1D* hProtonLambdaFrame_phi_YBins [NYBins ]; //rapidity binnig
    //loop for TH1D object instance 
    for(Int_t iHisto = 0; iHisto < NpTBins; iHisto++) hProtonLambdaFrame_phi_pTBins[iHisto] = new TH1D(TString::Format("hProtonLambdaFrame_phi_pT_%i",iHisto), " ; #Delta#phi, rad;Counts", 50,  0, 2*TMath::Pi());
    for(Int_t iHisto = 0; iHisto < NYBins;  iHisto++) hProtonLambdaFrame_phi_YBins [iHisto] = new TH1D(TString::Format("hProtonLambdaFrame_phi_Y_%i",iHisto), " ; #Delta#phi, rad;Counts", 50,  0, 2*TMath::Pi());

    Int_t dummy[2];
    UParticle *lambda = new UParticle(1, 1, 1, 1, 1, 1, 1, dummy, 1., 1., 1., 1., 1., 1., 1., 1., 1.); //lambda particle instance
    UParticle *proton = new UParticle(1, 1, 1, 1, 1, 1, 1, dummy, 1., 1., 1., 1., 1., 1., 1., 1., 1.); //proton particle instance
    UParticle *pion   = new UParticle(1, 1, 1, 1, 1, 1, 1, dummy, 1., 1., 1., 1., 1., 1., 1., 1., 1.); //pion particle instance

    
    TFile* inFile  = TFile::Open(InFileName, "READ"); //Open input file with polarized proton
    TFile *outFile = TFile::Open(OutFileName, "RECREATE");//Resulting out file

    TTree *inTree = (TTree*)inFile->Get("decays");
    // Setup branches
    UUEvent *inUEvent = new UUEvent();
    inTree->SetBranchAddress("user_event", &inUEvent);

    //Loop over events
    Long64_t nEvents = inTree->GetEntries();
    for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {
        inUEvent->Clear();
        inTree->GetEntry(iEvent);
        //Loop over particles within event
        for (Int_t iPart = 0; iPart < inUEvent->GetNLambdas();iPart++) {
            if(!inUEvent->GetLambda(iPart)) continue;
            lambda = inUEvent->GetLambda(iPart); //lambda in laboratory frame
            proton = inUEvent->GetProton(iPart); //proton in lambda frame
            pion   = inUEvent->GetPion(iPart);   //pion   in lambda frame

            Double_t phiStar = positivePhi(proton->GetMomentum().Phi()); //Get phi*
            hProtonLambdaFrame_phi->Fill(phiStar); //phi* distribution 

            TLorentzVector lambda_lab = lambda->GetMomentum();
            TVector3 beta = -lambda_lab.BoostVector(); //Get lorentz boost

            TLorentzVector proton_lab = proton->GetMomentum(); // getting laboratory proton for further properly binning over [pT;Y]
            proton_lab.Boost(beta);                             //boost proton to lab frame

            //if bin isn't uderflow or overflow then put phi*-distribution in histogram[number of bin]
            if(returnNumberOfBin(proton_lab.Pt(),       fpTMin, fpTMax, NpTBins) > -1) hProtonLambdaFrame_phi_pTBins[(returnNumberOfBin(proton_lab.Pt(),       fpTMin, fpTMax, NpTBins))]->Fill(phiStar); 
            if(returnNumberOfBin(proton_lab.Rapidity(), fYMin,  fYMax,  NYBins ) > -1) hProtonLambdaFrame_phi_YBins [(returnNumberOfBin(proton_lab.Rapidity(), fYMin,  fYMax,  NYBins ))]->Fill(phiStar); 

            //Filling pT and Y histograms
            hProtonLabFrame_pT->Fill(proton_lab.Pt());
            hProtonLabFrame_Y->Fill(proton_lab.Rapidity());
        }

    }
    //--pT binning--//
    TH1D* hphi_pTBins= new TH1D("hphi_pTBins", " p_{T} binning; p_{T}, GeV/c; P_{#lambda}", NpTBins,  fpTMin, fpTMax); // Histogram with pT bins
    for(Int_t iHisto = 0; iHisto < NpTBins; iHisto++) { //Loop over pT bins
        hProtonLambdaFrame_phi_pTBins[iHisto]->Fit(fitPhiDistro);
        hProtonLambdaFrame_phi_pTBins[iHisto]->Write();
        hphi_pTBins->SetBinContent(iHisto+1, (fitPhiDistro->GetParameter(1)*8)/(TMath::Pi()*anisotropy));
        hphi_pTBins->SetBinError  (iHisto+1, (fitPhiDistro->GetParError(1)*8)/(TMath::Pi()*anisotropy));
    }

    hphi_pTBins->SetMarkerStyle(23); // Filled circle (ROOT default)
    hphi_pTBins->SetMarkerSize(1.2); // Slightly larger markers
    hphi_pTBins->Write();

    //--Y binning--//
    //same logic like in previous one
    TH1D* hphi_YBins= new TH1D("hphi_YBins", " rapidity binning; Y; P_{#lambda}", NYBins,  fYMin, fYMax);
    for(Int_t iHisto = 0; iHisto < NYBins; iHisto++) {
        hProtonLambdaFrame_phi_YBins[iHisto]->Fit(fitPhiDistro);
        hProtonLambdaFrame_phi_YBins[iHisto]->Write();
        hphi_YBins->SetBinContent(iHisto+1, (fitPhiDistro->GetParameter(1)*8)/(TMath::Pi()*anisotropy));
        hphi_YBins->SetBinError(iHisto+1, (fitPhiDistro->GetParError(1)*8)/(TMath::Pi()*anisotropy));
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
    lProtonLambdaFrame_phi->AddEntry("", Form("P_{#Lambda} = %.3f #pm %.3f", (fitPhiDistro->GetParameter(1)*8)/(TMath::Pi()*anisotropy), (fitPhiDistro->GetParError(1)*8)/(TMath::Pi()*anisotropy)), "");
    lProtonLambdaFrame_phi->AddEntry("", Form("P_{#Lambda} = (8p_{1}) / (%.3f#pi)", anisotropy), "");

    // Plot results
    TCanvas *c1 = new TCanvas("c1", "proton", 1200, 800);
    c1->Divide(2,2);
    c1->cd(1); hProtonLambdaFrame_phi->Draw(); lProtonLambdaFrame_phi->Draw("same");
    c1->cd(2); hProtonLabFrame_pT->Draw();
    c1->cd(3); hphi_pTBins->Draw("PE");
    c1->cd(4); hphi_YBins->Draw("PE");

    c1->SaveAs("proton_polarization_plots.png");
    c1->Write();
    hProtonLabFrame_pT->Write();
    hProtonLabFrame_Y->Write();
    hProtonLambdaFrame_phi->Write();
    outFile->Close();
    inFile->Close();
}

Double_t GetRandomValue(Double_t fMean, Double_t fSigma)
{
    TRandom3* rand = new TRandom3(0);  // Seed with 0 for reproducibility
    Double_t randomVal = rand->Gaus(fMean, fSigma);
    delete rand;
    return randomVal;
}


Int_t returnNumberOfBin(Double_t fValue, Double_t fMinValue, Double_t fMaxValue, Int_t NBins) { //Function return number of bin
    if (NBins <= 0) return -1;
    if (fValue < fMinValue) return -1;  // Underflow
    if (fValue >= fMaxValue) return -1; // Overflow
    
    Double_t binWidth = (fMaxValue - fMinValue)/NBins;
    //std::cout << "bin_value" <<  (fValue - fMinValue)/binWidth << std::endl;
    Int_t bin = static_cast<Int_t>((fValue - fMinValue)/binWidth);
    
    // Protect against edge case where fValue == fMaxValue
    return (bin < NBins) ? bin : NBins-1;
}

