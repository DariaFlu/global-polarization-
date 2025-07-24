#ifdef __CLING__
#pragma link C++ class LambdaAnalysisData+;
#pragma link C++ class UUEvent+;
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
#include "UEvent.h"
#include "UParticle.h"

void get_particles(const char* filename, Int_t fPdg);
void simulate_lambda_decays(TString inputFile, TString outputFile, TString confInFile, Int_t flag = 1, Int_t enhanceStat = 1);
void calc_global_polarization(TString InFileName, TString OutFileName);
Double_t Lambda_v1_v2_P_parameterization(TFile* Lambda_yield, Double_t fBVal, UParticle &ULambda); //Valeriy's function for properly lambda generation

Double_t GetRandomValue   (Double_t fMean, Double_t fSigma);
Int_t    returnNumberOfBin(Double_t fValue, Double_t fMinValue, Double_t fMaxValue, Int_t NBins); 


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

        Double_t fSigmaPol;

        Double_t fpolx  = 0; 
        Double_t fpolSx = 0.07;  
        Double_t fpoly  = 0;
        Double_t fpolSy = 0.07; 
        Double_t fpolz  = 0;
        Double_t fpolSz = 0.07;


        //TClonesArray *fParticles;

        TClonesArray *fLambdas = new TClonesArray("UParticle");
        TClonesArray *fProtons = new TClonesArray("UParticle");
        TClonesArray *fPions = new TClonesArray("UParticle");

        TClonesArray* fPolarizations; // Stores polarization vectors

        TVector3 polarizationAxis;
        TVector3 polarizationVec;

        TRotation rotation;
    public :
        UUEvent() {
            fLambdas = new TClonesArray("UParticle");
            fProtons = new TClonesArray("UParticle");
            fPions   = new TClonesArray("UParticle");

            fPolarizations = new TClonesArray("TVector3");  // Add this line
        }
        // Copy constructor
         UUEvent(const UUEvent& right) {
            fLambdas = new TClonesArray("UParticle",100);
            fProtons = new TClonesArray("UParticle",100);
            fPions = new TClonesArray("UParticle",100);

            fPolarizations = new TClonesArray("TVector3", 100);

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

            delete fPolarizations;
            // delete fParticles;
        }
        Int_t GetNLambdas() const { return fLambdas->GetEntriesFast(); }
        Int_t GetNProtons() const { return fProtons->GetEntriesFast(); }
        Int_t GetNPions() const { return fPions->GetEntriesFast(); }
        Int_t GetNPol  () const { return fPolarizations->GetEntriesFast(); }

        Double_t GetBimp(){ return fB; }

        void SetPolarizationSigma(Double_t _fSigmaPol) { fSigmaPol = _fSigmaPol; }

        void AddPolarization(const TVector3& pol) {
            // std::cout<<"Polarization has been added to UUEvent"<< std::endl;
            new ((*fPolarizations)[GetNPol()]) TVector3(pol);
        }

        void AddLambda(UParticle* particle) {
            new ((*fLambdas)[GetNLambdas()]) UParticle(*particle);
        }

        void AddProton(UParticle* particle) {
            new ((*fProtons)[GetNProtons()]) UParticle(*particle);
        }
    
        void AddPion(UParticle* particle) {
            new ((*fPions)[GetNPions()]) UParticle(*particle);
        }

        TVector3 GetPolarization(Int_t index) const {
            if (index < 0 || index >= GetNLambdas()) return TVector3(0, 0, 0);
            // return *(TVector3*)(fPolarizations->At(index));
            return *static_cast<TVector3*>(fPolarizations->At(index));
        }

        UParticle *GetLambda(Int_t index) const {
            if (index < 0 || index >= GetNLambdas()){
                // std::cout << "OOB\n";
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

        TVector3 GetPolLambda(UParticle& lambda, Double_t fpoly, Double_t _fSigmaPol = 0.3){ 
            //std::cout << 1 << std::endl;
            SetPolarizationSigma(_fSigmaPol);
            
            Double_t polX = gRandom->Gaus(fpolx,fpolSx);  // generate polarization direction
            Double_t polY = gRandom->Gaus(fpoly,fpolSy);  // generate polarization direction /100
            Double_t polZ = gRandom->Gaus(fpolz,fpolSz);  // generate polarization direction
            polarizationVec = TVector3(polX, polY, polZ);

            // std::cout<<"TVector3(polX, polY, polZ) "<<polX<<" "<<polY<<" "<<polZ<<std::endl;

            ROOT::Math::RotationZ rotateRP(fPhi);  // set rotation transformation
            polarizationVec = rotateRP*polarizationVec;  // rotate reaction plane

            Float_t polmag = TMath::Sqrt(polarizationVec.Mag());
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
                  x1 = 1-2*gRandom->Rndm();
                  x2 = 1-2*gRandom->Rndm();
                  R2 = x1*x1+x2*x2;
               } while(R2 >= 1);
               R=2*TMath::Sqrt(1-R2);
               polarizationVec.SetXYZ(x1*R, x2*R, 1-2*R2);  // random unitary vector - background for signal
            }
   
            Float_t xxx = gRandom->Rndm();
            if (xxx < 1./2.*(1.-polmag)) {  // probability of spin flip
                polarizationVec *= -1.;  // spin flip according to mean polarization
            }

            // std::cout<<"polmag =  "<<polmag<<std::endl;

            // std::cout<<"Polarization generated (X,Y,Z) "<<polarizationVec.X()<<" "<<polarizationVec.Y()<<" "<<polarizationVec.Z()<<std::endl;
            AddPolarization(polarizationVec);  
            return polarizationVec;
        }


        void Clear(){
            UEvent::Clear();
            fLambdas->Clear();
            fProtons->Clear();
            fPions->Clear();
            fPolarizations->Clear();
        }

        ClassDef(UUEvent, 3); // Version number
    };



double generateCosTheta(double alpha, double pol = 0.6) {

    // TF1* randomCosTheta = new TF1("randomCosTheta", "1+(TMath::Pi()*[0]*[1])/8*x", -1.0, 1.0);
    TF1* randomCosTheta = new TF1("randomCosTheta", "1+[0]*[1]*x", -1.0, 1.0);
    // Set the parameters
    // randomCosTheta->SetParameters(alpha);

    randomCosTheta->SetParameters(alpha, pol);
    // Generate a random value from the distribution
    double value = randomCosTheta->GetRandom();
    // std::cout<<"random costh = "<<value<<std::endl;
    // Clean up
    delete randomCosTheta;
    
    return value;

}

Double_t positivePhi(const Double_t& phi){
    Double_t phi_rot = phi;
    if (phi < 0) phi_rot += 2.*TMath::Pi();
    return phi_rot;
}


void simulate_lambda_decays(TString inputFile, TString outputFile, TString confInFile, Int_t flag, Int_t enhanceStat) {
    // Open input file
    TFile *inFile = TFile::Open(inputFile, "READ");
    TTree *tree = (TTree*)inFile->Get("events");
    // Create output file and tree
    TFile *outFile;
    TTree *outTree;

    TFile* Lambda_yield = TFile::Open(confInFile,"READ"); // File with TH2F pt-y

    UEvent  *inEvent = nullptr;
    UEvent  *outEvent  = new UEvent ();  // Unigen UEvent format
    UUEvent *outUEvent = new UUEvent();  // Custom class which inherited from UEvent

    tree->SetBranchAddress("event", &inEvent);

    //Histograms
    TH1D *hLambdaPt;
    TH1D *hCosTheta;
    TH1D *hLambdaPhi;
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

        hLambdaPt = new TH1D("hLambdaPt", "Lambda pT;pT [GeV/c];Counts", 100, 0, 2);
        hCosTheta = new TH1D("hCosTheta", "Proton cos(#theta);cos(#theta);Counts", 100, -1, 1);
        hLambdaPhi= new TH1D("hLambdaPhi", "Lambda (#Delta#varphi);#Delta#varphi;Counts", 100, 0, 2.*TMath::Pi());

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
        hLambdaPhi= (TH1D*) outFile->Get("hLambdaPhi");
        hAngVsPt  = (TH2D*) outFile->Get("hAngVsPt");

        // Connect branches
        outTree->SetBranchAddress("cos_theta_p", &cos_theta_p);
    }   

    TRandom3* rand = new TRandom3(0);  
    // Lambda mass and decay parameters
    const double mLambda = 1.115683;  // GeV/c²
    const double mProton = 0.938272;  // GeV/c²
    const double mPion   = 0.139570;  // GeV/c²
    
    Double_t fSigmaPolVal = 1.5;

    const double anisotropy = 0.732;  // Strength of anisotropy 
    const double polPercent = 0.6;  // Polarization percent

    Double_t fBMin = 3.44;
    Double_t fBMax = 7.44;
    Int_t child_null[2];

    // Process events
    Long64_t nEvents = tree->GetEntries();
    std::cout << "Events: " << nEvents << std::endl;
    for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {
        tree->GetEntry(iEvent);
        outUEvent->Init(*inEvent);
        outEvent->Clear(); // Clear previous event
        outUEvent->Clear();
        Int_t lambdaCounter = 0;
        for (Int_t i = 0; i < inEvent->GetNpa(); i++) {
            UParticle* part = inEvent->GetParticle(i);
            if(i == inEvent->GetNpa()-1 && lambdaCounter == 0){
                Double_t pt  = gRandom->Uniform(0.5, 3.0);
                Double_t phi = gRandom->Uniform(0, 2*TMath::Pi());
                Double_t eta = gRandom->Uniform(-1, 1);

                TLorentzVector newLambdaPos(
                    GetRandomValue(0, 0.03),
                    GetRandomValue(0, 0.03),
                    GetRandomValue(0, 0.03),
                    GetRandomValue(0, 0.03));
                TLorentzVector newLambdaMom(
                    GetRandomValue(pt*cos(phi), 0.03),
                    GetRandomValue(pt*sin(phi), 0.03),
                    GetRandomValue(pt*sinh(eta), 0.03),
                    GetRandomValue(sqrt(pt*pt*cosh(eta)*cosh(eta) + mLambda*mLambda), 0.03)
                );
                part = new UParticle(i, 3122, 1, 1, 1, -15, -1,
                                    child_null, newLambdaMom, newLambdaPos, 0 );

            }
            if (part->GetPdg() != 3122) continue; // Select Lambdas (PDG code 3122)
            lambdaCounter++;
            // if (inEvent->GetB() < fBMin || inEvent->GetB() > fBMax)  continue;
            lambda = *part;
            //Adding to custom UEvent class:
            outEvent->AddParticle(*part);
            outUEvent->UEvent::AddParticle(*part);
            Double_t fPolY = Lambda_v1_v2_P_parameterization(Lambda_yield, inEvent->GetB(), lambda);

            outUEvent->AddLambda(&lambda);

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

            TVector3 pol = outUEvent->GetPolLambda(lambda, fPolY/100., fSigmaPolVal);
            outUEvent->AddPolarization(pol);
            // std::cout<<"PolVec for proton gen (Mag, X, Y, Z) "<<pol.Mag()<<" "<< pol.X()<<" "<< pol.Y()<<" "<<pol.Z()<<std::endl;

            double cos_theta = generateCosTheta(0.732, pol.Mag());
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

            TVector3 unit = TVector3(
                sin_theta * cos(phi),
                sin_theta * sin(phi),
                cos_theta
            );

            unit.RotateUz(pol.Unit()); // rotate Z to norm (with extra phi rotation which is random)

            phi = unit.Phi();
            cos_theta = TMath::Cos(unit.Theta());
            sin_theta = TMath::Sin(unit.Theta());

            if (TMath::Abs(cos_theta) >= 1.0) cos_theta = TMath::Sign(1.0, cos_theta);
            else sin_theta = TMath::Sqrt((1. - cos_theta) * (1. + cos_theta));
      
            // Generate proton momentum in lambda rest frame (aligned with z-axis)
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

            TLorentzVector proton_lab_pos = part->GetPosition();
            TLorentzVector pion_lab_pos   = part->GetPosition();
            
            Double_t fWeight = 0;
            Int_t enhancedFlag = 0;

            proton = UParticle(i, 2212, 0, i, i, enhancedFlag, -1, child_null, proton_rest, proton_lab_pos, fWeight); //   UParticle(index, pdg, status, parent, parentDecay, mate, decay, child[2], mom,  pos,weight);
            pion   = UParticle(i, -211, 0, i, i, enhancedFlag, -1, child_null, pion_rest,   pion_lab_pos,   fWeight);

            //Proton rotation
            // outUEvent->SetPolarizationAxis();
            // outUEvent->RotateProton(proton);

            outUEvent->AddProton(&proton);
            outUEvent->AddPion(&pion);

            // Fill histograms
            hCosTheta->Fill(cos_theta);
            hAngVsPt->Fill(lambda_lab.Pt(), cos_theta);
            // hLambdaPhi->Fill(positivePhi(proton.GetMomentum().Phi()));
            hLambdaPhi->Fill(positivePhi(lambda.GetMomentum().Phi()));

            Double_t fEnhanceStat = enhanceStat;
            while(fEnhanceStat > 1){ //enhancing of lambdas
                enhancedFlag = -9;

                //add random function
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

                fPolY = Lambda_v1_v2_P_parameterization(Lambda_yield, inEvent->GetB(), enhancedLambda); 

                // Process decay. Same logic
                TLorentzVector lambda_lab = enhancedLambda.GetMomentum();
                hLambdaPt->Fill(lambda_lab.Pt());
                
                TVector3 beta = -lambda_lab.BoostVector();
                double pStar = sqrt((mLambda*mLambda - (mProton + mPion)*(mProton + mPion)) *
                                 (mLambda*mLambda - (mProton - mPion)*(mProton - mPion))) /
                                 (2*mLambda);
                
                double phi = rand->Uniform(0, 2*TMath::Pi());

                outUEvent->AddLambda(&enhancedLambda);

                TVector3 pol = outUEvent->GetPolLambda(lambda, fPolY/100., fSigmaPolVal);
                outUEvent->AddPolarization(pol);
                // std::cout<<"PolVec for proton gen (Mag, X, Y, Z) "<<pol.Mag()<<" "<< pol.X()<<" "<< pol.Y()<<" "<<pol.Z()<<std::endl;
    
                double cos_theta = generateCosTheta(0.732, pol.Mag());
                double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
    
                TVector3 unit = TVector3(
                    sin_theta * cos(phi),
                    sin_theta * sin(phi),
                    cos_theta
                );
    
                unit.RotateUz(pol.Unit()); // rotate Z to norm (with extra phi rotation which is random)
    
                phi = unit.Phi();
                cos_theta = TMath::Cos(unit.Theta());
                sin_theta = TMath::Sin(unit.Theta());
    
                if (TMath::Abs(cos_theta) >= 1.0) cos_theta = TMath::Sign(1.0, cos_theta);
                else sin_theta = TMath::Sqrt((1. - cos_theta) * (1. + cos_theta));
          
                // Generate proton momentum in lambda rest frame (aligned with z-axis)
                TVector3 p_proton_rest = pStar * TVector3(
                    sin_theta * cos(phi),
                    sin_theta * sin(phi),
                    cos_theta
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

                // outUEvent->SetPolarizationAxis();
                // outUEvent->RotateProton(protonEnhanced);

                // Add to event
                outUEvent->AddProton(&protonEnhanced);
                outUEvent->AddPion(&pionEnhanced);

                // Fill histograms
                hCosTheta->Fill(cos_theta);
                hAngVsPt->Fill(lambda_lab.Pt(), cos_theta);
                // hLambdaPhi->Fill(positivePhi(protonEnhanced.GetMomentum().Phi()));
                hLambdaPhi->Fill(positivePhi(enhancedLambda.GetMomentum().Phi()));

                fEnhanceStat--;
            }

        }
        outTree->Fill();
    }

    outFile->cd();
    outTree->Write("",TObject::kOverwrite);
    hLambdaPt->Write("",TObject::kOverwrite);
    hCosTheta->Write("",TObject::kOverwrite);
    hLambdaPhi->Write("",TObject::kOverwrite);
    hAngVsPt->Write("",TObject::kOverwrite);

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
    delete outUEvent;
    delete outEvent;
    delete rand;

    outFile->Close();
    inFile->Close();
    Lambda_yield->Close(); //proper file operation

}



void calc_global_polarization(TString InFileName, TString OutFileName){
    gStyle->SetOptStat(0);  // Disable stats globally
    TF1* fitPhiDistro = new TF1("fitPhiDistro", "[0]*(1+2*[1]*TMath::Sin(x)+2*[2]*TMath::Cos(x))", 0, 2*TMath::Pi()); //Fitting function
    // TF1* fitPhiDistro = new TF1("fitPhiDistro", 
    //     "[0]*(1 + [1]*sin(x) + [2]*cos(x) + [3]*sin(2*x) + [4]*cos(2*x))", 
    //     0, 2*TMath::Pi());
    const double anisotropy = 0.732;  // Strength of anisotropy (-1 to 1)

    TF1* fitPhiCollectiveFlowDistro = new TF1("fitPhiCollectiveFlowDistro", "[0]*(1+2*[1]*TMath::Cos(2*x))", 0, 2*TMath::Pi()); //Fitting function
    // const double anisotropy = 0.732;  // Strength of anisotropy (-1 to 1)

    //--Parameters for binning
    Int_t NpTBins = 11;//number of pT bins
    Int_t NYBins = 11; //number of rapidity bins
    Double_t fpTMax = 1.6; // max pT [GeV/c]
    Double_t fpTMin = 0.2; // min pT [GeV/c]
    Double_t fYMax =  0.5; // max rapidity
    Double_t fYMin = -0.5; // min rapidity

    //--Histograms--//
    TH1D* hProtonLambdaFrame_phi = new TH1D("hProtonLambdaFrame_phi", " ; #Delta#phi, rad;Counts", 50,  0, 2*TMath::Pi()); //proton phi* distribution (lambda frame)
    TH1D* hLambdaLabFrame_pT = new TH1D("hLambdaLabFrame_pT", " lambda lab frame ; p_{T}, GeV/c; Counts", 50,  0., 2.0); //proton pT distribution (lab frame)
    TH1D* hlambdaLabFrame_Y  = new TH1D("hlambdaLabFrame_Y", " lambda lab frame ; Y; Counts", 50,  -1.5, 1.5); //lambda Y distribution (lab frame)
    TH1D* hlambdaLabFrame_phiDistr  = new TH1D("hlambdaLabFrame_phiDistr", " lambda lab frame ; #Delta#phi; Counts", 100, 0, 2*TMath::Pi()); //lambda Y distribution (lab frame)
    TProfile* fProfV2pT = new TProfile("fProfV2pT", "v_{2} vs pT", 20, 0, 2, -10, 10);
    TProfile* fProfV1Y = new TProfile("fProfV1Y", "v_{1} vs Y", 20, -1, 1, -2, 2);

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
    Double_t fBMin = 3.44;
    Double_t fBMax = 6.64;
    //Loop over events
    Long64_t nEvents = inTree->GetEntries();
    for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {
        inUEvent->Clear();
        inTree->GetEntry(iEvent);
        //Loop over particles within event
        for (Int_t iPart = 0; iPart < inUEvent->GetNLambdas();iPart++) {
            if(inUEvent->GetBimp() < fBMin || inUEvent->GetBimp() > fBMax) continue;
            if(!inUEvent->GetLambda(iPart)) continue;
            lambda = inUEvent->GetLambda(iPart); //lambda in laboratory frame
            proton = inUEvent->GetProton(iPart); //proton in lambda frame
            pion   = inUEvent->GetPion(iPart);   //pion   in lambda frame

            Double_t phiStar = positivePhi(proton->GetMomentum().Phi()); //Get phi*
            hProtonLambdaFrame_phi->Fill(phiStar); //phi* distribution 

            TLorentzVector lambda_lab = lambda->GetMomentum();
            TVector3 beta = -lambda_lab.BoostVector(); //Get lorentz boost

            TLorentzVector lambda_rest = lambda_lab;

            lambda_rest.Boost(beta);

            //if bin isn't uderflow or overflow then put phi*-distribution in histogram[number of bin]
            if(returnNumberOfBin(lambda_lab.Pt(),       fpTMin, fpTMax, NpTBins) > -1) hProtonLambdaFrame_phi_pTBins[(returnNumberOfBin(lambda_lab.Pt(),       fpTMin, fpTMax, NpTBins))]->Fill(phiStar); 
            if(returnNumberOfBin(lambda_lab.Rapidity(), fYMin,  fYMax,  NYBins ) > -1) hProtonLambdaFrame_phi_YBins [(returnNumberOfBin(lambda_lab.Rapidity(), fYMin,  fYMax,  NYBins ))]->Fill(phiStar); 


            //Filling pT and Y histograms
            hLambdaLabFrame_pT->Fill(lambda_lab.Pt());
            hlambdaLabFrame_Y ->Fill(lambda_lab.Rapidity());
            hlambdaLabFrame_phiDistr->Fill(positivePhi(lambda->GetMomentum().Phi()));
            // if(lambda->GetMomentum().Rapidity() > -0.75 && lambda->GetMomentum().Rapidity() > -0.75 ) hlambdaLabFrame_phiDistrBin[(returnNumberOfBin(lambda_lab.Pt(), fpTMin, fpTMax, NpTBins))]->Fill(phiStar); 
            //std::cout << "Cos[2 phi]: " << TMath::Cos(2*positivePhi(lambda->GetMomentum().Phi())) << std::endl;
            Double_t value = (TMath::Power(lambda_lab.Px(), 2) - TMath::Power(lambda_lab.Py(), 2))/TMath::Power(lambda_lab.Pt(), 2);
            Double_t valPhi = TMath::Cos(2*(lambda_lab.Phi()));
            //if (lambda->GetMomentum().Rapidity() > -0.75 && lambda->GetMomentum().Rapidity() < 0.75 ) fProfV2pT->Fill(lambda_lab.Pt(), TMath::Cos(2*positivePhi(lambda->GetMomentum().Phi())));
            if (lambda->GetMomentum().Rapidity() > -0.75 && lambda->GetMomentum().Rapidity() < 0.75 ){
                fProfV2pT->Fill(lambda_lab.Pt(), valPhi);
                std::cout << fProfV2pT->GetYmin() << " " << fProfV2pT->GetYmax() << " " << TMath::Cos(2*(lambda_lab.Phi())) << std::endl;
            }

            if (lambda->GetMomentum().Pt() > 0.4 && lambda->GetMomentum().Pt() < 2.) fProfV1Y->Fill(lambda_lab.Rapidity(),TMath::Cos(positivePhi(lambda->GetMomentum().Phi())));
            



        }

    }
    //--pT binning--//
    TH1D* hphi_pTBins= new TH1D("hphi_pTBins", " p_{T} binning; p_{T}, GeV/c; P_{#lambda} %", NpTBins,  fpTMin, fpTMax); // Histogram with pT bins
    for(Int_t iHisto = 0; iHisto < NpTBins; iHisto++) { //Loop over pT bins
        hProtonLambdaFrame_phi_pTBins[iHisto]->Fit(fitPhiDistro);
        hProtonLambdaFrame_phi_pTBins[iHisto]->Write();
        hphi_pTBins->SetBinContent(iHisto+1, (fitPhiDistro->GetParameter(1)*8)/(TMath::Pi()*anisotropy)*100);
        hphi_pTBins->SetBinError  (iHisto+1, (fitPhiDistro->GetParError(1)*8)/(TMath::Pi()*anisotropy)*100);
    }

    hphi_pTBins->SetMarkerStyle(23); // Filled circle (ROOT default)
    hphi_pTBins->SetMarkerSize(1.2); // Slightly larger markers
    hphi_pTBins->Write();

    //--Y binning--//
    //same logic like in previous one
    TH1D* hphi_YBins= new TH1D("hphi_YBins", " rapidity binning; Y; P_{#lambda} %", NYBins,  fYMin, fYMax);
    for(Int_t iHisto = 0; iHisto < NYBins; iHisto++) {
        hProtonLambdaFrame_phi_YBins[iHisto]->Fit(fitPhiDistro);
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
    lProtonLambdaFrame_phi->AddEntry("", Form("P_{#Lambda} = %.0f #pm %.0f %%", (fitPhiDistro->GetParameter(1)*8)/(TMath::Pi()*anisotropy)*100, (fitPhiDistro->GetParError(1)*8)/(TMath::Pi()*anisotropy)*100), "");
    lProtonLambdaFrame_phi->AddEntry("", Form("P_{#Lambda} = (8p_{1}) / (%.3f#pi)", anisotropy), "");

    // Plot results
    TCanvas *c1 = new TCanvas("c1", "proton", 1200, 800);
    c1->Divide(2,2);
    c1->cd(1); hProtonLambdaFrame_phi->Draw(); lProtonLambdaFrame_phi->Draw("same");
    c1->cd(2); hLambdaLabFrame_pT->Draw();
    c1->cd(3); hphi_pTBins->Draw("PE");
    c1->cd(4); hphi_YBins->Draw("PE");

    c1->SaveAs("proton_polarization_plots.png");
    c1->Write();
    hLambdaLabFrame_pT->Write();
    hlambdaLabFrame_Y->Write();
    hProtonLambdaFrame_phi->Write();
    hlambdaLabFrame_phiDistr->Fit(fitPhiCollectiveFlowDistro);
    hlambdaLabFrame_phiDistr->Write();
    fProfV2pT->Write();
    fProfV1Y->Write();
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

Double_t Lambda_v1_v2_P_parameterization(TFile* Lambda_yield, Double_t fBVal, UParticle &ULambda){ //Valeriy's function for properly lambda generation

	Double_t centrality; //centrality for parameterization (b->centrality for Xe+Xe below)
	Double_t sNN = 2.87; // Energy of the collision in center-of-mass system
	if(fBVal<3.44)
                centrality=5;
        else if(fBVal<4.88)
                centrality=15;
        else if(fBVal<5.84)
                centrality=25;
        else if(fBVal<6.64)
                centrality=35;
        else if(fBVal<7.44)
                centrality=45;
        else if(fBVal<8.08)
                centrality=55;
        else if(fBVal<8.72)
                centrality=65;
        else if(fBVal<9.36)
                centrality=75;
        else if(fBVal<9.84)
                centrality=85;
        else
                centrality=95;


	// TFile* Lambda_yield = new TFile("/lhep/users/dflusova/lambda/afterburner/v.6/qa_out_xexe.root","read"); // File with TH2F pt-y
	TH2F* h_pt_y;
	if(centrality<10)
		h_pt_y = (TH2F*) Lambda_yield->Get("h2PartYpT010");
	else if(centrality>10 && centrality <40)
		h_pt_y = (TH2F*) Lambda_yield->Get("h2PartYpT");
	else
		h_pt_y = (TH2F*) Lambda_yield->Get("h2PartYpT40100");


        Double_t lambda_pT; // pT of Lambda from pT-y TH2F
	Double_t lambda_y;  // rapidity of Lambda from pT-y TH2F

	h_pt_y->GetRandom2(lambda_y,lambda_pT);

	Double_t v1 = (28.8635/(TMath::Power(sNN,2.89092))) * (  (-0.0233*centrality+0.5413* TMath::Power(centrality,1./3) ) * (163.536/18.0188 - 163.536/(lambda_pT +18.0188) )* lambda_y + ( -0.0056*centrality+0.377*TMath::Power(centrality,1./3) ) * (0.6653* lambda_pT - 0.6172 * TMath::Power(lambda_pT,2) +0.1154 * TMath::Power(lambda_pT,3) ) * TMath::Power(lambda_y,3) ); // v1 cent-pT-y func
	Double_t v2 = 0.07; //uniform distribution

    if( v1  > 1 ) v1 =1;
	// generate phi according to v1 and v2
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
      Double_t phi=f.GetRandom();


      //Global polarization parameterization

    //   Double_t mean_P = (2.8569/ TMath::Power(sNN,0.955513) ) * (2.4702 - 0.0461*centrality + 0.0042 * TMath::Power(centrality, 2)); // Value in % 

      Double_t fEnergyLambda = ULambda.GetMomentum().E();
      TLorentzVector vec;
      vec.SetPtEtaPhiE(lambda_pT, lambda_y, phi, fEnergyLambda);
      ULambda.SetMomentum(vec); //ULambda

      return (2.8569/ TMath::Power(sNN,0.955513) ) * (2.4702 - 0.0461*centrality + 0.0042 * TMath::Power(centrality, 2)); // Value in % 

}
