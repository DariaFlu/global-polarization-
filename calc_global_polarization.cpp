
TString InFileName= "/lhep/users/dflusova/lambda/afterburner/v.6/out/processed/merged_result_urqmd_xexe_2.87gev_mf_6195240_.mcini.root";
TString OutFileName= "/lhep/users/dflusova/lambda/afterburner/v.6/out/processed/globalPol_merged_result_urqmd_xexe_2.87gev_mf_6195240.mcini.root";



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

        TVector3 polarizationAxis;
        TVector3 polarizationVec;

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

        // Double_t GetB(){ return fB;}

        void SetPolarizationSigma(Double_t _fSigmaPol) { fSigmaPol = _fSigmaPol; }

        void SetPolarizationAxis(){

            TVector3 b(fB, 0, 0); // Impact parameter direction
            TVector3 z(0, 0, 1);   // Beam axis
            polarizationAxis = z.Cross(b.Unit()); // y = z × b̂
            rotation.SetZAxis(polarizationAxis).RotateZ(0);
        }

        void SetPolarization(UParticle& lambda, Double_t polx, Double_t poly, Double_t polz)
        {
           if(polx || poly || polz) {
            TLorentzVector lambda_lab = lambda.GetMomentum();
            TVector3 beta = -lambda_lab.BoostVector();
            TLorentzVector lambda_rest = lambda_lab;
            lambda_rest.Boost(beta);

            lambda_rest.SetTheta(TMath::ACos(polz/TMath::Sqrt(polx*polx+poly*poly+polz*polz)));
            lambda_rest.SetPhi(TMath::Pi()+TMath::ATan2(-poly,-polx));
            lambda_rest.Boost(-beta);
            lambda.GetMomentum().SetTheta(lambda_rest.Theta());
            lambda.GetMomentum().SetPhi(lambda_rest.Phi());
            //std::cout << lamTV.Theta() << " = ";
                // lambda.GetMomentum().SetTheta(TMath::ACos(polz/TMath::Sqrt(polx*polx+poly*poly+polz*polz)));
                // lambda.GetMomentum().SetPhi(TMath::Pi()+TMath::ATan2(-poly,-polx));
            // lamTV.SetTheta(TMath::ACos(polz/TMath::Sqrt(polx*polx+poly*poly+polz*polz)));
            // lamTV.SetPhi(TMath::Pi()+TMath::ATan2(-poly,-polx));
            // //std::cout << lamTV.Theta() << std::endl;
            // lambda.SetMomentum(lamTV);
           } else {
              Double_t fPolarTheta = -99;
              Double_t fPolarPhi = -99;
           }
        }

        void RotateProtons() { //Add % polarization
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

        void RotateLambda(UParticle& lambda, Double_t fpoly, Double_t _fSigmaPol = 0.3){ 
            //std::cout << 1 << std::endl;
            SetPolarizationSigma(_fSigmaPol);
            
            Double_t polX = gRandom->Gaus(fpolx,fpolSx);  // generate polarization direction
            Double_t polY = gRandom->Gaus(fpoly,fpolSy);  // generate polarization direction
            Double_t polZ = gRandom->Gaus(fpolz,fpolSz);  // generate polarization direction
            // TVector3 polarizationVec(polX, polY, polZ);
            polarizationVec = TVector3(polX, polY, polZ);

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
            // part->SetPolarisation(polarizationVec.X(),polarizationVec.Y(),polarizationVec.Z());    
            SetPolarization(lambda, polarizationVec.X(),polarizationVec.Y(),polarizationVec.Z());
        }

        void Clear(){
            UEvent::Clear();
            fLambdas->Clear();
            fProtons->Clear();
            fPions->Clear();
        }

        ClassDef(UUEvent, 3); // Version number
    };



double generateCosTheta(double alpha, double PolarizationPercent = 0.6) {

    TF1* randomCosTheta = new TF1("randomCosTheta", "1+[0]*[1]*x", -1.0, 1.0);
    // Set the parameters
    randomCosTheta->SetParameters(alpha, PolarizationPercent);
    // Generate a random value from the distribution
    double value = randomCosTheta->GetRandom();
    
    // Clean up
    delete randomCosTheta;
    
    return value;

}

Double_t positivePhi(const Double_t& phi){
    Double_t phi_rot = phi;
    if (phi < 0) phi_rot += 2.*TMath::Pi();
    return phi_rot;
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

void calc_global_polarization(TString InFileName, TString OutFileName){
    gStyle->SetOptStat(0);  // Disable stats globally
    TF1* fitPhiDistro = new TF1("fitPhiDistro", "[0]*(1+2*[1]*TMath::Sin(x)+2*[2]*TMath::Cos(x))", 0, 2*TMath::Pi()); //Fitting function
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
    TProfile* fProfV2pT = new TProfile("fProfV2pT", "v_{2} vs pT", 20, 0, 2, 0, 10);
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
    Double_t fBMax = 7.44;
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

            TLorentzVector lambda_rest = lambda_lab;

            lambda_rest.Boost(beta);
            // hlambdaLabFrame_phiDistr->Fill((lambda_rest.Phi()));


            // TLorentzVector proton_lab = proton->GetMomentum(); // getting laboratory proton for further properly binning over [pT;Y]
            // proton_lab.Boost(beta);                             //boost proton to lab frame

            // //if bin isn't uderflow or overflow then put phi*-distribution in histogram[number of bin]
            // if(returnNumberOfBin(proton_lab.Pt(),       fpTMin, fpTMax, NpTBins) > -1) hProtonLambdaFrame_phi_pTBins[(returnNumberOfBin(proton_lab.Pt(),       fpTMin, fpTMax, NpTBins))]->Fill(phiStar); 
            // if(returnNumberOfBin(proton_lab.Rapidity(), fYMin,  fYMax,  NYBins ) > -1) hProtonLambdaFrame_phi_YBins [(returnNumberOfBin(proton_lab.Rapidity(), fYMin,  fYMax,  NYBins ))]->Fill(phiStar); 

            // //Filling pT and Y histograms
            // hProtonLabFrame_pT->Fill(proton_lab.Pt());
            // hlambdaLabFrame_Y->Fill(proton_lab.Rapidity());

            //if bin isn't uderflow or overflow then put phi*-distribution in histogram[number of bin]
            if(returnNumberOfBin(lambda_lab.Pt(),       fpTMin, fpTMax, NpTBins) > -1) hProtonLambdaFrame_phi_pTBins[(returnNumberOfBin(lambda_lab.Pt(),       fpTMin, fpTMax, NpTBins))]->Fill(phiStar); 
            if(returnNumberOfBin(lambda_lab.Rapidity(), fYMin,  fYMax,  NYBins ) > -1) hProtonLambdaFrame_phi_YBins [(returnNumberOfBin(lambda_lab.Rapidity(), fYMin,  fYMax,  NYBins ))]->Fill(phiStar); 

            // if(returnNumberOfBin(proton_lab.Pt(),       fpTMin, fpTMax, NpTBins) > -1) hProtonLambdaFrame_phi_pTBins[(returnNumberOfBin(proton_lab.Pt(),       fpTMin, fpTMax, NpTBins))]->Fill(phiStar); 
            // if(returnNumberOfBin(proton_lab.Rapidity(), fYMin,  fYMax,  NYBins ) > -1) hProtonLambdaFrame_phi_YBins [(returnNumberOfBin(proton_lab.Rapidity(), fYMin,  fYMax,  NYBins ))]->Fill(phiStar); 


            //Filling pT and Y histograms
            hLambdaLabFrame_pT->Fill(lambda_lab.Pt());
            hlambdaLabFrame_Y ->Fill(lambda_lab.Rapidity());
            hlambdaLabFrame_phiDistr->Fill(positivePhi(lambda->GetMomentum().Phi()));
            // if(lambda->GetMomentum().Rapidity() > -0.75 && lambda->GetMomentum().Rapidity() > -0.75 ) hlambdaLabFrame_phiDistrBin[(returnNumberOfBin(lambda_lab.Pt(), fpTMin, fpTMax, NpTBins))]->Fill(phiStar); 
            // if (inUEvent->GetB() > fBMin && inUEvent->GetB() < fBMax && lambda->GetMomentum().Rapidity() > -0.75 && lambda->GetMomentum().Rapidity() > -0.75 ){
            if ( lambda->GetMomentum().Rapidity() > -0.75 && lambda->GetMomentum().Rapidity() > -0.75 ){
                std::cout<<"fB = "<<inUEvent->GetB()<<" Y = "<< lambda->GetMomentum().Rapidity()<<std::endl;
                fProfV2pT->Fill(lambda_lab.Pt(), TMath::Cos(2*positivePhi(lambda->GetMomentum().Phi())));
            }
            // 
            // hlambdaLabFrame_phiDistr->Fill(positivePhi(lambda_rest.Phi()));
            // hlambdaLabFrame_phiDistr->Fill((lambda_lab.Phi()));
            // hlambdaLabFrame_phiDistr->Fill(positivePhi(lambda_rest.Theta()));
            // hlambdaLabFrame_phiDistr->Fill(positivePhi(lambda->GetMomentum().Theta()));


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
        hphi_YBins->SetBinError(iHisto+1, (fitPhiDistro->GetParError(1)*8)/(TMath::Pi()*anisotropy)*100);
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
    outFile->Close();
    inFile->Close();
}