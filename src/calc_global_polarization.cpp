#include "calc_global_polarization.hpp"

void calc_global_polarization(TString fileIn , Int_t NFiles, TString OutFileName, Int_t enhancedFlag){
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

    // Double_t fpTMax_bin = 1.6; // max pT [GeV/c]
    // Double_t fpTMin_bin = 0.2; // min pT [GeV/c]
    // Double_t fYMax_bin =  1.5; // max rapidity
    // Double_t fYMin_bin = -1.5; // min rapidity

    Double_t fpTMax_bin = 1.6; // max pT [GeV/c]
    Double_t fpTMin_bin = 0.2; // min pT [GeV/c]
    Double_t fYMax_bin =  0.75; // max rapidity
    Double_t fYMin_bin = -0.75; // min rapidity

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

    TChain *inChain = new TChain("decays");
    //inChain->Add(InFileName);
    // TString pathIn = "/scratch3/dflusova/afterburner/out/"; //directory with input files
    // TString fileIn = "result_urqmd_xexe_2.87gev_mf_6195240"; //file for storing produced protons 
    for (int i = 1; i < NFiles; i++){
        // if (!(TFile::Open(TString::Format("%s%s_%i.mcini.root",pathIn.Data(), fileIn.Data(), i), "READ")->IsZombie())){
        TString fname = TString::Format("%s_%i.mcini.root", fileIn.Data(), i);
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

    inChain->SetBranchAddress("LambdaPolarization", &vecPol);
    inChain->SetBranchAddress("UProton", &UProton);
    inChain->SetBranchAddress("ULambda", &ULambda);
    inChain->SetBranchAddress("UPion", &UPion);

    // Double_t fBMin = 3.44;
    // Double_t fBMax = 6.64;
    Double_t fCenMin = 10;//%
    Double_t fCenMax = 40;//%
    Double_t fpTMin = 0.2;//GeV/C
    Double_t fpTMax = 1.6;//GeV/C
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

        Int_t lambdaCounter = 0;
        Double_t Qx = 0, Qy = 0;
        Double_t QxA = 0, QyA = 0, QxB = 0, QyB = 0;
        Double_t fPsi = 0;
        for(size_t iTrack = 0; iTrack < vecPol->size(); iTrack++){
            
            if(get_centrality(inEvent->GetB()) < fCenMin || get_centrality(inEvent->GetB()) > fCenMax) continue;
            if(ULambda->at(iTrack).GetMate() == -9 || ULambda->at(iTrack).GetMate() == 15) continue;
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
            if (ULambda->at(iTrack).GetMomentum().Pt()       > fpTMin && ULambda->at(iTrack).GetMomentum().Pt()       < fpTMax) fProfV1Y ->Fill(ULambda->at(iTrack).GetMomentum().Rapidity(),TMath::Cos(get_positive_phi(ULambda->at(iTrack).GetMomentum().Phi())));
           
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
    TLegend* leg_pT[NpTBins];
    TCanvas* canv_pT[NpTBins];

    for(Int_t iHisto = 0; iHisto < NpTBins; iHisto++) { //Loop over pT bins
        hProtonLambdaFrame_phi_pTBins[iHisto]->Fit(fitPhiDistro ,"L");
        hProtonLambdaFrame_phi_pTBins[iHisto]->Write();
        hphi_pTBins->SetBinContent(iHisto+1, (fitPhiDistro->GetParameter(1)*8)/(TMath::Pi()*anisotropy)*100);
        hphi_pTBins->SetBinError  (iHisto+1, (fitPhiDistro->GetParError(1)*8)/(TMath::Pi()*anisotropy)*100);

        leg_pT[iHisto] = new TLegend(0.50, 0.58, 0.9, 0.89);
        leg_pT[iHisto]->SetHeader("p_{0}(1 + 2p_{1}sin(#phi) + 2p_{2}cos(#phi))", "C");
        leg_pT[iHisto]->SetTextSize(0.05); 
        leg_pT[iHisto]->SetBorderSize(0);   
        leg_pT[iHisto]->AddEntry("", Form("p_{1} = %.3f #pm %.3f", fitPhiDistro->GetParameter(1), fitPhiDistro->GetParError(1)), "");
        leg_pT[iHisto]->AddEntry("", Form("P_{#Lambda} = %.2f #pm %.2f %%", (fitPhiDistro->GetParameter(1)*8)/(TMath::Pi()*anisotropy)*100, (fitPhiDistro->GetParError(1)*8)/(TMath::Pi()*anisotropy)*100), "");
        leg_pT[iHisto]->AddEntry("", Form("P_{#Lambda} = (8p_{1}) / (%.3f#pi)", anisotropy), "");
        canv_pT[iHisto] = new TCanvas(TString::Format("hProtonLambdaFrame_phi_pTBin__%i", iHisto), TString::Format("hProtonLambdaFrame_phi_pTBin__%i", iHisto), 1200, 800);
        canv_pT[iHisto]->cd();
        hProtonLambdaFrame_phi_pTBins[iHisto]->Draw();
        leg_pT[iHisto]->Draw("same");
        canv_pT[iHisto]->Write();
         
    }

    hphi_pTBins->SetMarkerStyle(23); // Filled circle (ROOT default)
    hphi_pTBins->SetMarkerSize(1.2); // Slightly larger markers
    hphi_pTBins->Write();

    //--Y binning--//
    //same logic like in previous one
    TH1D* hphi_YBins= new TH1D("hphi_YBins", " rapidity binning; Y; P_{#lambda} %", NYBins,  fYMin_bin, fYMax_bin);
    TLegend* leg_Y[NYBins];
    TCanvas* canv_Y[NYBins];
    for(Int_t iHisto = 0; iHisto < NYBins; iHisto++) {
        hProtonLambdaFrame_phi_YBins[iHisto]->Fit(fitPhiDistro, "L");
        hProtonLambdaFrame_phi_YBins[iHisto]->Write();
        hphi_YBins->SetBinContent(iHisto+1, (fitPhiDistro->GetParameter(1)*8)/(TMath::Pi()*anisotropy)*100);
        hphi_YBins->SetBinError  (iHisto+1, (fitPhiDistro->GetParError(1)*8)/ (TMath::Pi()*anisotropy)*100);

        leg_Y[iHisto] = new TLegend(0.50, 0.58, 0.9, 0.89);
        leg_Y[iHisto]->SetHeader("p_{0}(1 + 2p_{1}sin(#phi) + 2p_{2}cos(#phi))", "C");
        leg_Y[iHisto]->SetTextSize(0.05); 
        leg_Y[iHisto]->SetBorderSize(0);   
        leg_Y[iHisto]->AddEntry("", Form("p_{1} = %.3f #pm %.3f", fitPhiDistro->GetParameter(1), fitPhiDistro->GetParError(1)), "");
        leg_Y[iHisto]->AddEntry("", Form("P_{#Lambda} = %.2f #pm %.2f %%", (fitPhiDistro->GetParameter(1)*8)/(TMath::Pi()*anisotropy)*100, (fitPhiDistro->GetParError(1)*8)/(TMath::Pi()*anisotropy)*100), "");
        leg_Y[iHisto]->AddEntry("", Form("P_{#Lambda} = (8p_{1}) / (%.3f#pi)", anisotropy), "");
        canv_Y[iHisto] = new TCanvas(TString::Format("hProtonLambdaFrame_phi_YBin__%i", iHisto), TString::Format("hProtonLambdaFrame_phi_YBin__%i", iHisto), 1200, 800);
        canv_Y[iHisto]->cd();
        hProtonLambdaFrame_phi_YBins[iHisto]->Draw();
        leg_Y[iHisto]->Draw("same");
        canv_Y[iHisto]->Write();

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
    cProton->SaveAs((enhancedFlag > -1 ) ? TString::Format("picture/proton_polarization_plots_%i.C",enhancedFlag) : "picture/proton_polarization_plots.C");

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