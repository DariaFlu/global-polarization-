#include "calc_pol_vs_Nenh.hpp"


void calc_pol_vs_Nenh(TString InFileName, TString OutFileName, std::vector<Int_t> &enhancedFlag){
    // TFile *inFile = TFile::Open(InFileName);
    gStyle->SetOptStat(0);
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

    inChain->SetBranchAddress("LambdaPolarization", &vecPol);
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

    Double_t fpTMin = 0.2;//GeV/C
    Double_t fpTMax = 1.6;//GeV/C
    Double_t fYmin = -0.75;
    Double_t fYmax = 0.75;

    TFile *outFile = TFile::Open(OutFileName, "RECREATE");//Resulting out file


    TF1* fitPhiDistro = new TF1("fitPhiDistro", "[0]*(1+2*[1]*TMath::Sin(x)+2*[2]*TMath::Cos(x))", 0, 2*TMath::Pi()); //Fitting function
    fitPhiDistro->SetParameter(1, 4.*(TMath::Pi()*0.732)/(8.*100)); //initial guess for fit
    const double anisotropy = 0.732;

    TH1D* hPolarisationBins = new TH1D("hPolarisationBins", "hPolarisationBins", 50, 0, 50);
    TH1D* hPolarisationBinsErrors = new TH1D("hPolarisationBinsErrors", "hPolarisationBinsErrors", 50, 0, 50);

    Double_t fCenMin = 10;//%
    Double_t fCenMax = 40;//%

    for(size_t iEnh = 0; iEnh < enhancedFlag.size(); iEnh++){
        Int_t enhFlag = enhancedFlag[iEnh];
        std::cout << "eng flag: " << enhFlag << std::endl;
        TH1D* hProtonLambdaFrame_phi = new TH1D("hProtonLambdaFrame_phi", " ; #Delta#phi, rad;Counts", 50,  0, 2*TMath::Pi()); //proton phi* distribution (lambda frame)

        for(Long64_t iEvent = 0; iEvent < nEvents; iEvent++){
            // inTree->GetEntry(iEvent);
            inChain->GetEntry(iEvent);
            int lambdas = 0;

            for (size_t pols_i = 0; pols_i < vecPol->size(); pols_i++){
                if(get_centrality(inEvent->GetB()) < fCenMin || get_centrality(inEvent->GetB()) > fCenMax) continue;

                hPx_Nlamb->Fill(enhFlag, vecPol->at(pols_i).X());
                hPy_Nlamb->Fill(enhFlag, vecPol->at(pols_i).Y());
                hPmag_lamb->Fill(enhFlag, TMath::Sqrt(vecPol->at(pols_i).X()*vecPol->at(pols_i).X()+vecPol->at(pols_i).Y()*vecPol->at(pols_i).Y()));
                Double_t phiStar = get_positive_phi(UProton->at(pols_i).GetMomentum().Phi()); //Get phi*
                //  hProtonLambdaFrame_phi->Fill(phiStar); //phi* distribution 
                if((ULambda->at(pols_i).GetMomentum().Pt() < 0.34 && ULambda->at(pols_i).GetMomentum().Pt() > 0.2) && (ULambda->at(pols_i).GetMomentum().Rapidity() < fYmax && ULambda->at(pols_i).GetMomentum().Rapidity() > fYmin)) hProtonLambdaFrame_phi->Fill(phiStar); //phi* distribution 
                if (ULambda->at(pols_i).GetMate() == -9){
                    lambdas++;


                    hPx_Nenh->Fill(enhFlag, vecPol->at(pols_i).X());
                    hPy_Nenh->Fill(enhFlag, vecPol->at(pols_i).Y());
                    hPmag_Nenh->Fill(enhFlag, TMath::Sqrt(vecPol->at(pols_i).X()*vecPol->at(pols_i).X()+vecPol->at(pols_i).Y()*vecPol->at(pols_i).Y()));
                    if (lambdas > enhFlag) break;
                }

            }
        }
        hProtonLambdaFrame_phi->Fit(fitPhiDistro);//Fitting phi* full distribution  (without pT or Y binning)
        Double_t polRaw = (fitPhiDistro->GetParameter(1)*8.)/(TMath::Pi()*anisotropy)*100.;
        Double_t polRawErr = (fitPhiDistro->GetParError(1)*8.)/(TMath::Pi()*anisotropy)*100.;

        hProtonLambdaFrame_phi->Write();

        std::cout << "pol percent: " << polRaw << " pol error: " << polRawErr << std::endl;

        hPolarisationBins->SetBinContent((enhancedFlag.at(iEnh) + 1), polRaw);
        hPolarisationBins->SetBinError((enhancedFlag.at(iEnh) + 1), polRawErr);

        hPolarisationBinsErrors->SetBinContent((enhancedFlag.at(iEnh) + 1), (polRawErr/polRaw)*100);

        delete hProtonLambdaFrame_phi;
    }
    gStyle->SetErrorX(0);

    TCanvas* cPolarisationBins = new TCanvas("cPolarisationBins","cPolarisationBins", 1200, 800);
    cPolarisationBins->cd();
    // Draw histogram with markers and without bin width lines
    hPolarisationBins->Draw("E1 P");  // "P" option for points/markers
    // hPolarisationBins->SetLineWidth(0);  // Switch off bin lines
    // Set axis labels
    hPolarisationBins->GetXaxis()->SetTitle("enhancement");
    hPolarisationBins->GetYaxis()->SetTitle("P,  %");

    // Optional: Customize marker style if needed
    hPolarisationBins->SetMarkerStyle(24);  // Circle markers
    hPolarisationBins->SetMarkerSize(1);    // Marker size
    cPolarisationBins->SaveAs("cPolarisationBins.C");

    TCanvas* cPolarisationBinsErrors = new TCanvas("cPolarisationBinsErrors","cPolarisationBinsErrors", 1200, 800);
    cPolarisationBinsErrors->cd();
    // Set axis labels
    hPolarisationBinsErrors->GetXaxis()->SetTitle("enhancement");
    hPolarisationBinsErrors->GetYaxis()->SetTitle("error P,  %");
    hPolarisationBinsErrors->Draw("P");
    cPolarisationBinsErrors->SaveAs("cPolarisationBinsErrors.C");
    hPolarisationBinsErrors->Write();

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
    hPolarisationBins->Write();

    outFile->Close();
}
