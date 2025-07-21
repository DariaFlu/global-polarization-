#include "UParticle.h"
#include "UEvent.h"

#include "read_unigen_root.cpp"


#include "MpdEvent.h"  
#include "MpdTrack.h"
#include "MpdMCTrack.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include <iostream>


void lambdaPolAnal() {

    TString pathIn = "/eos/nica/mpd/users/parfenov/SimData/UrQMD/xexe_2.87gev_mf/6195240/files/mcini/"; //directory with input files
    TString pathConf = "/lhep/users/dflusova/lambda/afterburner/v.6/"; //out folder
    TString pathOut = "/lhep/users/dflusova/lambda/afterburner/v.6/out/"; //out folder

    TString fileIn = "urqmd_xexe_2.87gev_mf_6195240_713.mcini.root"; //name of input file; lambda from UrQMD generator
    TString fileConfIn = "qa_out_xexe.root"; //file with pT, Y, centrality dependencies
    TString fileOut = "result_urqmd_xexe_2.87gev_mf_6195240.mcini.root"; //file for storing produced protons 
    TString OutFileName = "result_global_polarization_urqmd_xexe_2.87gev_mf_6195240.mcini.root"; //resulting .root file, it contains deltaPhi distributions

    Int_t nFile = 5; //number of input file
    Int_t enhancedValue = 6; //number of enhanced lambdass

    for(Int_t iFile = 1; iFile < nFile; ++iFile){// Uncomment the loop if it is needed for proton production
        //Notice! simulate_lambda_decays() produced protons and pions in lambda state frame! 
        //In the output file lambda in laboratory frame   
        simulate_lambda_decays(TString::Format("%surqmd_xexe_2.87gev_mf_6195240_%i.mcini.root", pathIn.Data(), iFile), 
                                TString::Format("%s%s",pathOut.Data(),fileOut.Data()),
                                TString::Format("%s%s",pathConf.Data(),fileConfIn.Data()), 
                                iFile, enhancedValue);
    }
    //Polarization proccesing 
    calc_global_polarization(TString::Format("%s%s",pathOut.Data(), fileOut.Data()),  TString::Format("%s%s",pathOut.Data(), OutFileName.Data()));

    std::cout << "ROOT file saved and closed." << std::endl;
    std::cout << "Analysis complete! Results saved to lambda_polarization_results.root" << std::endl;

}



