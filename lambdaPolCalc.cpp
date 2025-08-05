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


void lambdaPolCalc() {
    gSystem->Load("AutoDict_vector_TVector3__cxx.so");
    gSystem->Load("AutoDict_vector_UParticle__cxx.so");

    // TString pathIn = "/scratch3/dflusova/afterburner/out/processed/"; //directory with input files
    // TString pathIn = "/lhep/users/dflusova/lambda/afterburner/v.6/out/"; //directory with input files
    TString pathIn = "/scratch3/dflusova/afterburner/out/"; //directory with input files

    // TString pathOut = "/lhep/users/dflusova/lambda/afterburner/v.6/calcAnal/"; //out folder
    TString pathOut = "/lhep/users/dflusova/lambda/afterburner/v.6/out/"; //out folder

    // TString fileIn = "merged_result_urqmd_xexe_2.87gev_mf_6195240_.mcini.root"; //file for storing produced protons 
    // TString fileIn = "result_urqmd_xexe_2.87gev_mf_6195240.mcini.root"; //file for storing produced protons 
    TString fileIn = "result_urqmd_xexe_2.87gev_mf_6195240"; //file for storing produced protons 

    TString fileOutName = "result_global_polarization_urqmd_xexe_2.87gev_mf_6195240"; //resulting .root file, it contains deltaPhi distributions
    TString fileOutFormat = ".mcini.root"; //resulting .root file, it contains deltaPhi distributions
    // std::vector<Int_t> vecEnhancedValue = {0, 1, 2, 5, 10, 50};
    std::vector<Int_t> vecEnhancedValue = {0, 50};
    // TFile* file;
    // for (auto iMult : vecEnhancedValue) calc_global_polarization(TString::Format("%s%s",pathIn.Data(), fileIn.Data()),  TString::Format("%s%s_%i_%s",pathOut.Data(), fileOutName.Data(), iMult, fileOutFormat.Data()), iMult);
    // for (int i = 1; i < 2000; i++){
    //     file = TFile::Open(TString::Format("%s%s_%i.mcini.root",pathIn.Data(), fileIn.Data(), i), "READ");
    //     if (file && !file->IsZombie()) {
    //     // if (!(TFile::Open(TString::Format("%s%s_%i.mcini.root",pathIn.Data(), fileIn.Data(), i), "READ")->IsZombie())){
    //     // if (gSystem->AccessPathName(TString::Format("%s%s_%i.mcini.root",pathIn.Data(), fileIn.Data(), i))){
    //         for (auto iMult : vecEnhancedValue) calc_global_polarization(TString::Format("%s%s_%i.mcini.root",pathIn.Data(), fileIn.Data(), i),  TString::Format("%s%s_iMult%i%s",pathOut.Data(), fileOutName.Data(), iMult, fileOutFormat.Data()), iMult);
    //     // }
    //     }
    //     else {
    //         // Handle the case when file couldn't be opened
    //         std::cerr << "Error: Could not open file " 
    //                   << TString::Format("%s%s_%i.mcini.root", pathIn.Data(), fileIn.Data(), i)
    //                   << std::endl;
    //     }
    // }   
    calc_global_polarization(TString::Format("%s%s_%i.mcini.root",pathIn.Data(), fileIn.Data(), 1),  TString::Format("%s%s%s",pathOut.Data(), fileOutName.Data(), fileOutFormat.Data()), 50);
    std::cout << "ROOT file saved and closed." << std::endl;
    std::cout << "Analysis complete! Results saved to lambda_polarization_results.root" << std::endl;

}





