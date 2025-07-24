#include <iostream>
#include "TChain.h"
#include "TSystem.h"

void merger(TString nameInputFile = "result_urqmd_xexe_2.87gev_mf_6195240_") {
    // Define the base path and the tree name
    TString basePath = "/lhep/users/dflusova/lambda/afterburner/v.6/out/";
    TString treeName = "decays";
    TString outputPath = "/lhep/users/dflusova/lambda/afterburner/v.6/out/processed/";

    //Create directories if they don't exist
    if (gSystem->mkdir(outputPath, true) == -1) {
        std::cerr << "Error: Could not create output directory " << outputPath << std::endl;
    }

    // Create a TChain
    TChain chain(treeName);

    // Add all files for the specified nameInputFile using wildcards
    TString filePattern = basePath + nameInputFile + "*.mcini.root";
    std::cout << "Looking for files matching pattern: " << filePattern << std::endl;

    // Add files to the chain
    int numFiles = chain.Add(filePattern);
    if (numFiles == 0) {
        std::cerr << "Error: No files found matching the pattern " << filePattern << std::endl;
        return;
    }

    // Print the list of files being added
    TObjArray* fileList = chain.GetListOfFiles();
    for (int i = 0; i < fileList->GetEntries(); i++) {
        std::cout << "Adding file: " << fileList->At(i)->GetTitle() << std::endl;
    }

    // Merge the files into a single output file
    TString outputFile = outputPath + "merged_" + nameInputFile + ".mcini.root";
    std::cout << "Merging " << numFiles << " files into " << outputFile << std::endl;

    chain.Merge(outputFile);

    std::cout << "Merged files for nameInputFile " << nameInputFile << " into " << outputFile << std::endl;
}