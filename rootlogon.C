void rootlogon(){
    gROOT->ProcessLine("gSystem->Load(\"/lhep/users/dflusova/lambda/afterburner/v.7/global-polarization-/build/libgp_dict\")");
    gROOT->ProcessLine("gSystem->Load(\"/lhep/users/dflusova/lambda/afterburner/v.7/global-polarization-/build/libgp_macros\")");
    gROOT->ProcessLine(".include \"/lhep/users/dflusova/lambda/afterburner/v.7/global-polarization-/build/../include\"");
    gROOT->ProcessLine("#include \"/lhep/users/dflusova/lambda/afterburner/v.7/global-polarization-/build/../include/afterburner.hpp\"");
}
