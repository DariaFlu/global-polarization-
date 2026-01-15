// macros/run_sim.C
#include <string>
#include "simulate_lambda_decays.hpp"  

void run_sim()
{
  simulate_lambda_decays(
    "/eos/nica/mpd/users/parfenov/SimData/UrQMD/xexe_2.87gev_mf/6195240/files/mcini/urqmd_xexe_2.87gev_mf_6195240_5.mcini.root", 
    "/lhep/users/dflusova/lambda/dummy/global-polarization-/result_test_urqmd_xexe_2.87gev_mf_6195240_5.mcini.root", 
    "/lhep/users/dflusova/lambda/afterburner/MPD_afterburner/qa_out_xexe.root", 
    2);
}
