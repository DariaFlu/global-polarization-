#include <cstdlib>
#include <string>
#include "TSystem.h"
#include "add_enhanced_phi.hpp"

void run_phi_enh()
{
  auto getenv_s = [](const char* k, const char* def="") -> std::string {
    const char* v = gSystem->Getenv(k);
    return v ? std::string(v) : std::string(def);
  };

  std::string outFile = getenv_s("PHI_OUTFILE", "out.root");
  long long fileStart = std::stoll(getenv_s("PHI_FILESTART", "0"));
  long long fileCount = std::stoll(getenv_s("PHI_FILECOUNT", "1"));
  std::string pattern = getenv_s("PHI_PATTERN", "");
  std::string qaFile  = getenv_s("PHI_QA_FILE", "qa.root");
  int enhanceStat     = std::stoi(getenv_s("PHI_ENHANCE_STAT", "2"));
  unsigned seed       = (unsigned)std::stoul(getenv_s("PHI_SEED", "12345"));

  add_enhanced_phi(outFile.c_str(),
                   fileStart,
                   fileCount,
                   pattern.c_str(),
                   qaFile.c_str(),
                   enhanceStat,
                   seed);
}
