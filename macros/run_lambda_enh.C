#include <cstdlib>
#include <string>
#include "TSystem.h"
#include "add_enhanced_lambda.hpp"

void run_lambda_enh()
{
  auto getenv_s = [](const char* k, const char* def="") -> std::string {
    const char* v = gSystem->Getenv(k);
    return v ? std::string(v) : std::string(def);
  };

  std::string outFile = getenv_s("LAMBDA_OUTFILE", "out.root");
  long long fileStart = std::stoll(getenv_s("LAMBDA_FILESTART", "0"));
  long long fileCount = std::stoll(getenv_s("LAMBDA_FILESTART", "1"));
  std::string pattern = getenv_s("LAMBDA_PATTERN", "");
  std::string qaFile  = getenv_s("LAMBDA_QA_FILE", "qa_out_xexe.root");
  int enhanceStat     = std::stoi(getenv_s("LAMBDA_ENHANCE_STAT", "2"));
  unsigned seed       = (unsigned)std::stoul(getenv_s("LAMBDA_SEED", "12345"));

  add_enhanced_lambda(outFile.c_str(),
                   fileStart,
                   fileCount,
                   pattern.c_str(),
                   qaFile.c_str(),
                   enhanceStat,
                   seed);
}
