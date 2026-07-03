// macros/run_phi_flow_post.C
#include "post_phi_flow.h"   // 
#include <string>

void run_phi_flow_post(const char* outFile,
                       long long fileStart,
                       long long fileCount,
                       const char* pattern,
                       const char* inDirOrQA,   // если нужно
                       int onlyEnhanced,
                       unsigned seed)
{
  post_phi_flow(outFile, fileStart, fileCount, pattern, inDirOrQA, onlyEnhanced, seed);
}
