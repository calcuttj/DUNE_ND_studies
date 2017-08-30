#include "make_effs.c"

void Run_effs(int n_evt = -1,char * fInputFileName ="../DUNE_ND_ntps/full/fgt/ndtf_output_nu_fgt.dst.root" ,char * fOutFileName="try.root") { 
//int main(int n_evt = -1, char* tag = "FGT") { 
  gROOT->LoadMacro("make_effs.c");
  dune_dst t(0,fInputFileName);
  t.Loop(n_evt,fOutFileName);
//  return 0;
}

