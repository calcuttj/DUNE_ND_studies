#include "make_effs.c"

void Run_effs(int n_evt = -1, char* tag = "FGT",char * fOutFileName="try.root") { 
//int main(int n_evt = -1, char* tag = "FGT") { 
  gROOT->LoadMacro("make_effs.c");
  dune_dst t;
  t.Loop(n_evt,tag,fOutFileName);
//  return 0;
}

