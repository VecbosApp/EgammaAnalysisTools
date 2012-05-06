#include "estimateMuonFakeRateHzz4lTree.C"

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace std;

int main(int argc, char* argv[]) {

  char outname[300];

  if( argc < 2) {
    cout << "missing argument: fakerate <basename>" << endl;
    return 1;
  }
  strcpy(outname,argv[1]);

  cout << "run the FR for non triggering muons..." << endl;
  TFile *file = TFile::Open("hzzTree.root");
  TTree *tree = (TTree*)file->Get("zllmtree/probe_tree");
  estimateMuonFakeRateHzz4lTree analyzer(tree);
  TString outfileUnbias(outname);
  outfileUnbias += TString("_zllm");
  analyzer.Loop(outfileUnbias);
  cout << "DONE unbiased mu." << endl;

  return 0;
}
