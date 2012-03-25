#include "estimateFakeRate.C"

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

  cout << "run the FR for triggering electrons..." << endl;
  TFile *file = TFile::Open("../results_data/fakes.root");
  TTree *tree = (TTree*)file->Get("eleIDdir/T1");
  estimateFakeRate analyzer(tree);
  analyzer.addIsoFriend("../results_data/fakes_hzzisoFriend.root");
  TString outfileBias(outname);
  outfileBias += TString("_trigger");
  analyzer.Loop(outfileBias);
  cout << "DONE triggering electrons." << endl;
  file->Close();

  cout << "run the FR for non triggering electrons..." << endl;
  TFile *file2 = TFile::Open("../results_data/fakes-zeeOneFake.root");
  TTree *tree2 = (TTree*)file2->Get("eleIDdir/T1");
  estimateFakeRate analyzer2(tree2);
  analyzer2.addIsoFriend("../results_data/fakes-zeeOneFake_hzzisoFriend.root");
  TString outfileUnbias(outname);
  outfileUnbias += TString("_zee1fake");
  analyzer2.Loop(outfileUnbias);
  cout << "DONE unbiased electrons." << endl;

  return 0;
}
