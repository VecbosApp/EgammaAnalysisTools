#include "estimateFakeRate.C"

#include <iostream>
#include "TFile.h"
#include "TTree.h"

using namespace std;

int main(int argc, char* argv[]) {

  char outname[300];

  if( argc < 2) {
    cout << "missing argument: fakerate <basename>" << endl;
    return 1;
  }
  strcpy(outname,argv[1]);

  TFile *file = TFile::Open("../results_data_fakes/merged.root");
  TTree *tree = (TTree*)file->Get("eleIDdir/T1");

  estimateFakeRate analyzer(tree);
  analyzer.Loop(outname);

  return 0;
}
