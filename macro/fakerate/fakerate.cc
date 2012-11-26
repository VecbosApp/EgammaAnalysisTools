#include "estimateFakeRate.C"

#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace std;

int main(int argc, char* argv[]) {

  char inname[500], outname[500];

  if( argc < 3) {
    cout << "missing argument: fakerate <inputfilebase> <basename>" << endl;
    return 1;
  }
  strcpy(inname,argv[1]);
  strcpy(outname,argv[2]);

  stringstream infilename, isofriend, idfriend;
  infilename << "../results_data/" << string(inname) << ".root";
  isofriend << "../results_data/" << string(inname) << "_hzzisoFriend.root";
  idfriend << "../results_data/" << string(inname) << "_hzzidbitsFriend.root";

  cout << "run the FR for triggering electrons on " << infilename.str() << "..." << endl;
  cout << "   iso friend file = " << isofriend.str() << "..." << endl;
  cout << "   id friend file = " << idfriend.str() << "..." << endl;
  TFile *file = TFile::Open(infilename.str().c_str());
  TTree *tree = (TTree*)file->Get("eleIDdir/T1");
  estimateFakeRate analyzer(tree);
  analyzer.addIsoFriend(isofriend.str().c_str());
  analyzer.addIdBitsFriend(idfriend.str().c_str());
  TString outfileBias(outname);
  outfileBias += TString("_trigger");
  analyzer.Loop(outfileBias);
  cout << "DONE triggering electrons." << endl;
  file->Close();

  // cout << "run the FR for non triggering electrons..." << endl;
  // TFile *file2 = TFile::Open("../results_data/fakes-zeeOneFake.root");
  // TTree *tree2 = (TTree*)file2->Get("eleIDdir/T1");
  // estimateFakeRate analyzer2(tree2);
  // analyzer2.addIsoFriend("../results_data/fakes-zeeOneFake_hzzisoFriend.root");
  // analyzer2.addIdBitsFriend("../results_data/fakes-zeeOneFake_hzzidbitsFriend.root");
  // TString outfileUnbias(outname);
  // outfileUnbias += TString("_zee1fake");
  // analyzer2.Loop(outfileUnbias);
  // cout << "DONE unbiased electrons." << endl;

  return 0;
}
