// ====================================================================
// usage. In a ROOT session, do:
// .L src/createPdfsTreeFromRooDataset.cxx
// creatRootTree("sPlots/zTaP.root")
// where the file sPlots/zTaP.root is the one containing the RooDataSet
// ==================================================================== 

#include <TFile.h>
#include <RooDataSet.h>
#include <TTree.h>

void creatRootTree(const char *sPlotDatasetRoot) {

  TFile *fileIn = TFile::Open(sPlotDatasetRoot);
  RooDataSet *dataset = (RooDataSet*)fileIn->Get("dataset");
  fileIn->Close();

  TFile *fileOut = TFile::Open("sPlots/sPlotsZee_tree.root","recreate");
  dataset->tree().Write();
  fileOut->Close();

}
