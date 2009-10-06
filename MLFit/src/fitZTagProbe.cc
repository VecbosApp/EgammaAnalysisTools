// Set Fit Options
MLOptions GetDefaultOptions() {
  MLOptions opts;
  // Fit configuration
  opts.addBoolOption("useinvMass", "Use invMass", kTRUE);        
  return opts;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

MLFit theFit;

void myFit() {

  MLFit theFit;

  // Various fit options...
  MLOptions opts = GetDefaultOptions();
  
  // define the structure of the dataset
  RooRealVar *zmass = new RooRealVar("zmass","e^{+}e^{-} Mass [GeV/c^{2}]",90,40,110);
  theFit.AddFlatFileColumn(zmass);
  
  // define a fit model
  theFit.addModel("myFit", "Tag And Probe Zee");
  
  // define species
  theFit.addSpecies("myFit", "sig", "Signal Component");
  theFit.addSpecies("myFit", "bkg", "Bkg   Component");
  
  // mLL PDF
  if(opts.getBoolVal("useinvMass")) {
    theFit.addPdfWName("myFit", "sig" , "zmass",  "Cruijff", "sig_Mass");
    theFit.addPdfWName("myFit", "bkg" , "zmass",  "Poly2",  "bkg_Mass");
  }
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Fit a sample of Z events
void FitZeeTagAndProbe() {
  
  myFit();

  // Load the data
  //  theFit.addDataSetFromRootFile("T1", "T1", "/cmsrm/pc21/emanuele/Likelihood2.1.X/PdfsDatasets/ZJetsMADGRAPHPdfsDataset.root");
  //  theFit.addDataSetFromRootFile("T1", "T1", "/cmsrm/pc21/emanuele/Likelihood2.1.X/PdfsDatasets/WpTTJetsMADGRAPHPdfsDataset.root");
  // W,Z+jets: 300 pb-1, 2.2 fb-1 ttbar+jets 
  theFit.addDataSetFromRootFile("T1", "T1", "/cmsrm/pc21/emanuele/Likelihood2.1.X/PdfsDatasets/ZpWpTTJetsMADGRAPHPdfsDataset.root");
  RooDataSet *data = theFit.getDataSet("T1");

  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");
  
  // Initialize the fit...
  theFit.initialize("MLFit/fitconfig/fit-tagprobe.config");
  
  // Print Fit configuration 
  myPdf->getParameters(data)->selectByAttrib("Constant",kTRUE)->Print("V");  
  myPdf->getParameters(data)->selectByAttrib("Constant",kFALSE)->Print("V");
  
  RooFitResult *fitres =  myPdf->fitTo(*data,theFit.getNoNormVars("myFit"),"MHTER");
  fitres->Print("V");
  
  // write the config file corresponding to the fit minimum
  theFit.writeConfigFile("MLFit/fitres/fitMinimum-tagprobe.config");  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PlotZeeTagAndProbe(int nbins=19) {

  myFit();

  // Load the data
  //  theFit.addDataSetFromRootFile("T1", "T1", "/cmsrm/pc21/emanuele/Likelihood2.1.X/PdfsDatasets/ZJetsMADGRAPHPdfsDataset.root");
  //  theFit.addDataSetFromRootFile("T1", "T1", "/cmsrm/pc21/emanuele/Likelihood2.1.X/PdfsDatasets/WpTTJetsMADGRAPHPdfsDataset.root");
  // W,Z+jets: 300 pb-1, 2.2 fb-1 ttbar+jets
  theFit.addDataSetFromRootFile("T1", "T1", "/cmsrm/pc21/emanuele/Likelihood2.1.X/PdfsDatasets/ZpWpTTJetsMADGRAPHPdfsDataset.root");

  RooDataSet *data = theFit.getDataSet("T1");

  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");

  // Initialize the fit...
  theFit.initialize("MLFit/fitres/fitMinimum-tagprobe.config");

  TCanvas *c = new TCanvas("c","fitResult");
  TFile *output = new TFile("MLFit/fitoutput/data-tagprobe.root","RECREATE");

  RooPlot* MassPlot = MakePlot("zmass", &theFit, data, nbins, false);    

  MassPlot->SetYTitle("Events");
  MassPlot->Draw();
  c->SaveAs("MLFit/fitoutput/data-tagprobe.eps");
  MassPlot->Write();
  //  output->Close();

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make the plot for a given variable
RooPlot *MakePlot(TString VarName, MLFit* theFit, RooDataSet* theData, int nbins, bool poissonError=true)
{
  RooRealVar* Var = theFit->RealObservable(VarName);
  double min=Var->getMin();
  double max=Var->getMax();
  RooPlot *plot = Var->frame(min,max,nbins);
  
  // plot the data
  if(poissonError)
    theData->plotOn(plot);
  else 
    theData->plotOn(plot, RooFit::DataError(RooAbsData::SumW2) );

  // plot the total likelihood
  RooAbsPdf *thePdf = theFit->getPdf("myFit");
  thePdf->plotOn(plot, RooFit::LineColor(kBlack));

  double Ns = theFit->getRealPar("N_sig")->getVal();
  double Nb = theFit->getRealPar("N_bkg")->getVal();

  // plot (dashed) the bkg component
  theFit->getRealPar("N_sig")->setVal(0.);
  thePdf->plotOn(plot, RooFit::Normalization(Nb/(Ns+Nb)),RooFit::LineColor(kBlack),RooFit::LineStyle(kDashed));

  
  return plot;
}



