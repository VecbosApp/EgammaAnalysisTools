// Set Fit Options
MLOptions GetDefaultOptions() {
  MLOptions opts;
  // Fit configuration
  opts.addBoolOption("useinvMass", "Use invMass", kTRUE);        
  opts.addBoolOption("AllFit",          "Fit all species",        kTRUE);
  opts.addBoolOption("ZOnlyFit",        "Fit Z species only",     kFALSE);
  opts.addBoolOption("bkgOnlyFit",      "Fit bkg species only",   kFALSE);
  return opts;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

MLFit theFit;

void myFit() {

  MLFit theFit;

  // Various fit options...
  MLOptions opts = GetDefaultOptions();
  
  // define the structure of the dataset
  RooRealVar *zmass = new RooRealVar("zmass","e^{+}e^{-} Mass [GeV/c^{2}]",90,60,110);
  RooRealVar* weight = new RooRealVar("weight", "weight",1);
  theFit.AddFlatFileColumn(zmass);
  theFit.AddFlatFileColumn(weight);
  
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

  // Various fit options...
  MLOptions opts = GetDefaultOptions();

  // Load the data
  char datasetname[200];
  if(opts.getBoolVal("AllFit")) sprintf(datasetname,"datasets_ZTaP/data.root");
  if(opts.getBoolVal("ZOnlyFit")) sprintf(datasetname,"datasets_ZTaP/zee.root");
  if(opts.getBoolVal("bkgOnlyFit")) sprintf(datasetname,"datasets_ZTaP/background.root");
  theFit.addDataSetFromRootFile("T1", "T1", datasetname);

  RooDataSet *data = theFit.getDataSet("T1");
  data->setWeightVar("weight");

  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");
  
  // Initialize the fit...
  if(opts.getBoolVal("AllFit")) theFit.initialize("fitconfig/fit-tagprobe.config");
  if(opts.getBoolVal("ZOnlyFit")) theFit.initialize("fitconfig/fit-tagprobe-zonly.config");
  if(opts.getBoolVal("bkgOnlyFit")) theFit.initialize("fitconfig/fit-tagprobe-bkgonly.config");
  
  // Print Fit configuration 
  myPdf->getParameters(data)->selectByAttrib("Constant",kTRUE)->Print("V");
  myPdf->getParameters(data)->selectByAttrib("Constant",kFALSE)->Print("V");
  
  RooFitResult *fitres =  myPdf->fitTo(*data,theFit.getNoNormVars("myFit"),"MHTER");
  fitres->Print("V");
  
  // write the config file corresponding to the fit minimum
  if(opts.getBoolVal("AllFit")) theFit.writeConfigFile("fitres/fitMinimum-tagprobe.config");  
  if(opts.getBoolVal("ZOnlyFit")) theFit.writeConfigFile("fitres/fitMinimum-tagprobe-zonly.config");  
  if(opts.getBoolVal("bkgOnlyFit")) theFit.writeConfigFile("fitres/fitMinimum-tagprobe-bkgonly.config");  

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PlotZeeTagAndProbe(int nbins=19) {

  myFit();

  // Various fit options...
  MLOptions opts = GetDefaultOptions();

  // Load the data
  char datasetname[200];
  if(opts.getBoolVal("AllFit")) sprintf(datasetname,"datasets_ZTaP/data.root");
  if(opts.getBoolVal("ZOnlyFit")) sprintf(datasetname,"datasets_ZTaP/zee.root");
  if(opts.getBoolVal("bkgOnlyFit")) sprintf(datasetname,"datasets_ZTaP/background.root");
  theFit.addDataSetFromRootFile("T1", "T1", datasetname);

  RooDataSet *data = theFit.getDataSet("T1");
  data->setWeightVar("weight");

  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");

  // Initialize the fit...
  if(opts.getBoolVal("AllFit")) theFit.initialize("fitres/fitMinimum-tagprobe.config");
  if(opts.getBoolVal("ZOnlyFit")) theFit.initialize("fitres/fitMinimum-tagprobe-zonly.config");
  if(opts.getBoolVal("bkgOnlyFit")) theFit.initialize("fitres/fitMinimum-tagprobe-bkgonly.config");

  TCanvas *c = new TCanvas("c","fitResult");

  RooPlot* MassPlot = MakePlot("zmass", &theFit, data, nbins, false);    

  MassPlot->SetYTitle("electrons in 10 pb^{-1}");
  MassPlot->Draw();
  if(opts.getBoolVal("AllFit")) c->SaveAs("fitoutput/mass-tagprobe.eps");
  if(opts.getBoolVal("ZOnlyFit")) c->SaveAs("fitoutput/mass-tagprobe-zonly.eps");
  if(opts.getBoolVal("bkgOnlyFit")) c->SaveAs("fitoutput/mass-tagprobe-bkgonly.eps");

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make the plot for a given variable
RooPlot *MakePlot(TString VarName, MLFit* theFit, RooDataSet* theData, int nbins, bool poissonError=true)
{

  // Various fit options...
  MLOptions opts = GetDefaultOptions();

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
  thePdf->plotOn(plot, RooFit::LineColor(kRed+3));

  double Ns = theFit->getRealPar("N_sig")->getVal();
  double Nb = theFit->getRealPar("N_bkg")->getVal();

  // plot (dashed) the bkg component
  if(opts.getBoolVal("AllFit")) {
    theFit->getRealPar("N_sig")->setVal(0.);
    thePdf->plotOn(plot, RooFit::Normalization(Nb/(Ns+Nb)),RooFit::LineColor(kAzure+3),RooFit::LineStyle(kDashed));
  }
  
  return plot;
}



