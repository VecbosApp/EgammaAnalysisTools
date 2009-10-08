// Set Fit Options
MLOptions GetDefaultOptions() {
  MLOptions opts;
  // Fit configuration
  opts.addBoolOption("useDeltaPhi",     "Use deltaphi", kTRUE);        
  opts.addBoolOption("useMET",          "Use MET", kFALSE);        
  opts.addBoolOption("useMass",         "Use Mass", kFALSE);        
  opts.addBoolOption("AllFit",          "Fit all species",        kFALSE);
  opts.addBoolOption("QCDOnlyFit",      "Fit QCD species only",   kTRUE);
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
  RooRealVar *deltaphi = new RooRealVar("qcdDeltaphi","di-jet #Delta #phi",1,0.3,TMath::Pi());
  RooRealVar *met = new RooRealVar("qcdMet","MET",10,0,100, "GeV");
  theFit.AddFlatFileColumn(deltaphi);
  theFit.AddFlatFileColumn(met);
  
  // define a fit model
  theFit.addModel("myFit", "Tag And Probe Zee");
  
  // define species
  theFit.addSpecies("myFit", "sig", "Signal Component");
  theFit.addSpecies("myFit", "bkg", "Bkg   Component");

  // deltaphi PDF
  if(opts.getBoolVal("useDeltaPhi")) {
    theFit.addPdfWName("myFit", "sig" , "qcdDeltaphi",  "Poly3", "sig_deltaphi");
    theFit.addPdfWName("myFit", "bkg" , "qcdDeltaphi",  "Poly3", "bkg_deltaphi");
  }
  
  // MET PDF
  if(opts.getBoolVal("useMET")) {
    theFit.addPdfWName("myFit", "sig" , "qcdMet",  "Cruijff", "sig_met");
    theFit.addPdfWName("myFit", "bkg" , "qcdMet",  "Cruijff",  "bkg_met");
  }
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Fit a sample of Z events
void FitQCDTagAndProbe() {
  
  myFit();

  // Various fit options...
  MLOptions opts = GetDefaultOptions();

  // Load the data
  char datasetname[200];
  if(opts.getBoolVal("AllFit")) sprintf(datasetname,"datasets/merged_data.root");
  if(opts.getBoolVal("QCDOnlyFit")) sprintf(datasetname,"datasets/qcd.root");
  if(opts.getBoolVal("bkgOnlyFit")) sprintf(datasetname,"datasets/qcd-bkg.root");
  theFit.addDataSetFromRootFile("T1", "T1", datasetname);
  RooDataSet *data = theFit.getDataSet("T1");

  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");
  
  // Initialize the fit...
  if(opts.getBoolVal("AllFit")) theFit.initialize("fitconfig/fit-qcdtagprobe.config");
  if(opts.getBoolVal("QCDOnlyFit")) theFit.initialize("fitconfig/fit-qcdtagprobe-qcdonly.config");
  if(opts.getBoolVal("bkgOnlyFit")) theFit.initialize("fitconfig/fit-qcdtagprobe-bkgonly.config");
  
  // Print Fit configuration 
  myPdf->getParameters(data)->selectByAttrib("Constant",kTRUE)->Print("V");
  myPdf->getParameters(data)->selectByAttrib("Constant",kFALSE)->Print("V");
  
  RooFitResult *fitres =  myPdf->fitTo(*data,theFit.getNoNormVars("myFit"),"MHTER");
  fitres->Print("V");
  
  // write the config file corresponding to the fit minimum
  if(opts.getBoolVal("AllFit")) theFit.writeConfigFile("fitres/fitMinimum-qcdtagprobe.config");  
  if(opts.getBoolVal("QCDOnlyFit")) theFit.writeConfigFile("fitres/fitMinimum-qcdtagprobe-qcdonly.config");  
  if(opts.getBoolVal("bkgOnlyFit")) theFit.writeConfigFile("fitres/fitMinimum-qcdtagprobe-bkgonly.config");  

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PlotQCDTagAndProbe(int nbins=19) {

  myFit();

  // Various fit options...
  MLOptions opts = GetDefaultOptions();

  // Load the data
  char datasetname[200];
  if(opts.getBoolVal("AllFit")) sprintf(datasetname,"datasets/merged_data.root");
  if(opts.getBoolVal("QCDOnlyFit")) sprintf(datasetname,"datasets/qcd.root");
  if(opts.getBoolVal("bkgOnlyFit")) sprintf(datasetname,"datasets/qcd-bkg.root");
  theFit.addDataSetFromRootFile("T1", "T1", datasetname);

  RooDataSet *data = theFit.getDataSet("T1");

  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");

  // Initialize the fit...
  if(opts.getBoolVal("AllFit")) theFit.initialize("fitres/fitMinimum-qcdtagprobe.config");
  if(opts.getBoolVal("QCDOnlyFit")) theFit.initialize("fitres/fitMinimum-qcdtagprobe-qcdonly.config");
  if(opts.getBoolVal("bkgOnlyFit")) theFit.initialize("fitres/fitMinimum-qcdtagprobe-bkgonly.config");

  if(opts.getBoolVal("useDeltaPhi")) {
    TCanvas *c = new TCanvas("c","fitResult");

    RooPlot* MassPlot = MakePlot("qcdDeltaphi", &theFit, data, nbins, false);    
    
    MassPlot->SetYTitle("Events");
    MassPlot->Draw();
    if(opts.getBoolVal("AllFit")) c->SaveAs("fitoutput/deltaphi-data-qcdtagprobe.eps");
    if(opts.getBoolVal("QCDOnlyFit")) c->SaveAs("fitoutput/deltaphi-qcdtagprobe-qcdonly.eps");
    if(opts.getBoolVal("bkgOnlyFit")) c->SaveAs("fitoutput/deltaphi-qcdtagprobe-bkgonly.eps");
  }
  if(opts.getBoolVal("useMET")) {
    TCanvas *c = new TCanvas("c","fitResult");

    RooPlot* MassPlot = MakePlot("qcdMet", &theFit, data, nbins, false);    
    
    MassPlot->SetYTitle("Events");
    MassPlot->Draw();
    if(opts.getBoolVal("AllFit")) c->SaveAs("fitoutput/met-data-qcdtagprobe.eps");
    if(opts.getBoolVal("QCDOnlyFit")) c->SaveAs("fitoutput/met-qcdtagprobe-qcdonly.eps");
    if(opts.getBoolVal("bkgOnlyFit")) c->SaveAs("fitoutput/met-qcdtagprobe-bkgonly.eps");
  }

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



