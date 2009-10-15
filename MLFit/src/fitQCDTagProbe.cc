// Set Fit Options
MLOptions GetDefaultOptions() {
  MLOptions opts;
  // Fit configuration
  opts.addBoolOption("useDeltaPhi",      "Use deltaphi",            kTRUE);        
  opts.addBoolOption("useMET",           "Use MET",                 kTRUE);        
  opts.addBoolOption("useMass",          "Use Mass",                kFALSE);        
  opts.addBoolOption("AllFit",           "Fit all species",         kFALSE);
  opts.addBoolOption("QCDOnlyFit",       "Fit QCD species only",    kFALSE);
  opts.addBoolOption("wenuOnlyFit",      "Fit W->wnu species only", kFALSE);
  opts.addBoolOption("ttbarOnlyFit",     "Fit TTbar species only",  kTRUE);
  return opts;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

MLFit theFit;

void myFit() {

  MLFit theFit;

  // Various fit options...
  MLOptions opts = GetDefaultOptions();

  // define the structure of the dataset
  RooRealVar *deltaphi = new RooRealVar("qcdDeltaphi","di-jet #Delta #phi",     1,1.0,3.1415);
  RooRealVar *met      = new RooRealVar("qcdMet",     "MET",                    10,0,100, "GeV");  
  RooRealVar *mass     = new RooRealVar("qcdInvmass", "tag-probeinvariant mass",10,0,600, "GeV");
  RooRealVar *weight   = new RooRealVar("weight",     "weight",                 1);
  theFit.AddFlatFileColumn(deltaphi);
  theFit.AddFlatFileColumn(met);
  theFit.AddFlatFileColumn(mass);      
  theFit.AddFlatFileColumn(weight);
  
  // define a fit model
  theFit.addModel("myFit", "Tag And Probe Zee");
  
  // define species
  theFit.addSpecies("myFit", "sig",    "Signal Component");
  theFit.addSpecies("myFit", "wenu",   "Wenu   Component");
  theFit.addSpecies("myFit", "ttbar",  "TTbar  Component");

  // deltaphi PDF
  if(opts.getBoolVal("useDeltaPhi")) {
    theFit.addPdfWName("myFit", "sig" ,    "qcdDeltaphi",  "Cruijff", "sig_deltaphi");
    theFit.addPdfWName("myFit", "wenu" ,   "qcdDeltaphi",  "Poly2", "wenu_deltaphi");
    theFit.addPdfWName("myFit", "ttbar" ,  "qcdDeltaphi",  "Poly3",   "ttbar_deltaphi");
  }
  
  // MET PDF
  if(opts.getBoolVal("useMET")) {
    theFit.addPdfWName("myFit", "sig" ,        "qcdMet",  "Cruijff",  "sig_met");
    theFit.addPdfWName("myFit", "wenu" ,       "qcdMet",  "Cruijff",  "wenu_met");
    theFit.addPdfWName("myFit", "ttbar" ,      "qcdMet",  "Cruijff",  "ttbar_met");
  }  

  // mass PDF
  if(opts.getBoolVal("useMass")) {
    theFit.addPdfWName("myFit", "sig" ,        "qcdInvmass",  "Cruijff",  "sig_mass");
    theFit.addPdfWName("myFit", "wenu" ,       "qcdInvmass",  "Cruijff",  "wenu_mass");
    theFit.addPdfWName("myFit", "ttbar" ,      "qcdInvmass",  "Cruijff",  "ttbar_mass");
  }  
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Fit a sample of QCD events
void FitQCDTagAndProbe() {
  
  myFit();
  
  // Various fit options...
  MLOptions opts = GetDefaultOptions();

  // Load the data
  char datasetname[200];
  if(opts.getBoolVal("AllFit"))             sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/data.root");
  if(opts.getBoolVal("QCDOnlyFit"))         sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/qcdSignal.root"); // QCD+gammaPlusJet 
  if(opts.getBoolVal("wenuOnlyFit"))        sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/wenu.root");
  if(opts.getBoolVal("ttbarOnlyFit"))       sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/ttbar.root");
  theFit.addDataSetFromRootFile("T1", "T1", datasetname);
  RooDataSet *data = theFit.getDataSet("T1");
  RooDataSet *reddata = data->reduce("qcdDeltaphi>1.0 && qcdDeltaphi<3.1415");
  reddata->setWeightVar("weight");
  
  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");
  
  // Initialize the fit...
  if(opts.getBoolVal("AllFit"))            theFit.initialize("fitconfig/fit-qcdtagprobe.config");
  if(opts.getBoolVal("QCDOnlyFit"))        theFit.initialize("fitconfig/fit-qcdtagprobe-qcdonly.config");
  if(opts.getBoolVal("wenuOnlyFit"))       theFit.initialize("fitconfig/fit-qcdtagprobe-wenuonly.config");
  if(opts.getBoolVal("ttbarOnlyFit"))      theFit.initialize("fitconfig/fit-qcdtagprobe-ttbaronly.config");
  
  // Print Fit configuration 
  myPdf->getParameters(reddata)->selectByAttrib("Constant",kTRUE)->Print("V");
  myPdf->getParameters(reddata)->selectByAttrib("Constant",kFALSE)->Print("V");
  
  RooFitResult *fitres =  myPdf->fitTo(*reddata,theFit.getNoNormVars("myFit"),"MHTER");
  fitres->Print("V");
  
  // write the config file corresponding to the fit minimum
  if(opts.getBoolVal("AllFit"))            theFit.writeConfigFile("fitres/fitMinimum-qcdtagprobe.config");  
  if(opts.getBoolVal("QCDOnlyFit"))        theFit.writeConfigFile("fitres/fitMinimum-qcdtagprobe-qcdonly.config");  
  if(opts.getBoolVal("wenuOnlyFit"))       theFit.writeConfigFile("fitres/fitMinimum-qcdtagprobe-wenuonly.config");  
  if(opts.getBoolVal("ttbarOnlyFit"))      theFit.writeConfigFile("fitres/fitMinimum-qcdtagprobe-ttbaronly.config");  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PlotQCDTagAndProbe(int nbins=19) {

  myFit();

  // Various fit options...
  MLOptions opts = GetDefaultOptions();

  // Load the data
  char datasetname[200];
  if(opts.getBoolVal("AllFit"))             sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/data.root");
  if(opts.getBoolVal("QCDOnlyFit"))         sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/qcdSignal.root");
  if(opts.getBoolVal("wenuOnlyFit"))        sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/wenu.root"); 
  if(opts.getBoolVal("ttbarOnlyFit"))       sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/ttbar.root");
  theFit.addDataSetFromRootFile("T1", "T1", datasetname);

  RooDataSet *data = theFit.getDataSet("T1");
  RooDataSet *reddata = data->reduce("qcdDeltaphi>1.0 && qcdDeltaphi<3.1415");
  reddata->setWeightVar("weight");

  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");

  // Initialize the fit...
  if(opts.getBoolVal("AllFit"))            theFit.initialize("fitres/fitMinimum-qcdtagprobe.config");  
  if(opts.getBoolVal("QCDOnlyFit"))        theFit.initialize("fitres/fitMinimum-qcdtagprobe-qcdonly.config");  
  if(opts.getBoolVal("wenuOnlyFit"))       theFit.initialize("fitres/fitMinimum-qcdtagprobe-wenuonly.config");  
  if(opts.getBoolVal("ttbarOnlyFit"))      theFit.initialize("fitres/fitMinimum-qcdtagprobe-ttbaronly.config");  

  if(opts.getBoolVal("useDeltaPhi")) {
    TCanvas *c = new TCanvas("c","fitResult");

    RooPlot* MassPlot = MakePlot("qcdDeltaphi", &theFit, reddata, nbins, false);    
    
    MassPlot->SetYTitle("Events");
    MassPlot->Draw();
    if(opts.getBoolVal("AllFit"))            c->SaveAs("fitoutput/deltaphi-data-qcdtagprobe.eps");
    if(opts.getBoolVal("QCDOnlyFit"))        c->SaveAs("fitoutput/deltaphi-qcdtagprobe-qcdonly.eps");
    if(opts.getBoolVal("wenuOnlyFit"))       c->SaveAs("fitoutput/deltaphi-qcdtagprobe-wenuonly.eps");
    if(opts.getBoolVal("ttbarOnlyFit"))      c->SaveAs("fitoutput/deltaphi-qcdtagprobe-ttbaronly.eps");
  }

  if(opts.getBoolVal("useMET")) {
    TCanvas *c = new TCanvas("c","fitResult");

    RooPlot* MassPlot = MakePlot("qcdMet", &theFit, reddata, nbins, false);    
    
    MassPlot->SetYTitle("Events");
    MassPlot->Draw();
    if(opts.getBoolVal("AllFit"))            c->SaveAs("fitoutput/met-data-qcdtagprobe.eps");
    if(opts.getBoolVal("QCDOnlyFit"))        c->SaveAs("fitoutput/met-qcdtagprobe-qcdonly.eps");
    if(opts.getBoolVal("wenuOnlyFit"))       c->SaveAs("fitoutput/met-qcdtagprobe-wenuonly.eps");
    if(opts.getBoolVal("ttbarOnlyFit"))      c->SaveAs("fitoutput/met-qcdtagprobe-ttbaronly.eps");
  }

  if(opts.getBoolVal("useMass")) {

    TCanvas *c = new TCanvas("c","fitResult");

    RooPlot* MassPlot = MakePlot("qcdInvmass", &theFit, reddata, nbins, false);    
    
    MassPlot->SetYTitle("Events");
    MassPlot->Draw();
    if(opts.getBoolVal("AllFit"))            c->SaveAs("fitoutput/mass-data-qcdtagprobe.eps");
    if(opts.getBoolVal("QCDOnlyFit"))        c->SaveAs("fitoutput/mass-qcdtagprobe-qcdonly.eps");
    if(opts.getBoolVal("wenuOnlyFit"))       c->SaveAs("fitoutput/mass-qcdtagprobe-wenuonly.eps");
    if(opts.getBoolVal("ttbarOnlyFit"))      c->SaveAs("fitoutput/mass-qcdtagprobe-ttbaronly.eps");
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

  double Ns     = theFit->getRealPar("N_sig")->getVal();
  double Nwenu  = theFit->getRealPar("N_wenu")->getVal();
  double Nttbar = theFit->getRealPar("N_ttbar")->getVal();
  double Nb     = Nwenu + Nttbar;

  // plot (dashed) the bkg component
  theFit->getRealPar("N_sig")->setVal(0.);
  thePdf->plotOn(plot, RooFit::Normalization(Nb/(Ns+Nb)),RooFit::LineColor(kBlack),RooFit::LineStyle(kDashed));

  
  return plot;
}



