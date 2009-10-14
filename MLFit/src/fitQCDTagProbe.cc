// Set Fit Options
MLOptions GetDefaultOptions() {
  MLOptions opts;
  // Fit configuration
  opts.addBoolOption("useDeltaPhi",      "Use deltaphi",            kTRUE);        
  opts.addBoolOption("useMET",           "Use MET",                 kFALSE);        
  opts.addBoolOption("useMass",          "Use Mass",                kFALSE);        
  opts.addBoolOption("AllFit",           "Fit all species",         kFALSE);
  opts.addBoolOption("QCDOnlyFit",       "Fit QCD species only",    kFALSE);
  opts.addBoolOption("wenuOnlyFit",      "Fit W->wnu species only", kFALSE);
  opts.addBoolOption("ttbarOnlyFit",     "Fit TTbar species only",  kFALSE);
  opts.addBoolOption("zeeOnlyFit",       "Fit Z species only",      kTRUE);
  // opts.addBoolOption("zAndGammaJOnlyFit", "Fit Z and gamma+jet species only", kTRUE);
  return opts;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

MLFit theFit;

void myFit() {

  MLFit theFit;

  // Various fit options...
  MLOptions opts = GetDefaultOptions();

  // define the structure of the dataset
  RooRealVar *deltaphi = new RooRealVar("qcdDeltaphi","di-jet #Delta #phi",     1,0.5,3.1415);
  RooRealVar *met      = new RooRealVar("qcdMet",     "MET",                    10,0,100, "GeV");
  RooRealVar *mass     = new RooRealVar("qcdInvmass", "tag-probeinvariant mass",10,0,600, "GeV");
  // RooRealVar *weight   = new RooRealVar("weight",     "weight",                 1);
  theFit.AddFlatFileColumn(deltaphi);
  theFit.AddFlatFileColumn(met);
  theFit.AddFlatFileColumn(mass);
  // theFit.AddFlatFileColumn(weight);
  
  // define a fit model
  theFit.addModel("myFit", "Tag And Probe Zee");
  
  // define species
  theFit.addSpecies("myFit", "sig",        "Signal Component");
  theFit.addSpecies("myFit", "wenu",       "Wenu   Component");
  theFit.addSpecies("myFit", "ttbar",      "TTbar  Component");
  // theFit.addSpecies("myFit", "zAndGammaJ", "Zee and photon + jet   Component");
  theFit.addSpecies("myFit", "zee",        "Zee Component");

  // deltaphi PDF
  if(opts.getBoolVal("useDeltaPhi")) {
    theFit.addPdfWName("myFit", "sig" ,        "qcdDeltaphi",  "Poly3",   "sig_deltaphi");
    theFit.addPdfWName("myFit", "wenu" ,       "qcdDeltaphi",  "Cruijff", "wenu_deltaphi");
    theFit.addPdfWName("myFit", "ttbar" ,      "qcdDeltaphi",  "Poly3",   "ttbar_deltaphi");
    // theFit.addPdfWName("myFit", "zee" ,        "qcdDeltaphi",  "Poly3",   "zee_deltaphi");
    theFit.addPdfWName("myFit", "zee" ,        "qcdDeltaphi",  "Cruijff",  "zee_deltaphi");
    // theFit.addPdfWName("myFit", "zAndGammaJ" , "qcdDeltaphi",  "Poly3", "zAndGammaJ_deltaphi");
  }
  
  // MET PDF
  if(opts.getBoolVal("useMET")) {
    theFit.addPdfWName("myFit", "sig" ,        "qcdMet",  "Cruijff",  "sig_met");
    theFit.addPdfWName("myFit", "wenu" ,       "qcdMet",  "Cruijff",  "wenu_met");
    theFit.addPdfWName("myFit", "ttbar" ,      "qcdMet",  "Cruijff",  "ttbar_met");
    theFit.addPdfWName("myFit", "zee" ,        "qcdMet",  "Cruijff",  "zee_met");
    // theFit.addPdfWName("myFit", "zAndGammaJ" , "qcdMet",  "Cruijff",  "zAndGammaJ_met");
  }  

  // mass PDF
  if(opts.getBoolVal("useMass")) {
    theFit.addPdfWName("myFit", "sig" ,        "qcdInvmass",  "Cruijff",  "sig_mass");
    theFit.addPdfWName("myFit", "wenu" ,       "qcdInvmass",  "Cruijff",  "wenu_mass");
    theFit.addPdfWName("myFit", "ttbar" ,      "qcdInvmass",  "Cruijff",  "ttbar_mass");
    theFit.addPdfWName("myFit", "zee" ,        "qcdInvmass",  "Cruijff",  "zee_mass");
    // theFit.addPdfWName("myFit", "zAndGammaJ" , "qcdInvmass",  "Cruijff",  "zAndGammaJ_mass");
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
  if(opts.getBoolVal("QCDOnlyFit"))         sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/qcd_pt80.root");
  if(opts.getBoolVal("wenuOnlyFit"))        sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/wenu_5jobs.root");
  if(opts.getBoolVal("ttbarOnlyFit"))       sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/ttbar.root");
  if(opts.getBoolVal("zeeOnlyFit"))         sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/zee.root");
  // if(opts.getBoolVal("zAndGammaJOnlyFit"))  sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/zAndPhotonJet.root");
  theFit.addDataSetFromRootFile("T1", "T1", datasetname);
  RooDataSet *data = theFit.getDataSet("T1");
  // data->setWeightVar("weight");
  
  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");
  
  // Initialize the fit...
  if(opts.getBoolVal("AllFit"))            theFit.initialize("fitconfig/fit-qcdtagprobe.config");
  if(opts.getBoolVal("QCDOnlyFit"))        theFit.initialize("fitconfig/fit-qcdtagprobe-qcdonly.config");
  if(opts.getBoolVal("wenuOnlyFit"))       theFit.initialize("fitconfig/fit-qcdtagprobe-wenuonly.config");
  if(opts.getBoolVal("ttbarOnlyFit"))      theFit.initialize("fitconfig/fit-qcdtagprobe-ttbaronly.config");
  if(opts.getBoolVal("zeeOnlyFit"))        theFit.initialize("fitconfig/fit-qcdtagprobe-zeeonly.config");
  // if(opts.getBoolVal("zAndGammaJOnlyFit")) theFit.initialize("fitconfig/fit-qcdtagprobe-zAndGammaJonly.config");
  
  // Print Fit configuration 
  myPdf->getParameters(data)->selectByAttrib("Constant",kTRUE)->Print("V");
  myPdf->getParameters(data)->selectByAttrib("Constant",kFALSE)->Print("V");
  
  RooFitResult *fitres =  myPdf->fitTo(*data,theFit.getNoNormVars("myFit"),"MHTER");
  fitres->Print("V");
  
  // write the config file corresponding to the fit minimum
  if(opts.getBoolVal("AllFit"))            theFit.writeConfigFile("fitres/fitMinimum-qcdtagprobe.config");  
  if(opts.getBoolVal("QCDOnlyFit"))        theFit.writeConfigFile("fitres/fitMinimum-qcdtagprobe-qcdonly.config");  
  if(opts.getBoolVal("wenuOnlyFit"))       theFit.writeConfigFile("fitres/fitMinimum-qcdtagprobe-wenuonly.config");  
  if(opts.getBoolVal("ttbarOnlyFit"))      theFit.writeConfigFile("fitres/fitMinimum-qcdtagprobe-ttbaronly.config");  
  if(opts.getBoolVal("zeeOnlyFit"))        theFit.writeConfigFile("fitres/fitMinimum-qcdtagprobe-zeeonly.config");  
  // if(opts.getBoolVal("zAndGammaJOnlyFit")) theFit.writeConfigFile("fitres/fitMinimum-qcdtagprobe-zAndGammaJonly.config");  

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PlotQCDTagAndProbe(int nbins=19) {

  myFit();

  // Various fit options...
  MLOptions opts = GetDefaultOptions();

  // Load the data
  char datasetname[200];
  if(opts.getBoolVal("AllFit"))             sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/data.root");
  if(opts.getBoolVal("QCDOnlyFit"))         sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/qcd_pt80.root");
  if(opts.getBoolVal("wenuOnlyFit"))        sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/wenu_5jobs.root");
  if(opts.getBoolVal("ttbarOnlyFit"))       sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/ttbar.root");
  if(opts.getBoolVal("zeeOnlyFit"))         sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/zee.root");
  // if(opts.getBoolVal("zAndGammaJOnlyFit"))  sprintf(datasetname,"/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/zAndPhotonJet.root");
  theFit.addDataSetFromRootFile("T1", "T1", datasetname);

  RooDataSet *data = theFit.getDataSet("T1");
  // data->setWeightVar("weight");

  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");

  // Initialize the fit...
  if(opts.getBoolVal("AllFit"))            theFit.initialize("fitres/fitMinimum-qcdtagprobe.config");  
  if(opts.getBoolVal("QCDOnlyFit"))        theFit.initialize("fitres/fitMinimum-qcdtagprobe-qcdonly.config");  
  if(opts.getBoolVal("wenuOnlyFit"))       theFit.initialize("fitres/fitMinimum-qcdtagprobe-wenuonly.config");  
  if(opts.getBoolVal("ttbarOnlyFit"))      theFit.initialize("fitres/fitMinimum-qcdtagprobe-ttbaronly.config");  
  if(opts.getBoolVal("zeeOnlyFit"))        theFit.initialize("fitres/fitMinimum-qcdtagprobe-zeeonly.config");  
  // if(opts.getBoolVal("zAndGammaJOnlyFit")) theFit.initialize("fitres/fitMinimum-qcdtagprobe-zAndGammaJonly.config");  

  if(opts.getBoolVal("useDeltaPhi")) {
    TCanvas *c = new TCanvas("c","fitResult");

    RooPlot* MassPlot = MakePlot("qcdDeltaphi", &theFit, data, nbins, false);    
    
    MassPlot->SetYTitle("Events");
    MassPlot->Draw();
    if(opts.getBoolVal("AllFit"))            c->SaveAs("fitoutput/deltaphi-data-qcdtagprobe.eps");
    if(opts.getBoolVal("QCDOnlyFit"))        c->SaveAs("fitoutput/deltaphi-qcdtagprobe-qcdonly.eps");
    if(opts.getBoolVal("wenuOnlyFit"))       c->SaveAs("fitoutput/deltaphi-qcdtagprobe-wenuonly.eps");
    if(opts.getBoolVal("ttbarOnlyFit"))      c->SaveAs("fitoutput/deltaphi-qcdtagprobe-ttbaronly.eps");
    if(opts.getBoolVal("zeeOnlyFit"))        c->SaveAs("fitoutput/deltaphi-qcdtagprobe-zeeonly.eps");
    // if(opts.getBoolVal("zAndGammaJOnlyFit")) c->SaveAs("fitoutput/deltaphi-qcdtagprobe-zAndGammaJonly.eps");
  }

  if(opts.getBoolVal("useMET")) {
    TCanvas *c = new TCanvas("c","fitResult");

    RooPlot* MassPlot = MakePlot("qcdMet", &theFit, data, nbins, false);    
    
    MassPlot->SetYTitle("Events");
    MassPlot->Draw();
    if(opts.getBoolVal("AllFit"))            c->SaveAs("fitoutput/met-data-qcdtagprobe.eps");
    if(opts.getBoolVal("QCDOnlyFit"))        c->SaveAs("fitoutput/met-qcdtagprobe-qcdonly.eps");
    if(opts.getBoolVal("wenuOnlyFit"))       c->SaveAs("fitoutput/met-qcdtagprobe-wenuonly.eps");
    if(opts.getBoolVal("ttbarOnlyFit"))      c->SaveAs("fitoutput/met-qcdtagprobe-ttbaronly.eps");
    if(opts.getBoolVal("zeeOnlyFit"))        c->SaveAs("fitoutput/met-qcdtagprobe-zeeonly.eps");
    // if(opts.getBoolVal("zAndGammaJOnlyFit")) c->SaveAs("fitoutput/met-qcdtagprobe-zAndGammaJonly.eps");
  }

  if(opts.getBoolVal("useMass")) {

    TCanvas *c = new TCanvas("c","fitResult");

    RooPlot* MassPlot = MakePlot("qcdInvmass", &theFit, data, nbins, false);    
    
    MassPlot->SetYTitle("Events");
    MassPlot->Draw();
    if(opts.getBoolVal("AllFit"))            c->SaveAs("fitoutput/mass-data-qcdtagprobe.eps");
    if(opts.getBoolVal("QCDOnlyFit"))        c->SaveAs("fitoutput/mass-qcdtagprobe-qcdonly.eps");
    if(opts.getBoolVal("wenuOnlyFit"))       c->SaveAs("fitoutput/mass-qcdtagprobe-wenuonly.eps");
    if(opts.getBoolVal("ttbarOnlyFit"))      c->SaveAs("fitoutput/mass-qcdtagprobe-ttbaronly.eps");
    if(opts.getBoolVal("zeeOnlyFit"))        c->SaveAs("fitoutput/mass-qcdtagprobe-zeeonly.eps");
    //if(opts.getBoolVal("zAndGammaJOnlyFit")) c->SaveAs("fitoutput/mass-qcdtagprobe-zAndGammaJonly.eps");
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

  double Ns          = theFit->getRealPar("N_sig")->getVal();
  double Nwenu       = theFit->getRealPar("N_wenu")->getVal();
  double Nttbar      = theFit->getRealPar("N_ttbar")->getVal();
  // double NzAndGammaJ = theFit->getRealPar("N_zAndGammaJ")->getVal();
  double Nzee        = theFit->getRealPar("N_zee")->getVal();
  // double Nb = Nwenu + Nttbar + NzAndGammaJ;
  double Nb = Nwenu + Nttbar + Nzee;

  // plot (dashed) the bkg component
  theFit->getRealPar("N_sig")->setVal(0.);
  thePdf->plotOn(plot, RooFit::Normalization(Nb/(Ns+Nb)),RooFit::LineColor(kBlack),RooFit::LineStyle(kDashed));

  
  return plot;
}



