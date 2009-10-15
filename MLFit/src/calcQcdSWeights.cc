void calcQcdSWeight_metDeltaphi() {

  MLFit theFit;
  
  // define the structure of the dataset
  RooRealVar *deltaphi = new RooRealVar("qcdDeltaphi","di-jet #Delta #phi",     1,1.0,3.1415);
  RooRealVar *met      = new RooRealVar("qcdMet",     "MET",                    10,0,100, "GeV");
  RooRealVar *weight   = new RooRealVar("weight",     "weight",                 1);
  theFit.AddFlatFileColumn(deltaphi);
  theFit.AddFlatFileColumn(met);
  theFit.AddFlatFileColumn(weight);

  // define a fit model
  theFit.addModel("myFit", "Tag And Probe Zee");

  // define species                                                                                                                                   
  theFit.addSpecies("myFit", "sig",    "Signal Component");
  theFit.addSpecies("myFit", "wenu",   "Wenu   Component");
  theFit.addSpecies("myFit", "ttbar",  "TTbar  Component");
  
  // deltaphi PDF                                                                                                                                     
  theFit.addPdfWName("myFit", "sig" ,    "qcdDeltaphi",  "Cruijff", "sig_deltaphi");
  theFit.addPdfWName("myFit", "wenu" ,   "qcdDeltaphi",  "Poly2", "wenu_deltaphi");
  theFit.addPdfWName("myFit", "ttbar" ,  "qcdDeltaphi",  "Poly3",   "ttbar_deltaphi");

  // MET PDF                                                                                                                                          
  theFit.addPdfWName("myFit", "sig" ,        "qcdMet",  "Cruijff",  "sig_met");
  theFit.addPdfWName("myFit", "wenu" ,       "qcdMet",  "Cruijff",  "wenu_met");
  theFit.addPdfWName("myFit", "ttbar" ,      "qcdMet",  "Cruijff",  "ttbar_met");

  theFit.addDataSetFromRootFile("T1", "T1", "/cmsrm/pc21/crovelli/data/Like3.2.X/datasets_QCDTaP/data.root");

  RooDataSet *data = theFit.getDataSet("T1");
  
  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");
  
  // Initialize the fit...
  theFit.initialize("fitres/fitMinimum-qcdtagprobe.config");
  
  // fix all parameters, float the yields and fit
  theFit._parameterSet.selectByName("*")->setAttribAll("Constant",kTRUE);
  (static_cast<RooRealVar*>(theFit._parameterSet.find("N_sig")))->setConstant(kFALSE) ;
  //  (static_cast<RooRealVar*>(theFit._parameterSet.find("N_ttbar")))->setConstant(kFALSE) ;
  (static_cast<RooRealVar*>(theFit._parameterSet.find("N_wenu")))->setConstant(kFALSE) ;
  
  RooFitResult *fitres =  myPdf->fitTo(*data,theFit.getNoNormVars("myFit"),"MHTER");
  fitres->Print("V");

  // add appropriate column to dataset
  RooArgList yieldsList;
  yieldsList.add(*theFit._fracList.find("N_sig"));
  //  yieldsList.add(*theFit._fracList.find("N_ttbar"));
  yieldsList.add(*theFit._fracList.find("N_wenu"));
  cout << "number of entries in set to write: " << data->numEntries() << endl ;
  RooArgSet nonormvars;
  RooDataSet* dsnew = MLSPlot::addSWeightToData((RooSimultaneous*)(myPdf), yieldsList, *data, nonormvars) ;

  TFile sPlots("/cmsrm/pc21/crovelli/data/Like3.2.X/sPlots/qcdTaP.root","recreate");
  cout << "number of entries in set to write: " << dsnew->numEntries() << endl ;
  dsnew->Write();
  sPlots.Close();
}
