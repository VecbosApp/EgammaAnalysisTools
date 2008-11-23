void calcSWeight_mee() {

  MLFit theFit;
  
  // define the structure of the dataset
  RooRealVar* zmass = new RooRealVar("zmass",  "e^{+}e^{-} Mass [GeV/c^{2}]" , 90., 40., 110.);

  theFit.AddFlatFileColumn(zmass);
  
  // define a fit model
  theFit.addModel("myFit", "Tag And Probe Zee");
  
  // define species
  theFit.addSpecies("myFit", "sig", "Signal Component");
  theFit.addSpecies("myFit", "bkg", "Bkg   Component");
  
  // mLL PDF
  theFit.addPdfWName("myFit", "sig" , "zmass",  "Cruijff",  "sig_Mass");
  theFit.addPdfWName("myFit", "bkg" , "zmass",  "Poly2",  "bkg_Mass");

  theFit.addDataSetFromRootFile("T1", "T1", "/cmsrm/pc21/emanuele/Likelihood2.1.X/PdfsDatasets/ZpWpTTJetsMADGRAPHPdfsDataset.root");

  RooDataSet *data = theFit.getDataSet("T1");
  
  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");
  
  // Initialize the fit...
  theFit.initialize("MLFit/fitres/fitMinimum-tagprobe.config");
  
  // fix all parameters, float the yields and fit
  theFit._parameterSet.selectByName("*")->setAttribAll("Constant",kTRUE);
  (static_cast<RooRealVar*>(theFit._parameterSet.find("N_sig")))->setConstant(kFALSE) ;
  (static_cast<RooRealVar*>(theFit._parameterSet.find("N_bkg")))->setConstant(kFALSE) ;
  
  RooFitResult *fitres =  myPdf->fitTo(*data,theFit.getNoNormVars("myFit"),"MHTER");
  fitres->Print("V");

  // add appropriate column to dataset
  RooArgList yieldsList;
  yieldsList.add(*theFit._fracList.find("N_sig"));
  yieldsList.add(*theFit._fracList.find("N_bkg"));
  cout << "number of entries in set to write: " << data->numEntries() << endl ;
  RooArgSet nonormvars;
  RooDataSet* dsnew = MLSPlot::addSWeightToData((RooSimultaneous*)(myPdf), yieldsList, *data, nonormvars) ;

  TFile sPlots("MLFit/fitoutput/sPlotsTagAndProbe.root","recreate");
  cout << "number of entries in set to write: " << dsnew->numEntries() << endl ;
  dsnew->Write();
  sPlots.Close();
}
