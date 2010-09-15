{
// style
TStyle *tesiStyle = new TStyle("tesiStyle","");
tesiStyle->SetCanvasColor(0);
tesiStyle->SetFrameFillColor(0);
tesiStyle->SetStatColor(0);
tesiStyle->SetOptStat(0);
tesiStyle->SetTitleFillColor(0);
tesiStyle->SetCanvasBorderMode(0);
tesiStyle->SetPadBorderMode(0);
tesiStyle->SetFrameBorderMode(0);
tesiStyle->cd();
 
 TFile fileWcandle("pdfsQCD_data.root");
 TFile fileJets   ("pdfsQCD_data.root");

 TH1F *dPhiUnsplitEle[2][2][2];   // [2][2][2] = [EB,EE][pt<20; pt>20][Wcandle,Jet]
 TH1F *dEtaUnsplitEle[2][2][2];
 TH1F *EoPUnsplitEle[2][2][2];
 TH1F *HoEUnsplitEle[2][2][2];
 TH1F *sigmaIEtaIEtaUnsplitEle[2][2][2];
 TH1F *sigmaIPhiIPhiUnsplitEle[2][2][2];
 TH1F *fBremUnsplitEle[2][2][2];


 // taking histos
 char buffer[400];
 for (int iecal=0; iecal<2; iecal++) {
   for(int iptbin=0; iptbin<2; iptbin++) {

     sprintf(buffer,"dPhiUnsplit_hadrons_%d_%d",iecal,iptbin);
     std::cout << buffer << std::endl;
     dPhiUnsplitEle[iecal][iptbin][0] = (TH1F*)fileWcandle->Get(buffer);
     dPhiUnsplitEle[iecal][iptbin][1] = (TH1F*)fileJets->Get(buffer);

     sprintf(buffer,"dEtaUnsplit_hadrons_%d_%d",iecal,iptbin);
     dEtaUnsplitEle[iecal][iptbin][0] = (TH1F*)fileWcandle->Get(buffer);
     dEtaUnsplitEle[iecal][iptbin][1] = (TH1F*)fileJets->Get(buffer);
     
     sprintf(buffer,"EoPUnsplit_hadrons_%d_%d",iecal,iptbin);
     EoPUnsplitEle[iecal][iptbin][0] = (TH1F*)fileWcandle->Get(buffer);
     EoPUnsplitEle[iecal][iptbin][1] = (TH1F*)fileJets->Get(buffer);

     sprintf(buffer,"HoEUnsplit_hadrons_%d_%d",iecal,iptbin);
     HoEUnsplitEle[iecal][iptbin][0] = (TH1F*)fileWcandle->Get(buffer);
     HoEUnsplitEle[iecal][iptbin][1] = (TH1F*)fileJets->Get(buffer);

     sprintf(buffer,"sigmaIEtaIEtaUnsplit_hadrons_%d_%d",iecal,iptbin);
     sigmaIEtaIEtaUnsplitEle[iecal][iptbin][0] = (TH1F*)fileWcandle->Get(buffer);
     sigmaIEtaIEtaUnsplitEle[iecal][iptbin][1] = (TH1F*)fileJets->Get(buffer);

     sprintf(buffer,"sigmaIPhiIPhiUnsplit_hadrons_%d_%d",iecal,iptbin);
     sigmaIPhiIPhiUnsplitEle[iecal][iptbin][0] = (TH1F*)fileWcandle->Get(buffer);
     sigmaIPhiIPhiUnsplitEle[iecal][iptbin][1] = (TH1F*)fileJets->Get(buffer);

     sprintf(buffer,"fBremUnsplit_hadrons_%d_%d",iecal,iptbin);
     fBremUnsplitEle[iecal][iptbin][0] = (TH1F*)fileWcandle->Get(buffer);
     fBremUnsplitEle[iecal][iptbin][1] = (TH1F*)fileJets->Get(buffer);
   }
 }


 // scaling 
 for (int iecal=0; iecal<2; iecal++) {
   for(int iptbin=0; iptbin<2; iptbin++) {
     for(int isample=0;isample<2;isample++) {
       dPhiUnsplitEle[iecal][iptbin][isample]  -> Sumw2();
       dEtaUnsplitEle[iecal][iptbin][isample]  -> Sumw2();
       EoPUnsplitEle[iecal][iptbin][isample]   -> Sumw2();
       HoEUnsplitEle[iecal][iptbin][isample]   -> Sumw2();
       fBremUnsplitEle[iecal][iptbin][isample] -> Sumw2();
       sigmaIEtaIEtaUnsplitEle[iecal][iptbin][isample] -> Sumw2();
       sigmaIPhiIPhiUnsplitEle[iecal][iptbin][isample] -> Sumw2();

       dPhiUnsplitEle[iecal][iptbin][isample]  -> 
	 Scale(dPhiUnsplitEle[iecal][iptbin][isample]->Integral());
       dEtaUnsplitEle[iecal][iptbin][isample]  -> 
	 Scale(dEtaUnsplitEle[iecal][iptbin][isample]->Integral());
       EoPUnsplitEle[iecal][iptbin][isample]   -> 
	 Scale(EoPUnsplitEle[iecal][iptbin][isample]->Integral());
       HoEUnsplitEle[iecal][iptbin][isample]   -> 
	 Scale(HoEUnsplitEle[iecal][iptbin][isample]->Integral());
       fBremUnsplitEle[iecal][iptbin][isample] -> 
	 Scale(fBremUnsplitEle[iecal][iptbin][isample]->Integral());
       sigmaIEtaIEtaUnsplitEle[iecal][iptbin][isample] -> 
	 Scale(sigmaIEtaIEtaUnsplitEle[iecal][iptbin][isample]->Integral());
       sigmaIPhiIPhiUnsplitEle[iecal][iptbin][isample] -> 
	 Scale(sigmaIPhiIPhiUnsplitEle[iecal][iptbin][isample]->Integral());
     }}}
 

 // cosmetics
 for (int iecal=0; iecal<2; iecal++) {
   for(int iptbin=0; iptbin<2; iptbin++) {
     dPhiUnsplitEle[iecal][iptbin][0]  -> SetFillColor(kPink);
     dEtaUnsplitEle[iecal][iptbin][0]  -> SetFillColor(kPink);
     EoPUnsplitEle[iecal][iptbin][0]   -> SetFillColor(kPink);
     HoEUnsplitEle[iecal][iptbin][0]   -> SetFillColor(kPink);
     fBremUnsplitEle[iecal][iptbin][0] -> SetFillColor(kPink);
     sigmaIEtaIEtaUnsplitEle[iecal][iptbin][0] -> SetFillColor(kPink);
     sigmaIPhiIPhiUnsplitEle[iecal][iptbin][0] -> SetFillColor(kPink);
   }}



 // plots
 TLegend a(0.7,0.8,0.85,0.92);
 a.AddEntry(dPhiUnsplitEle[0][1][0],"W candle","f");
 a.AddEntry(dPhiUnsplitEle[0][1][1],"Jet HLT", "l");
 a.Draw();
 a.SetFillColor(0);
 a.SetBorderSize(0.4);
 
 TCanvas *c11 = new TCanvas("c11", "",1);  
 dPhiUnsplitEle[0][1][0]->Draw("hist");
 dPhiUnsplitEle[0][1][1]->Draw("same");
 a.Draw();
 c11->SaveAs("dPhiUnsplitEle_EB.eps");
 c11->SaveAs("dPhiUnsplitEle_EB.root");

 TCanvas *c21 = new TCanvas("c21", "",1);  
 dEtaUnsplitEle[0][1][0]->Draw("hist");
 dEtaUnsplitEle[0][1][1]->Draw("same");
 a.Draw();
 c21->SaveAs("dEtaUnsplitEle_EB.eps");
 c21->SaveAs("dEtaUnsplitEle_EB.root");

 TCanvas *c31 = new TCanvas("c31", "",1);  
 EoPUnsplitEle[0][1][0]->Draw("hist");
 EoPUnsplitEle[0][1][1]->Draw("same");
 a.Draw();
 c31->SaveAs("EoPUnsplitEle_EB.eps");
 c31->SaveAs("EoPUnsplitEle_EB.root");

 TCanvas *c41 = new TCanvas("c41", "",1);  
 HoEUnsplitEle[0][1][0]->Draw("hist");
 HoEUnsplitEle[0][1][1]->Draw("same");
 a.Draw();
 c41->SaveAs("HoEUnsplitEle_EB.eps");
 c41->SaveAs("HoEUnsplitEle_EB.root");

 TCanvas *c51 = new TCanvas("c51", "",1);  
 fBremUnsplitEle[0][1][0]->Draw("hist");
 fBremUnsplitEle[0][1][1]->Draw("same");
 a.Draw();
 c51->SaveAs("fBremUnsplitEle_EB.eps");
 c51->SaveAs("fBremUnsplitEle_EB.root");

 TCanvas *c61 = new TCanvas("c61", "",1);  
 sigmaIEtaIEtaUnsplitEle[0][1][0]->Draw("hist");
 sigmaIEtaIEtaUnsplitEle[0][1][1]->Draw("same");
 a.Draw();
 c61->SaveAs("sigmaIEtaIEtaUnsplitEle_EB.eps");
 c61->SaveAs("sigmaIEtaIEtaUnsplitEle_EB.root");

 TCanvas *c71 = new TCanvas("c71", "",1);  
 sigmaIPhiIPhiUnsplitEle[0][1][0]->Draw("hist");
 sigmaIPhiIPhiUnsplitEle[0][1][1]->Draw("same");
 a.Draw();
 c71->SaveAs("sigmaIPhiIPhiUnsplitEle_EB.eps");
 c71->SaveAs("sigmaIPhiIPhiUnsplitEle_EB.root");

 TCanvas *c12 = new TCanvas("c12", "",1);  
 dPhiUnsplitEle[1][1][0]->Draw("hist");
 dPhiUnsplitEle[1][1][1]->Draw("same");
 a.Draw();
 c12->SaveAs("dPhiUnsplitEle_EE.eps");
 c12->SaveAs("dPhiUnsplitEle_EE.root");

 TCanvas *c22 = new TCanvas("c22", "",1);  
 dEtaUnsplitEle[1][1][0]->Draw("hist");
 dEtaUnsplitEle[1][1][1]->Draw("same");
 a.Draw();
 c22->SaveAs("dEtaUnsplitEle_EE.eps");
 c22->SaveAs("dEtaUnsplitEle_EE.root");

 TCanvas *c32 = new TCanvas("c32", "",1);  
 EoPUnsplitEle[1][1][0]->Draw("hist");
 EoPUnsplitEle[1][1][1]->Draw("same");
 a.Draw();
 c32->SaveAs("EoPUnsplitEle_EE.eps");
 c32->SaveAs("EoPUnsplitEle_EE.root");

 TCanvas *c42 = new TCanvas("c42", "",1);  
 HoEUnsplitEle[1][1][0]->Draw("hist");
 HoEUnsplitEle[1][1][1]->Draw("same");
 a.Draw();
 c42->SaveAs("HoEUnsplitEle_EE.eps");
 c42->SaveAs("HoEUnsplitEle_EE.root");

 TCanvas *c52 = new TCanvas("c52", "",1);  
 fBremUnsplitEle[1][1][0]->Draw("hist");
 fBremUnsplitEle[1][1][1]->Draw("same");
 a.Draw();
 c52->SaveAs("fBremUnsplitEle_EE.eps");
 c52->SaveAs("fBremUnsplitEle_EE.root");

 TCanvas *c62 = new TCanvas("c62", "",1);  
 sigmaIEtaIEtaUnsplitEle[1][1][0]->Draw("hist");
 sigmaIEtaIEtaUnsplitEle[1][1][1]->Draw("same");
 a.Draw();
 c62->SaveAs("sigmaIEtaIEtaUnsplitEle_EE.eps");
 c62->SaveAs("sigmaIEtaIEtaUnsplitEle_EE.root");

 TCanvas *c72 = new TCanvas("c72", "",1);  
 sigmaIPhiIPhiUnsplitEle[1][1][0]->Draw("hist");
 sigmaIPhiIPhiUnsplitEle[1][1][1]->Draw("same");
 a.Draw();
 c72->SaveAs("sigmaIPhiIPhiUnsplitEle_EE.eps");
 c72->SaveAs("sigmaIPhiIPhiUnsplitEle_EE.root");
}
