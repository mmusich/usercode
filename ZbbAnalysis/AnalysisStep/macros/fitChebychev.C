void fitChebychev(){
  gROOT->ProcessLine(".x ./rootlogon.C");
  
  //TFile f("analyzePAT_MC_ZlJets_All.root");
  //TFile f("analyzePAT_DATA2011.root");
  //TFile f("analyzePAT_DATA2011.root");
  TFile f("analyzePAT_MC_TTJets_All.root");
  f.cd("finaldistros_ssvhem/JetBTag");
  Z_mass=(TH1F*)JetBTag_MllNOMASSCUT->Clone();
  Z_mass->Rebin(2);
  
  // double hmin = Z_mass->GetXaxis()->GetXmin();
  // double hmax = Z_mass->GetXaxis()->GetXmax();
  
  float hmin = 60.0;
  float hmax = 120.0;
  
  // Declare observable x
  RooRealVar x("x","x",hmin,hmax) ;
  RooDataHist dh("dh","dh",x,Import(*Z_mass)) ;
  
  RooPlot* frame = x.frame(Title("Z mass")) ;
  dh.plotOn(frame,MarkerColor(2),MarkerSize(0.9),MarkerStyle(21));  //this will show histogram data points on canvas 
  dh.statOn(frame);  //this will display hist stat on canvas
  
  //RooRealVar mean("mean","mean",91.0, 70.0, 120.0);
  //RooRealVar width("width","width",10.0, 0.0, 120.0);
  //RooRealVar sigma("sigma","sigma",5.0, 0.0, 120.0);
  //RooGaussian gauss("gauss","gauss",x,mean,sigma);
  //RooBreitWigner gauss("gauss","gauss",x,mean,sigma);
  //RooVoigtian gauss("gauss","gauss",x,mean,width,sigma);
  
  //polinomio bkg tt
  RooRealVar c0("c0", "", -0.5,-1., 1.);
  RooRealVar c1("c1", "", -0.3,-1., 1.);
  RooRealVar c2("c2", "", -0.3,-1., 1.);
  RooChebychev pdf_bkg("pdf_bkg", "", x, RooArgList(c0, c1));

  RooFitResult* filters = pdf_bkg.fitTo(dh,"qr");
  pdf_bkg.plotOn(frame,LineColor(4));//this will show fit overlay on canvas 
  pdf_bkg.paramOn(frame); //this will display the fit parameters on canvas
  //filters->Print("v");
  
  // Draw all frames on a canvas
  TCanvas* c = new TCanvas("ZmassHisto","ZmassHisto",800,400) ;
  c->cd() ; gPad->SetLeftMargin(0.15);
  
  frame->GetXaxis()->SetTitle("Z mass (in GeV/c^{2})");  frame->GetXaxis()->SetTitleOffset(1.2);
  float binsize = Z_mass->GetBinWidth(1); char Bsize[50]; 
  //sprintf(Bsize,"Events per %2.2f",binsize);
  // frame->GetYaxis()->SetTitle(Bsize);  
  //frame->GetYaxis()->SetTitleOffset(1.2);
  frame->Draw() ;
  c->Print("myZmaa.png");  
  
}


