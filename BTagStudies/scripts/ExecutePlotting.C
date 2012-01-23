{
  
  gSystem->mkdir("figures_jetbyjet");
  gROOT->LoadMacro("./JetByJetComparison.C++g");
  JetByJetComparison t("JetByJetComparison_origVsnewgeom.root");
  t.Loop(); 
  
  TString processline = Form(".! mv *.png figures_jetbyjet") ;
  gROOT->ProcessLine(processline.Data());
  gSystem->Sleep(100);
}
