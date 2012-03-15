void ExecutePlotting(TString foldername){
  
  /* available files
     
     JetByJetComparisonTree_DeltaZ500noAPEVsIdealnoAPE.root
     JetByJetComparisonTree_DeltaZ500VsDeltaZ500noAPE.root
     JetByJetComparisonTree_DeltaZ500VsIdeal.root
     JetByJetComparisonTree_IdealVsIdealnoAPE.root
     JetByJetComparisonTree_Startup2011VsIdeal.root
     JetByJetComparisonTree_DeltaZ40VsDeltaZ40noAPE.root                                          
     JetByJetComparisonTree_DeltaZ40VsIdeal.root     
  */

  std::string foldername_s = foldername;
  
  gSystem->mkdir(foldername);
  gROOT->LoadMacro("./JetInfo.cxx+");

  gROOT->LoadMacro("./tdrStyle.C+");
  setTDRStyle("gray");
  //setTDRStyle("logredblue");
  gROOT->LoadMacro("./JetByJetComparisonHistos.cxx++g");
  gROOT->LoadMacro("./JetByJetComparison.C++g");
  JetByJetComparison t("../data/JetByJetComparisonTree_"+foldername+".root");
  t.Loop(); 
  
  gROOT->LoadMacro("./diow.C+");
  diow(".","index.html");

  TString processline = Form(".! mv *.png index.html %s",foldername_s) ;
  gROOT->ProcessLine(processline.Data());
  gSystem->Sleep(100);
  processline.Clear();

  processline = Form(".! cp diow.C %s",foldername_s) ;
  gROOT->ProcessLine(processline.Data());
 
}
