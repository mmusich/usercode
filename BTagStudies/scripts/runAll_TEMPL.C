{
  
  gROOT->LoadMacro("./EventByEventNtuple.C+");
  EventByEventNtuple("FILE1TEMPLATE=NAME1TEMPLATE","FILE2TEMPLATE=NAME2TEMPLATE");
  
  gROOT->LoadMacro("./JetInfo.cxx+");
 
  gROOT->LoadMacro("./MatchTheTree.C+");
  MatchTheTree(true,"MatrixOfMatches.root",-1);

  // this will work only for TTree and not for flat Ntuple
  //gROOT->LoadMacro("./JetByJetComparisonHistos.cxx++g");
  //gROOT->LoadMacro("./JetByJetComparison.C++g");
  //JetByJetComparison t("JetByJetComparisonTree_origVsnewgeom.root");
  //t.Loop(); 
  
  // TString processline = Form(".! mv *.png figures_jetbyjet") ;
  // gROOT->ProcessLine(processline.Data());
  // gSystem->Sleep(100);  

}
