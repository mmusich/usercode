{
  
  gROOT->LoadMacro("./EventByEventNtuple.C+");
  EventByEventNtuple("../data/standardPFNtuple_ttbar_orig.root=orig","../data/standardPFNtuple_ttbar_newgeom.root=newgeom");
  
  /* list of available files for analysis
     
     => DATA:
     
     rfio:/castor/cern.ch/user/m/musich/BtaggingNtuples/DATA/standardPFNtuple_GR10v4.root
     rfio:/castor/cern.ch/user/m/musich/BtaggingNtuples/DATA/standardPFNtuple_GR10v4_TkAlAPEzero.root
     rfio:/castor/cern.ch/user/m/musich/BtaggingNtuples/DATA/standardPFNtuple_ideal.root
     rfio:/castor/cern.ch/user/m/musich/BtaggingNtuples/DATA/standardPFNtuple_ideal_TkAlAPEzero.root
     rfio:/castor/cern.ch/user/m/musich/BtaggingNtuples/DATA/standardPFNtuple_newFlat.root
     rfio:/castor/cern.ch/user/m/musich/BtaggingNtuples/DATA/standardPFNtuple_newFlat_TkAlAPEzero.root
     rfio:/castor/cern.ch/user/m/musich/BtaggingNtuples/DATA/standardPFNtuple_newKB.root
     rfio:/castor/cern.ch/user/m/musich/BtaggingNtuples/DATA/standardPFNtuple_newKB_TkAlAPEzero.root
     
     => MC:
     
     rfio:/castor/cern.ch/user/m/musich/BtaggingNtuples/MC/standardPFNtuple_ttbar_orig_1.root
     rfio:/castor/cern.ch/user/m/musich/BtaggingNtuples/MC/standardPFNtuple_ttbar_newgeom_1.root
     
  */

  gROOT->LoadMacro("./JetInfo.cxx+");
 
  gROOT->LoadMacro("./MatchTheTree.C+");
  MatchTheTree(true,"MatrixOfMatches.root",-1);

  // this will work only for TTree and not for flat Ntuple

  gROOT->LoadMacro("./tdrStyle.C+");
  //  setTDRStyle("blues");
  setTDRStyle("logredblue");
  gROOT->LoadMacro("./JetByJetComparisonHistos.cxx++g");
  gROOT->LoadMacro("./JetByJetComparison.C++g");
  JetByJetComparison t("JetByJetComparisonTree_origVsnewgeom.root");
  t.Loop(); 
  
  // TString processline = Form(".! mv *.png figures_jetbyjet") ;
  // gROOT->ProcessLine(processline.Data());
  // gSystem->Sleep(100);  


}
