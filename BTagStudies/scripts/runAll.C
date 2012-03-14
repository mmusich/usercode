{
  
  //===================================================
  // Disclaimer:
  //
  // histogram of differences vs some third variable are always booked as:
  // valueA - valueB vs valueB ======> valueB is ALWAYS the reference!
  //===================================================


  // match the events in the ntuples
  
  gROOT->LoadMacro("./EventByEventNtuple.C+");
  EventByEventNtuple("../data/standardPFNtuple_ttbar_newgeom.root=newgeom","../data/standardPFNtuple_ttbar_orig.root=orig");
  
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

     rfio:/castor/cern.ch/user/e/emiglior/Alignment/bTag/CMSSW_42X/TTJets_Summer11-PU_S4/500k/standardPFNtuple_IDEAL_PF_NTPL.root
     rfio:/castor/cern.ch/user/e/emiglior/Alignment/bTag/CMSSW_42X/TTJets_Summer11-PU_S4/500k/standardPFNtuple_IDEAL_ZEROAPE_PF_NTPL.root
     rfio:/castor/cern.ch/user/e/emiglior/Alignment/bTag/CMSSW_42X/TTJets_Summer11-PU_S4/500k/standardPFNtuple_STARTUP2011_BS44_PF_NTPL.root
     rfio:/castor/cern.ch/user/e/emiglior/Alignment/bTag/CMSSW_42X/TTJets_Summer11-PU_S4/500k/standardPFNtuple_bpix500um_PF_NTPL.root
     rfio:/castor/cern.ch/user/e/emiglior/Alignment/bTag/CMSSW_42X/TTJets_Summer11-PU_S4/500k/standardPFNtuple_bpix500um_ZEROAPE_PF_NTPL.root 
     
  */

  // needed auxiliary class
  gROOT->LoadMacro("./JetInfo.cxx+");
 
  // match the ntuples (bool -> tree or flat ntuple,outfile name,maxEvents) 
  gROOT->LoadMacro("./MatchTheTree.C+");
  MatchTheTree(true,"MatrixOfMatches.root",-1);

  // graphics
  gROOT->LoadMacro("./tdrStyle.C+");
  setTDRStyle("gray");
  //setTDRStyle("logredblue");

  // this will work only for TTree and not for flat Ntuple
  gROOT->LoadMacro("./JetByJetComparisonHistos.cxx++g");
  gROOT->LoadMacro("./JetByJetComparison.C++g");
  JetByJetComparison t("JetByJetComparisonTree_newgeomVsorig.root");
  t.Loop(); 
  
  // TString processline = Form(".! mv *.png figures_jetbyjet") ;
  // gROOT->ProcessLine(processline.Data());
  // gSystem->Sleep(100);  

}
