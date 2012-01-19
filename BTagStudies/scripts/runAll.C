{
  
  gROOT->LoadMacro("./EventByEventNtuple.C+");
  //EventByEventNtuple("rfio:/castor/cern.ch/user/a/atravers/StudiesJetByJet/standardPFNtuple_ideal.root=Ideal","rfio:/castor/cern.ch/user/a/atravers/StudiesJetByJet/standardPFNtuple_newKB.root=KinksAndBows");
  EventByEventNtuple("rfio:/castor/cern.ch/user/m/musich/BtaggingNtuples/MC/standardPFNtuple_ttbar_orig_1.root=orig","rfio:/castor/cern.ch/user/m/musich/BtaggingNtuples/MC/standardPFNtuple_ttbar_newgeom_1.root=newgeom");

  /* list of available files

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

  gROOT->LoadMacro("./MatchTheTree.C+");
  MatchTheTree("MatrixOfMatches.root");
  
}
