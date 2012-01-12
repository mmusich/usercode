{
  gROOT->LoadMacro("./EventByEventNtuple.C+");
  EventByEventNtuple("../data/RunRangeIdeal-165364.root=ideal","../data/RunRangeNEWKB-165364.root=KinksAndBows");
  gROOT->LoadMacro("./MatchTheTree.C+");
  MatchTheTree("MatrixOfMatches.root");
}
