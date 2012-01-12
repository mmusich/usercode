void FilterIntoRunRanges() {


  TChain *chain = new TChain("t");
  chain->Add("standardPFNtuple_newKB.root");
  // chain->Add("/storage/5/jyothsna/BTagCommissioning2010_Sep18/standardPFNtuple_JetMETTau_Run2010A_PromptReco_v4_101_1_vXO.root");
  chain->GetEntry(0); 
  
  //Create a new file + a clone of old tree in new file
  TFile *newfile = new TFile("RunRangeNEWKB-165364.root","recreate");

   /*
     (A) 138562 - 139790
     (B) 139965 - 140399
     (C) 142928 - 144011
     (D) 146428 - 146807
   */

  Int_t nentries = Int_t(chain->GetEntries());
  
  std::cout << "+++++ No. of entries in the chain: " << nentries << std::endl;

  TTree *newtree = t->CopyTree("runNumber == 165364");

  newtree->Print();
  newtree->AutoSave();

  delete newfile;
  delete chain;
}
