//////////////////////////////////////////////////////////
//
// Simple script to run ROOT analysis job:
//
// Usage:
//
// root -b
// .x Execute("yourFavoriteFolder")
//
// will copy your plots to ./yourFavoriteFolder
//
// Original Author: M. Musich INFN Torino
//
//////////////////////////////////////////////////////////

void Execute(TString outdirname){
  
  gSystem->mkdir(outdirname);
  gROOT->LoadMacro("./ZBB.C");
  ZBB("mu");
  ZBB("ele"); 
  gSystem->Sleep(500); 
  gROOT->LoadMacro("./FastSuperImposeHistos.C++");
  FastSuperImposeHistos("plots_datasetMu.root=Mu dataset,plots_datasetEle.root=Ele dataset",2);
  TString processline = Form(".! mv *.eps %s",outdirname.Data()) ;
  gROOT->ProcessLine(processline.Data());

}
