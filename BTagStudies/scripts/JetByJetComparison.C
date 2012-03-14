#define JetByJetComparison_cxx

#include "JetByJetComparisonHistos.h"
#include "JetByJetComparison.h"

#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>


using namespace std; 

bool isUDSG(int mctruth) {
  return mctruth == 21 || mctruth < 4;
}

bool isB(int mctruth) {
  return mctruth == 5;
}

void JetByJetComparison::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L JetByJetComparison.C
//      Root > JetByJetComparison t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


  TString file_out_name = "JetByJetComparisonPlots_"+CompNames[0]+"_vs_"+CompNames[1]+".root";
  TFile *file_out=new TFile(file_out_name,"recreate");  
  file_out->cd();

  JetByJetComparisonHistos jetbyjethistos_notdefnotdef_all("NotDefault_NotDefault_all",file_out);
  JetByJetComparisonHistos jetbyjethistos_notdefnotdef_b("NotDefault_NotDefault_b",file_out);
  JetByJetComparisonHistos jetbyjethistos_notdefnotdef_udsg("NotDefault_NotDefault_udsg",file_out);
//  JetByJetComparisonHistos jetbyjethistos_notdefdef("NotDefault_Default",file_out);
//  JetByJetComparisonHistos jetbyjethistos_defnotdef("Default_NotDefault",file_out);

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if ( jentry%10000 == 0 )  cout << "JetByJetComprison Processing Event: " << jentry << endl;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    jetbyjethistos_notdefnotdef_all.fillAllHistos(*JetInfoA,*JetInfoB,file_out);
    if (isB(JetInfoA->MCTrueFlavor) && isB(JetInfoB->MCTrueFlavor)) {
      jetbyjethistos_notdefnotdef_b.fillAllHistos(*JetInfoA,*JetInfoB,file_out);
    } else if (isUDSG(JetInfoA->MCTrueFlavor) && isUDSG(JetInfoB->MCTrueFlavor) ) {
      jetbyjethistos_notdefnotdef_udsg.fillAllHistos(*JetInfoA,*JetInfoB,file_out);
    }
    //   jetbyjethistos_notdefdef.fillAllHistos(*JetInfoA,*JetInfoB,file_out);
    //  jetbyjethistos_defnotdef.fillAllHistos(*JetInfoA,*JetInfoB,file_out);
     if (Cut(ientry) < 0) continue;
  }
  
   jetbyjethistos_notdefnotdef_all.drawNice2dHistos(file_out);
   jetbyjethistos_notdefnotdef_b.drawNice2dHistos(file_out);
   jetbyjethistos_notdefnotdef_udsg.drawNice2dHistos(file_out);
   // jetbyjethistos_notdefdef.drawNice2dHistos(file_out);
   // jetbyjethistos_defnotdef.drawNice2dHistos(file_out);


  // end of job stuff
   file_out->cd();
   file_out->Write();
   file_out->Close();

}
