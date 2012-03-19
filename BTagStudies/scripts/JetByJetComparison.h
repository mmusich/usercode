//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 21 08:03:11 2012 by ROOT version 5.28/00
// from TTree JetByJetComparisonTree/tree file1 vs file2
// found on file: JetByJetComparison_origVsnewgeom.root
//////////////////////////////////////////////////////////

#ifndef JetByJetComparison_h
#define JetByJetComparison_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TObjArray.h>
#include "JetInfo.h"

#include <iostream>

using namespace std; 

class JetByJetComparison {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   TString         CompNames[2];

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           evt;
   Float_t         pthat;
   Float_t         mcweight;
   Char_t          isBGluonSplitting;
   Char_t          isCGluonSplitting;

   JetInfo        *JetInfoA;
   JetInfo        *JetInfoB;

   TBranch        *b_Trun;   //!
   TBranch        *b_Tlumi;   //!
   TBranch        *b_Tevent;   //!
   TBranch        *b_Tpthat;   //!
   TBranch        *b_Tmcweight;   //!
   TBranch        *b_TisBGluonSplitting;   //!
   TBranch        *b_TisCGluonSplitting;   //!
 
   JetByJetComparison(TString filename, TTree *tree=0);
   virtual ~JetByJetComparison();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef JetByJetComparison_cxx

JetByJetComparison::JetByJetComparison(TString filename,TTree *tree)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    // TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("JetByJetComparison_origVsnewgeom.root");
    // TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
    TFile *f = TFile::Open(filename)
    if (!f) {
      //    f = new TFile("JetByJetComparison_origVsnewgeom.root");
      f = new TFile(filename);
    }
    tree = (TTree*)gDirectory->Get("JetByJetComparisonTree");
    
  }
  Init(tree);
  filename.ReplaceAll(".root","");
  TObjArray *allInfos = filename.Tokenize("_");
  TString filestoCompare_s = allInfos->At(1)->GetName();
  TObjArray *filestoCompare = (filestoCompare_s).Tokenize("Vs");
  
  CompNames[0] = filestoCompare->At(0)->GetName();
  CompNames[1] = filestoCompare->At(1)->GetName();
}

JetByJetComparison::~JetByJetComparison()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JetByJetComparison::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t JetByJetComparison::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void JetByJetComparison::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   //  fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run",   &run,   &b_Trun);
   fChain->SetBranchAddress("lumi",  &lumi,  &b_Tlumi);
   fChain->SetBranchAddress("evt",   &evt,   &b_Tevent);
   fChain->SetBranchAddress("pthat", &pthat, &b_Tpthat);
   fChain->SetBranchAddress("mcweight", &mcweight, &b_Tmcweight);
   fChain->SetBranchAddress("isBGluonSplitting", &isBGluonSplitting, &b_TisBGluonSplitting);
   fChain->SetBranchAddress("isCGluonSplitting", &isCGluonSplitting, &b_TisCGluonSplitting);

   JetInfoA=new JetInfo();
   JetInfoB=new JetInfo();   
   fChain->SetBranchAddress("JetInfoA", &JetInfoA);
   fChain->SetBranchAddress("JetInfoB", &JetInfoB);

   Notify();
}

Bool_t JetByJetComparison::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void JetByJetComparison::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t JetByJetComparison::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef JetByJetComparison_cxx
