//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 11 12:37:06 2013 by ROOT version 5.34/03
// from TTree TkOffVal/TkOffVal
// found on file: AlignmentValidation_WeeklyValidation_shiftPlots.root
//////////////////////////////////////////////////////////

#ifndef OfflineValidationTreeAnalysis_h
#define OfflineValidationTreeAnalysis_h

#include <TROOT.h>
#include <TProfile.h>
#include <TChain.h>
#include <TFile.h>
#include <TColor.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TObjArray.h>
#include <TPaveText.h>
#include <iostream>
#include <TCanvas.h>
//#include <TIter.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class OfflineValidationTreeAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TString         fChainName;
   TFile          *fOutputFile;
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //TkOffTreeVariables *TkOffTreeVariables;
   Float_t         meanLocalX;
   Float_t         meanNormLocalX;
   Float_t         meanX;
   Float_t         meanNormX;
   Float_t         meanY;
   Float_t         meanNormY;
   Float_t         medianX;
   Float_t         medianY;
   Float_t         chi2PerDofX;
   Float_t         chi2PerDofY;
   Float_t         rmsLocalX;
   Float_t         rmsNormLocalX;
   Float_t         rmsX;
   Float_t         rmsNormX;
   Float_t         rmsY;
   Float_t         rmsNormY;
   Float_t         sigmaX;
   Float_t         sigmaNormX;
   Float_t         fitMeanX;
   Float_t         fitSigmaX;
   Float_t         fitMeanNormX;
   Float_t         fitSigmaNormX;
   Float_t         fitMeanY;
   Float_t         fitSigmaY;
   Float_t         fitMeanNormY;
   Float_t         fitSigmaNormY;
   Float_t         posR;
   Float_t         posPhi;
   Float_t         posEta;
   Float_t         posX;
   Float_t         posY;
   Float_t         posZ;
   Float_t         numberOfUnderflows;
   Float_t         numberOfOverflows;
   Float_t         numberOfOutliers;
   Float_t         rDirection;
   Float_t         phiDirection;
   Float_t         zDirection;
   Float_t         rOrZDirection;
   UInt_t          entries;
   UInt_t          moduleId;
   UInt_t          subDetId;
   UInt_t          layer;
   UInt_t          side;
   UInt_t          half;
   UInt_t          rod;
   UInt_t          ring;
   UInt_t          petal;
   UInt_t          blade;
   UInt_t          panel;
   UInt_t          outerInner;
   UInt_t          module;
   Bool_t          isDoubleSide;
   Bool_t          isStereo;
   string          histNameLocalX;
   string          histNameNormLocalX;
   string          histNameLocalY;
   string          histNameX;
   string          histNameNormX;
   string          histNameY;
   string          histNameNormY;
   Float_t         meanResXvsX;
   Float_t         meanResXvsY;
   Float_t         meanResYvsX;
   Float_t         meanResYvsY;
   Float_t         rmsResXvsX;
   Float_t         rmsResXvsY;
   Float_t         rmsResYvsX;
   Float_t         rmsResYvsY;
   string          profileNameResXvsX;
   string          profileNameResXvsY;
   string          profileNameResYvsX;
   string          profileNameResYvsY;

   // List of branches
   TBranch        *b_TkOffTreeVariables_meanLocalX;   //!
   TBranch        *b_TkOffTreeVariables_meanNormLocalX;   //!
   TBranch        *b_TkOffTreeVariables_meanX;   //!
   TBranch        *b_TkOffTreeVariables_meanNormX;   //!
   TBranch        *b_TkOffTreeVariables_meanY;   //!
   TBranch        *b_TkOffTreeVariables_meanNormY;   //!
   TBranch        *b_TkOffTreeVariables_medianX;   //!
   TBranch        *b_TkOffTreeVariables_medianY;   //!
   TBranch        *b_TkOffTreeVariables_chi2PerDofX;   //!
   TBranch        *b_TkOffTreeVariables_chi2PerDofY;   //!
   TBranch        *b_TkOffTreeVariables_rmsLocalX;   //!
   TBranch        *b_TkOffTreeVariables_rmsNormLocalX;   //!
   TBranch        *b_TkOffTreeVariables_rmsX;   //!
   TBranch        *b_TkOffTreeVariables_rmsNormX;   //!
   TBranch        *b_TkOffTreeVariables_rmsY;   //!
   TBranch        *b_TkOffTreeVariables_rmsNormY;   //!
   TBranch        *b_TkOffTreeVariables_sigmaX;   //!
   TBranch        *b_TkOffTreeVariables_sigmaNormX;   //!
   TBranch        *b_TkOffTreeVariables_fitMeanX;   //!
   TBranch        *b_TkOffTreeVariables_fitSigmaX;   //!
   TBranch        *b_TkOffTreeVariables_fitMeanNormX;   //!
   TBranch        *b_TkOffTreeVariables_fitSigmaNormX;   //!
   TBranch        *b_TkOffTreeVariables_fitMeanY;   //!
   TBranch        *b_TkOffTreeVariables_fitSigmaY;   //!
   TBranch        *b_TkOffTreeVariables_fitMeanNormY;   //!
   TBranch        *b_TkOffTreeVariables_fitSigmaNormY;   //!
   TBranch        *b_TkOffTreeVariables_posR;   //!
   TBranch        *b_TkOffTreeVariables_posPhi;   //!
   TBranch        *b_TkOffTreeVariables_posEta;   //!
   TBranch        *b_TkOffTreeVariables_posX;   //!
   TBranch        *b_TkOffTreeVariables_posY;   //!
   TBranch        *b_TkOffTreeVariables_posZ;   //!
   TBranch        *b_TkOffTreeVariables_numberOfUnderflows;   //!
   TBranch        *b_TkOffTreeVariables_numberOfOverflows;   //!
   TBranch        *b_TkOffTreeVariables_numberOfOutliers;   //!
   TBranch        *b_TkOffTreeVariables_rDirection;   //!
   TBranch        *b_TkOffTreeVariables_phiDirection;   //!
   TBranch        *b_TkOffTreeVariables_zDirection;   //!
   TBranch        *b_TkOffTreeVariables_rOrZDirection;   //!
   TBranch        *b_TkOffTreeVariables_entries;   //!
   TBranch        *b_TkOffTreeVariables_moduleId;   //!
   TBranch        *b_TkOffTreeVariables_subDetId;   //!
   TBranch        *b_TkOffTreeVariables_layer;   //!
   TBranch        *b_TkOffTreeVariables_side;   //!
   TBranch        *b_TkOffTreeVariables_half;   //!
   TBranch        *b_TkOffTreeVariables_rod;   //!
   TBranch        *b_TkOffTreeVariables_ring;   //!
   TBranch        *b_TkOffTreeVariables_petal;   //!
   TBranch        *b_TkOffTreeVariables_blade;   //!
   TBranch        *b_TkOffTreeVariables_panel;   //!
   TBranch        *b_TkOffTreeVariables_outerInner;   //!
   TBranch        *b_TkOffTreeVariables_module;   //!
   TBranch        *b_TkOffTreeVariables_isDoubleSide;   //!
   TBranch        *b_TkOffTreeVariables_isStereo;   //!
   TBranch        *b_TkOffTreeVariables_histNameLocalX;   //!
   TBranch        *b_TkOffTreeVariables_histNameNormLocalX;   //!
   TBranch        *b_TkOffTreeVariables_histNameLocalY;   //!
   TBranch        *b_TkOffTreeVariables_histNameX;   //!
   TBranch        *b_TkOffTreeVariables_histNameNormX;   //!
   TBranch        *b_TkOffTreeVariables_histNameY;   //!
   TBranch        *b_TkOffTreeVariables_histNameNormY;   //!
   TBranch        *b_TkOffTreeVariables_meanResXvsX;   //!
   TBranch        *b_TkOffTreeVariables_meanResXvsY;   //!
   TBranch        *b_TkOffTreeVariables_meanResYvsX;   //!
   TBranch        *b_TkOffTreeVariables_meanResYvsY;   //!
   TBranch        *b_TkOffTreeVariables_rmsResXvsX;   //!
   TBranch        *b_TkOffTreeVariables_rmsResXvsY;   //!
   TBranch        *b_TkOffTreeVariables_rmsResYvsX;   //!
   TBranch        *b_TkOffTreeVariables_rmsResYvsY;   //!
   TBranch        *b_TkOffTreeVariables_profileNameResXvsX;   //!
   TBranch        *b_TkOffTreeVariables_profileNameResXvsY;   //!
   TBranch        *b_TkOffTreeVariables_profileNameResYvsX;   //!
   TBranch        *b_TkOffTreeVariables_profileNameResYvsY;   //!

   OfflineValidationTreeAnalysis(TString filename="",TString theName="");
   virtual ~OfflineValidationTreeAnalysis();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual TPaveText* addNiceTLegend(TObjArray *histos);
   virtual void     makeNiceProfileStyleAndColor(TProfile *hist,int color);
   virtual void     makeNicePlotStyleAndColor(TH1F *hist,int color);
   virtual void     makeNicePlotStyle(TH1F *hist);
   virtual void     makeNice2DPlotStyle(TH2F *hist);
   virtual void     makeNiceCanv(TCanvas *c);
   virtual void     setStyle(TString palettename);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef OfflineValidationTreeAnalysis_cxx
OfflineValidationTreeAnalysis::OfflineValidationTreeAnalysis(TString filename,TString theName) : fChain(0) 
{
 
  TTree *tree=0;
  if (tree == 0) {  
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
    if (!f || !f->IsOpen()) {
   
      f = new TFile(filename);
    }
    TDirectory * dir = (TDirectory*)f->Get(filename+":/TrackerOfflineValidationStandalone");
    dir->GetObject("TkOffVal",tree);
  }
  Init(tree);
  theName.ReplaceAll(" ","_");
  fChainName=theName;
  fOutputFile = TFile::Open(Form("TrackerOfflineValidationStandalone_%s.root",fChainName.Data()),"RECREATE");
  fOutputFile->mkdir("Pulls");
  fOutputFile->mkdir("rmsDMR");
}

OfflineValidationTreeAnalysis::~OfflineValidationTreeAnalysis()
{
  fOutputFile->Close();
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

/*--------------------------------------------------------------------*/
Int_t OfflineValidationTreeAnalysis::GetEntry(Long64_t entry)
/*--------------------------------------------------------------------*/
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

/*--------------------------------------------------------------------*/
Long64_t OfflineValidationTreeAnalysis::LoadTree(Long64_t entry)
/*--------------------------------------------------------------------*/
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

/*--------------------------------------------------------------------*/
void OfflineValidationTreeAnalysis::Init(TTree *tree)
/*--------------------------------------------------------------------*/
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
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("meanLocalX", &meanLocalX, &b_TkOffTreeVariables_meanLocalX);
   fChain->SetBranchAddress("meanNormLocalX", &meanNormLocalX, &b_TkOffTreeVariables_meanNormLocalX);
   fChain->SetBranchAddress("meanX", &meanX, &b_TkOffTreeVariables_meanX);
   fChain->SetBranchAddress("meanNormX", &meanNormX, &b_TkOffTreeVariables_meanNormX);
   fChain->SetBranchAddress("meanY", &meanY, &b_TkOffTreeVariables_meanY);
   fChain->SetBranchAddress("meanNormY", &meanNormY, &b_TkOffTreeVariables_meanNormY);
   fChain->SetBranchAddress("medianX", &medianX, &b_TkOffTreeVariables_medianX);
   fChain->SetBranchAddress("medianY", &medianY, &b_TkOffTreeVariables_medianY);
   fChain->SetBranchAddress("chi2PerDofX", &chi2PerDofX, &b_TkOffTreeVariables_chi2PerDofX);
   fChain->SetBranchAddress("chi2PerDofY", &chi2PerDofY, &b_TkOffTreeVariables_chi2PerDofY);
   fChain->SetBranchAddress("rmsLocalX", &rmsLocalX, &b_TkOffTreeVariables_rmsLocalX);
   fChain->SetBranchAddress("rmsNormLocalX", &rmsNormLocalX, &b_TkOffTreeVariables_rmsNormLocalX);
   fChain->SetBranchAddress("rmsX", &rmsX, &b_TkOffTreeVariables_rmsX);
   fChain->SetBranchAddress("rmsNormX", &rmsNormX, &b_TkOffTreeVariables_rmsNormX);
   fChain->SetBranchAddress("rmsY", &rmsY, &b_TkOffTreeVariables_rmsY);
   fChain->SetBranchAddress("rmsNormY", &rmsNormY, &b_TkOffTreeVariables_rmsNormY);
   fChain->SetBranchAddress("sigmaX", &sigmaX, &b_TkOffTreeVariables_sigmaX);
   fChain->SetBranchAddress("sigmaNormX", &sigmaNormX, &b_TkOffTreeVariables_sigmaNormX);
   fChain->SetBranchAddress("fitMeanX", &fitMeanX, &b_TkOffTreeVariables_fitMeanX);
   fChain->SetBranchAddress("fitSigmaX", &fitSigmaX, &b_TkOffTreeVariables_fitSigmaX);
   fChain->SetBranchAddress("fitMeanNormX", &fitMeanNormX, &b_TkOffTreeVariables_fitMeanNormX);
   fChain->SetBranchAddress("fitSigmaNormX", &fitSigmaNormX, &b_TkOffTreeVariables_fitSigmaNormX);
   fChain->SetBranchAddress("fitMeanY", &fitMeanY, &b_TkOffTreeVariables_fitMeanY);
   fChain->SetBranchAddress("fitSigmaY", &fitSigmaY, &b_TkOffTreeVariables_fitSigmaY);
   fChain->SetBranchAddress("fitMeanNormY", &fitMeanNormY, &b_TkOffTreeVariables_fitMeanNormY);
   fChain->SetBranchAddress("fitSigmaNormY", &fitSigmaNormY, &b_TkOffTreeVariables_fitSigmaNormY);
   fChain->SetBranchAddress("posR", &posR, &b_TkOffTreeVariables_posR);
   fChain->SetBranchAddress("posPhi", &posPhi, &b_TkOffTreeVariables_posPhi);
   fChain->SetBranchAddress("posEta", &posEta, &b_TkOffTreeVariables_posEta);
   fChain->SetBranchAddress("posX", &posX, &b_TkOffTreeVariables_posX);
   fChain->SetBranchAddress("posY", &posY, &b_TkOffTreeVariables_posY);
   fChain->SetBranchAddress("posZ", &posZ, &b_TkOffTreeVariables_posZ);
   fChain->SetBranchAddress("numberOfUnderflows", &numberOfUnderflows, &b_TkOffTreeVariables_numberOfUnderflows);
   fChain->SetBranchAddress("numberOfOverflows", &numberOfOverflows, &b_TkOffTreeVariables_numberOfOverflows);
   fChain->SetBranchAddress("numberOfOutliers", &numberOfOutliers, &b_TkOffTreeVariables_numberOfOutliers);
   fChain->SetBranchAddress("rDirection", &rDirection, &b_TkOffTreeVariables_rDirection);
   fChain->SetBranchAddress("phiDirection", &phiDirection, &b_TkOffTreeVariables_phiDirection);
   fChain->SetBranchAddress("zDirection", &zDirection, &b_TkOffTreeVariables_zDirection);
   fChain->SetBranchAddress("rOrZDirection", &rOrZDirection, &b_TkOffTreeVariables_rOrZDirection);
   fChain->SetBranchAddress("entries", &entries, &b_TkOffTreeVariables_entries);
   fChain->SetBranchAddress("moduleId", &moduleId, &b_TkOffTreeVariables_moduleId);
   fChain->SetBranchAddress("subDetId", &subDetId, &b_TkOffTreeVariables_subDetId);
   fChain->SetBranchAddress("layer", &layer, &b_TkOffTreeVariables_layer);
   fChain->SetBranchAddress("side", &side, &b_TkOffTreeVariables_side);
   fChain->SetBranchAddress("half", &half, &b_TkOffTreeVariables_half);
   fChain->SetBranchAddress("rod", &rod, &b_TkOffTreeVariables_rod);
   fChain->SetBranchAddress("ring", &ring, &b_TkOffTreeVariables_ring);
   fChain->SetBranchAddress("petal", &petal, &b_TkOffTreeVariables_petal);
   fChain->SetBranchAddress("blade", &blade, &b_TkOffTreeVariables_blade);
   fChain->SetBranchAddress("panel", &panel, &b_TkOffTreeVariables_panel);
   fChain->SetBranchAddress("outerInner", &outerInner, &b_TkOffTreeVariables_outerInner);
   fChain->SetBranchAddress("module", &module, &b_TkOffTreeVariables_module);
   fChain->SetBranchAddress("isDoubleSide", &isDoubleSide, &b_TkOffTreeVariables_isDoubleSide);
   fChain->SetBranchAddress("isStereo", &isStereo, &b_TkOffTreeVariables_isStereo);
   fChain->SetBranchAddress("histNameLocalX", &histNameLocalX, &b_TkOffTreeVariables_histNameLocalX);
   fChain->SetBranchAddress("histNameNormLocalX", &histNameNormLocalX, &b_TkOffTreeVariables_histNameNormLocalX);
   fChain->SetBranchAddress("histNameLocalY", &histNameLocalY, &b_TkOffTreeVariables_histNameLocalY);
   fChain->SetBranchAddress("histNameX", &histNameX, &b_TkOffTreeVariables_histNameX);
   fChain->SetBranchAddress("histNameNormX", &histNameNormX, &b_TkOffTreeVariables_histNameNormX);
   fChain->SetBranchAddress("histNameY", &histNameY, &b_TkOffTreeVariables_histNameY);
   fChain->SetBranchAddress("histNameNormY", &histNameNormY, &b_TkOffTreeVariables_histNameNormY);
   fChain->SetBranchAddress("meanResXvsX", &meanResXvsX, &b_TkOffTreeVariables_meanResXvsX);
   fChain->SetBranchAddress("meanResXvsY", &meanResXvsY, &b_TkOffTreeVariables_meanResXvsY);
   fChain->SetBranchAddress("meanResYvsX", &meanResYvsX, &b_TkOffTreeVariables_meanResYvsX);
   fChain->SetBranchAddress("meanResYvsY", &meanResYvsY, &b_TkOffTreeVariables_meanResYvsY);
   fChain->SetBranchAddress("rmsResXvsX", &rmsResXvsX, &b_TkOffTreeVariables_rmsResXvsX);
   fChain->SetBranchAddress("rmsResXvsY", &rmsResXvsY, &b_TkOffTreeVariables_rmsResXvsY);
   fChain->SetBranchAddress("rmsResYvsX", &rmsResYvsX, &b_TkOffTreeVariables_rmsResYvsX);
   fChain->SetBranchAddress("rmsResYvsY", &rmsResYvsY, &b_TkOffTreeVariables_rmsResYvsY);
   fChain->SetBranchAddress("profileNameResXvsX", &profileNameResXvsX, &b_TkOffTreeVariables_profileNameResXvsX);
   fChain->SetBranchAddress("profileNameResXvsY", &profileNameResXvsY, &b_TkOffTreeVariables_profileNameResXvsY);
   fChain->SetBranchAddress("profileNameResYvsX", &profileNameResYvsX, &b_TkOffTreeVariables_profileNameResYvsX);
   fChain->SetBranchAddress("profileNameResYvsY", &profileNameResYvsY, &b_TkOffTreeVariables_profileNameResYvsY);
   Notify();
}

/*--------------------------------------------------------------------*/
Bool_t OfflineValidationTreeAnalysis::Notify()
/*--------------------------------------------------------------------*/
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

/*--------------------------------------------------------------------*/
void OfflineValidationTreeAnalysis::Show(Long64_t entry)
/*--------------------------------------------------------------------*/
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

/*--------------------------------------------------------------------*/
void OfflineValidationTreeAnalysis::makeNicePlotStyle(TH1F *hist)
/*--------------------------------------------------------------------*/
{
  hist->SetLineWidth(2);
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  hist->GetXaxis()->SetTitleFont(42); 
  hist->GetYaxis()->SetTitleFont(42);  
  hist->GetXaxis()->SetTitleSize(0.045);
  hist->GetYaxis()->SetTitleSize(0.045);
  hist->GetXaxis()->SetTitleOffset(1.);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(.035);
  hist->GetXaxis()->SetLabelSize(.035);
  hist->SetMarkerSize(1.4);
}

/*--------------------------------------------------------------------*/
void OfflineValidationTreeAnalysis::makeNice2DPlotStyle(TH2F *hist)
/*--------------------------------------------------------------------*/
{
  hist->SetLineWidth(3);
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  hist->GetXaxis()->SetTitleFont(42); 
  hist->GetYaxis()->SetTitleFont(42);  
  hist->GetXaxis()->SetTitleSize(0.045);
  hist->GetYaxis()->SetTitleSize(0.045);
  hist->GetXaxis()->SetTitleOffset(1.);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(.035);
  hist->GetXaxis()->SetLabelSize(.035);
}

/*--------------------------------------------------------------------*/
void OfflineValidationTreeAnalysis::makeNiceProfileStyleAndColor(TProfile *hist,int color)
/*--------------------------------------------------------------------*/
{

  Int_t colors[8]={kRed,kBlack,kBlue,kMagenta,kGreen,kCyan,kViolet,kPink};
  Int_t styles[8]={kFullCircle,kOpenCircle,kFullSquare,kOpenSquare,kFullTriangleUp,kOpenTriangleUp,kFullTriangleDown,kOpenTriangleDown};

  hist->SetLineWidth(3);
  hist->SetLineColor(colors[color]);
  //hist->SetFillColor(colors[color]);
  hist->SetMarkerColor(colors[color]);
  hist->SetMarkerStyle(styles[color]);
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  hist->GetXaxis()->SetTitleFont(42); 
  hist->GetYaxis()->SetTitleFont(42);  
  hist->GetXaxis()->SetTitleSize(0.045);
  hist->GetYaxis()->SetTitleSize(0.045);
  hist->GetXaxis()->SetTitleOffset(1.);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(.045);
  hist->GetXaxis()->SetLabelSize(.045);
  hist->SetMarkerSize(1.8);
}

/*--------------------------------------------------------------------*/
void OfflineValidationTreeAnalysis::makeNicePlotStyleAndColor(TH1F *hist,int color)
/*--------------------------------------------------------------------*/
{

  Int_t colors[8]={kRed,kBlack,kBlue,kMagenta,kGreen,kCyan,kViolet,kPink};
  Int_t styles[8]={kFullCircle,kOpenCircle,kFullSquare,kOpenSquare,kFullTriangleUp,kOpenTriangleUp,kFullTriangleDown,kOpenTriangleDown};

  hist->SetLineWidth(3);
  hist->SetLineColor(colors[color]);
  hist->SetMarkerColor(colors[color]);
  hist->SetMarkerStyle(styles[color]);
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  hist->GetXaxis()->SetTitleFont(42); 
  hist->GetYaxis()->SetTitleFont(42);  
  hist->GetXaxis()->SetTitleSize(0.045);
  hist->GetYaxis()->SetTitleSize(0.045);
  hist->GetXaxis()->SetTitleOffset(1.);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(.045);
  hist->GetXaxis()->SetLabelSize(.045);
  hist->SetMarkerSize(1.8);
}

/*--------------------------------------------------------------------*/
void OfflineValidationTreeAnalysis::setStyle(TString palettename)
/*--------------------------------------------------------------------*/
{

  TH1::StatOverflows(kTRUE);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.1);
  //gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadRightMargin(0.0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFont(42);
  //gStyle->SetTitleColor(kBlue);
  //gStyle->SetTitleTextColor(kBlue);
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.05);///---> gStyle->SetStatFontSize(0.025);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptFit(1);
  gStyle->SetNdivisions(510);

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  Double_t stops[NRGBs];
  Double_t red[NRGBs];
  Double_t green[NRGBs];
  Double_t blue[NRGBs];
  
  if (palettename == "halfgray"){
    Double_t stops1[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red1[5] = {1.00, 0.91, 0.80, 0.67, 1.00};
    Double_t green1[5] = {1.00, 0.91, 0.80, 0.67, 1.00};
    Double_t blue1[5] = {1.00, 0.91, 0.80, 0.67, 1.00};
    
    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }

  } else if(palettename == "gray"){
    Double_t stops1[5] = {0.00, 0.01, 0.05, 0.09, 0.1};
    Double_t red1[5] = {1.00, 0.84, 0.61, 0.34, 0.00};
    Double_t green1[5] = {1.00, 0.84, 0.61, 0.34, 0.00};
    Double_t blue1[5] = {1.00, 0.84, 0.61, 0.34, 0.00};

    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }
    
  } else if(palettename == "blues"){
    Double_t stops1[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red1[5] = {1.00, 0.84, 0.61, 0.34, 0.00};
    Double_t green1[5] = {1.00, 0.84, 0.61, 0.34, 0.00};
    Double_t blue1[5] = {1.00, 1.00, 1.00, 1.00, 1.00};
    
    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }
    
  } else if(palettename == "reds"){
    Double_t stops1[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red1[5] = {1.00, 1.00, 1.00, 1.00, 1.00};
    Double_t green1[5] = {1.00, 0.84, 0.61, 0.34, 0.00};
    Double_t blue1[5] = {1.00, 0.84, 0.61, 0.34, 0.00};
    
    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }

  } else if(palettename == "antigray"){
    Double_t stops1[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red1[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t green1[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t blue1[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
    
    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }

  } else if(palettename == "fire"){
   
    Double_t stops1[5] = {0.00, 0.25, 0.50, 0.75, 1.00};
    Double_t red1[5]   = {1.00, 1.00, 1.00, 0.50, 0.00};
    Double_t green1[5] = {1.00, 1.00, 0.00, 0.00, 0.00};
    Double_t blue1[5]  = {0.30, 0.30, 0.00, 0.00, 0.00};

    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }

  } else if(palettename == "antifire"){
    Double_t stops1[5] = {0.00, 0.20, 0.80, 1.00};
    Double_t red1[5] = {0.50, 1.00, 1.00, 1.00};
    Double_t green1[5] = {0.00, 0.00, 1.00, 1.00};
    Double_t blue1[5] = {0.00, 0.00, 0.00, 0.20};
    
    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }

  } else if(palettename == "logredblue") {
    
    Double_t stops1[5] = {0.0001, 0.0010, 0.0100, 0.1000, 1.0000};
    Double_t red1[5] = {1.00, 0.75, 0.50, 0.25, 0.00};
    Double_t green1[5] = {0.00, 0.00, 0.00, 0.00, 0.00};
    Double_t blue1[5] = {0.00, 0.25, 0.50, 0.75, 1.00};
    
    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }
    
  } else if(palettename == "logbluered") {

    Double_t stops1[5] = {0.0001, 0.0010, 0.0100, 0.1000, 1.0000};
    Double_t red1[5] = {0.00, 0.25, 0.50, 0.75, 1.00};
    Double_t green1[5] = {0.00, 0.00, 0.00, 0.00, 0.00};
    Double_t blue1[5] = {1.00, 0.75, 0.50, 0.25, 0.00};
    
    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }
    
  } else{
    // default palette, looks cool
    Double_t stops1[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red1[5] = {0.00, 0.00, 0.87, 1.00, 0.51};
    Double_t green1[5] = {0.00, 0.81, 1.00, 0.20, 0.00};
    Double_t blue1[5] = {0.51, 1.00, 0.12, 0.00, 0.00};
    
    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }
  }
  
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

}

/*--------------------------------------------------------------------*/
TPaveText* OfflineValidationTreeAnalysis::addNiceTLegend(TObjArray *histos)
/*--------------------------------------------------------------------*/
{

  const Double_t cm_To_um = 10000;
  TPaveText *leg = new TPaveText(0.096,0.723,0.406,0.955,"NDC"); 
  TIter myiter(histos);
  TH1F* h;
  while((h=(TH1F*)myiter.Next())){
    TString LEGO;
    TString tmp = "#splitline{"+TString(h->GetTitle())+" (#mum)}{mean = %.2f, rms = %.2f}";

    LEGO.Form(tmp,cm_To_um*(h->GetMean()),cm_To_um*(h->GetRMS()));

    //std::cout<<h->GetMean()<<std::endl;

    TText *le1 = leg->AddText(LEGO); 
    le1->SetTextFont(42);
    le1->SetTextAlign(11);
    le1->SetTextSize(0.035);
    le1->SetTextColor(h->GetLineColor());

    //std::cout<<h->GetTitle()<<std::endl;
    //LEGO.Clear();
    //tmp.Clear();
  }
   
  leg->SetFillColor(10);
  leg->SetLineColor(10);
  leg->SetShadowColor(0);
  
  return leg;

}

/*--------------------------------------------------------------------*/
void OfflineValidationTreeAnalysis::makeNiceCanv(TCanvas *c)
/*--------------------------------------------------------------------*/
{
  
  c->cd()->SetBottomMargin(0.11);
  c->cd()->SetLeftMargin(0.11);
  c->cd()->SetRightMargin(0.07);
  c->cd()->SetTopMargin(0.10);  

}


#endif // #ifdef OfflineValidationTreeAnalysis_Cxx
