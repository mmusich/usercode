//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 17 11:11:48 2012 by ROOT version 5.27/06b
// from TTree JetByJetComparison/file1 vs file2
// found on file: JetByJetComparison.root
//////////////////////////////////////////////////////////

#ifndef JetByJetComparison_h
#define JetByJetComparison_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <TH2.h>
#include <TH1.h>

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
  Float_t         PVx[2];
  Float_t         PVy[2];
  Float_t         PVz[2];
  Float_t         PVChi2[2];
  Float_t         PVndof[2];
  Float_t         PVNormalizedChi2[2];
  Char_t          isBGluonSplitting;
  Char_t          isCGluonSplitting;
  Float_t         pt[2];
  Float_t         eta[2];
  Float_t         phi[2];
  Int_t           jetId[2];
  Int_t           jetnTracks[2];
  Int_t           MCTrueFlavor[2];
  Float_t         SV3dDistance[2];
  Float_t         SV3dDistanceError[2];
  Float_t         SV2dDistance[2];
  Float_t         SV2dDistanceError[2];
  Float_t         SVChi2[2];
  Float_t         SVDegreesOfFreedom[2];
  Float_t         SVNormChi2[2];
  Float_t         SVMass[2];
  Int_t           SVtotCharge[2];
  Int_t           SVnVertices[2];
  Int_t           SVnVertexTracks[2];
  Int_t           SVnVertexTracksAll[2];
  Float_t         IP3d1[2];
  Float_t         IP3dError1[2];
  Float_t         IP3dDecayLength1[2];
  Float_t         IP3dTransverseMomentum1[2];
  Float_t         IP3dEta1[2];
  Float_t         IP3dPhi1[2];
  Float_t         IP3d2[2];
  Float_t         IP3dError2[2];
  Float_t         IP3dDecayLength2[2];
  Float_t         IP3dTransverseMomentum2[2];
  Float_t         IP3dEta2[2];
  Float_t         IP3dPhi2[2];
  Float_t         IP3d3[2];
  Float_t         IP3dError3[2];
  Float_t         IP3dDecayLength3[2];
  Float_t         IP3dTransverseMomentum3[2];
  Float_t         IP3dEta3[2];
  Float_t         IP3dPhi3[2];
  Float_t         IP3d4[2];
  Float_t         IP3dError4[2];
  Float_t         IP3dDecayLength4[2];
  Float_t         IP3dTransverseMomentum4[2];
  Float_t         IP3dEta4[2];
  Float_t         IP3dPhi4[2];                                   
  Float_t         tche[2];
  Float_t         tchp[2];
  Float_t         ssvhe[2];
  Float_t         ssvhp[2];
  Float_t         csv[2];
  Float_t         jp[2];
  Float_t         jbp[2];
  
   // List of branches
  TBranch        *b_Trun;   //!
  TBranch        *b_Tlumi;   //!
  TBranch        *b_Tevent;   //!
  TBranch        *b_Tpthat;   //!
  TBranch        *b_Tmcweight;   //!
  TBranch        *b_TPVx;   //!
  TBranch        *b_TPVy;   //!
  TBranch        *b_TPVz;   //!
  TBranch        *b_TisBGluonSplitting;   //!
  TBranch        *b_TisCGluonSplitting;   //!
  TBranch        *b_Tpt;   //!
  TBranch        *b_Teta;   //!
  TBranch        *b_Tphi;   //!
  TBranch        *b_TjetId;   //!
  TBranch        *b_TnTracks;   //!
  TBranch        *b_TMCTrueFlavor;   //!
  TBranch        *b_TSV3dDistance;   //!
  TBranch        *b_TSV3dDistanceError;   //!
  TBranch        *b_TSV2dDistance;   //!
  TBranch        *b_TSV2dDistanceError;   //!
  TBranch        *b_TSVChi2;   //!
  TBranch        *b_TSVDegreesOfFreedom;   //!
  TBranch        *b_TSVNormChi2;   //!
  TBranch        *b_TSVMass;   //!
  TBranch        *b_TSVtotCharge;   //!
  TBranch        *b_TSVnVertices;   //!
  TBranch        *b_TSVnVertexTracks;   //!
  TBranch        *b_TSVnVertexTracksAll;   //!
  TBranch        *b_TIP3dFirst;   //!
  TBranch        *b_TIP3dErrorFirst;   //!
  TBranch        *b_TIP3dDecayLengthFirst;   //!
  TBranch        *b_TIP3dTransverseMomentumFirst;   //!
  TBranch        *b_TIP3dEtaFirst;   //!
  TBranch        *b_TIP3dPhiFirst;   //!
  TBranch        *b_TIP3dSecond;   //!
  TBranch        *b_TIP3dErrorSecond;   //!
  TBranch        *b_TIP3dDecayLengthSecond;   //!
  TBranch        *b_TIP3dTransverseMomentumSecond;   //!
  TBranch        *b_TIP3dEtaSecond;   //!
  TBranch        *b_TIP3dPhiSecond;   //!
  TBranch        *b_TIP3dThird;   //!
  TBranch        *b_TIP3dErrorThird;   //!
  TBranch        *b_TIP3dDecayLengthThird;   //!
  TBranch        *b_TIP3dTransverseMomentumThird;   //!
  TBranch        *b_TIP3dEtaThird;   //!
  TBranch        *b_TIP3dPhiThird;   //!
  TBranch        *b_TIP3dFourth;   //!
  TBranch        *b_TIP3dErrorFourth;   //!
  TBranch        *b_TIP3dDecayLengthFourth;   //!
  TBranch        *b_TIP3dTransverseMomentumFourth;   //!
  TBranch        *b_TIP3dEtaFourth;   //!
  TBranch        *b_TIP3dPhiFourth;   //!
  TBranch        *b_Ttche;   //!
  TBranch        *b_Ttchp;   //!
  TBranch        *b_Tssvhe;   //!
  TBranch        *b_Tssvhp;   //!
  TBranch        *b_Tcsv;   //!
  TBranch        *b_Tjp;   //!
  TBranch        *b_Tjbp;   //!
  
  JetByJetComparison(TString filename,TTree *tree=0);
  virtual ~JetByJetComparison();
  //virtual Int_t  Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
 
  virtual void     setTDRStyle();
  virtual void     cmsPrel(const double& intLumi);

  virtual void     AddHisto(vector<TH1F*> &Histo1DB, string name, string title,const int& nbins, const Float_t& min, const Float_t& max);
  virtual void     AddHisto2D(vector<TH2F*> &Histo2DB, string name, string title,TString firstCond,TString secondCond ,const int& nbins, const Float_t& min, const Float_t& max, const int& nbinsy, const Float_t& miny, const Float_t& maxy);
  virtual TH1F*    findHisto(vector<TH1F*> &Histo1DB,TString keyword);
  virtual TH2F*    findHisto2D(vector<TH2F*> &Histo2DB,TString keyword);

};

#endif

#ifdef JetByJetComparison_cxx
JetByJetComparison::JetByJetComparison(TString filename,TTree *tree)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("JetByJetComparison.root");
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (!f) {
	f = new TFile(filename);
      }
    tree = (TTree*)gDirectory->Get("JetByJetComparison");
    
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
   fChain->SetMakeClass(1);
   
   fChain->SetBranchAddress("run", &run, &b_Trun);
   fChain->SetBranchAddress("lumi", &lumi, &b_Tlumi);
   fChain->SetBranchAddress("evt", &evt, &b_Tevent);
   fChain->SetBranchAddress("pthat", &pthat, &b_Tpthat);
   fChain->SetBranchAddress("mcweight", &mcweight, &b_Tmcweight);
   fChain->SetBranchAddress("PVx", &PVx, &b_TPVx);
   fChain->SetBranchAddress("PVy", &PVy, &b_TPVy);
   fChain->SetBranchAddress("PVz", &PVz, &b_TPVz);
   fChain->SetBranchAddress("PVChi2", &PVChi2, &b_TPVx);
   fChain->SetBranchAddress("PVndof", &PVndof, &b_TPVy);
   fChain->SetBranchAddress("PVNormalizedChi2", &PVNormalizedChi2, &b_TPVz);
   fChain->SetBranchAddress("isBGluonSplitting", &isBGluonSplitting, &b_TisBGluonSplitting);
   fChain->SetBranchAddress("isCGluonSplitting", &isCGluonSplitting, &b_TisCGluonSplitting);
   fChain->SetBranchAddress("isCGluonSplitting", &isCGluonSplitting, &b_TisCGluonSplitting);
   fChain->SetBranchAddress("pt", pt, &b_Tpt);
   fChain->SetBranchAddress("eta", eta, &b_Teta);
   fChain->SetBranchAddress("phi", phi, &b_Tphi);
   fChain->SetBranchAddress("jetId", jetId, &b_TjetId);
   fChain->SetBranchAddress("jetnTracks", jetnTracks, &b_TnTracks);
   fChain->SetBranchAddress("MCTrueFlavor", MCTrueFlavor, &b_TMCTrueFlavor);
   fChain->SetBranchAddress("SV3dDistance", SV3dDistance, &b_TSV3dDistance);
   fChain->SetBranchAddress("SV3dDistanceError", SV3dDistanceError, &b_TSV3dDistanceError);
   fChain->SetBranchAddress("SV2dDistance", SV2dDistance, &b_TSV2dDistance);
   fChain->SetBranchAddress("SV2dDistanceError", SV2dDistanceError, &b_TSV2dDistanceError);
   fChain->SetBranchAddress("SVChi2", SVChi2, &b_TSVChi2);
   fChain->SetBranchAddress("SVDegreesOfFreedom", SVDegreesOfFreedom, &b_TSVDegreesOfFreedom);
   fChain->SetBranchAddress("SVNormChi2", SVNormChi2, &b_TSVNormChi2);
   fChain->SetBranchAddress("SVMass", SVMass, &b_TSVMass);
   fChain->SetBranchAddress("SVtotCharge", SVtotCharge, &b_TSVtotCharge);
   fChain->SetBranchAddress("SVnVertices", SVnVertices, &b_TSVnVertices);
   fChain->SetBranchAddress("SVnVertexTracks", SVnVertexTracks, &b_TSVnVertexTracks);
   fChain->SetBranchAddress("SVnVertexTracksAll", SVnVertexTracksAll, &b_TSVnVertexTracksAll);
   fChain->SetBranchAddress("IP3d1", IP3d1, &b_TIP3dFirst);
   fChain->SetBranchAddress("IP3dError1", IP3dError1, &b_TIP3dErrorFirst);
   fChain->SetBranchAddress("IP3dDecayLength1", IP3dDecayLength1, &b_TIP3dDecayLengthFirst);
   fChain->SetBranchAddress("IP3dTransverseMomentum1", IP3dTransverseMomentum1, &b_TIP3dTransverseMomentumFirst);
   fChain->SetBranchAddress("IP3dEta1", IP3dEta1, &b_TIP3dEtaFirst);
   fChain->SetBranchAddress("IP3dPhi1", IP3dPhi1, &b_TIP3dPhiFirst);
   fChain->SetBranchAddress("IP3d2", IP3d2, &b_TIP3dSecond);
   fChain->SetBranchAddress("IP3dError2", IP3dError2, &b_TIP3dErrorSecond);
   fChain->SetBranchAddress("IP3dDecayLength2", IP3dDecayLength2, &b_TIP3dDecayLengthSecond);
   fChain->SetBranchAddress("IP3dTransverseMomentum2", IP3dTransverseMomentum2, &b_TIP3dTransverseMomentumSecond);
   fChain->SetBranchAddress("IP3dEta2", IP3dEta2, &b_TIP3dEtaSecond);
   fChain->SetBranchAddress("IP3dPhi2", IP3dPhi2, &b_TIP3dPhiSecond);
   fChain->SetBranchAddress("IP3d3", IP3d3, &b_TIP3dThird);
   fChain->SetBranchAddress("IP3dError3", IP3dError3, &b_TIP3dErrorThird);
   fChain->SetBranchAddress("IP3dDecayLength3", IP3dDecayLength3, &b_TIP3dDecayLengthThird);
   fChain->SetBranchAddress("IP3dTransverseMomentum3", IP3dTransverseMomentum3, &b_TIP3dTransverseMomentumThird);
   fChain->SetBranchAddress("IP3dEta3", IP3dEta3, &b_TIP3dEtaThird);
   fChain->SetBranchAddress("IP3dPhi3", IP3dPhi3, &b_TIP3dPhiThird);
   fChain->SetBranchAddress("IP3d4", IP3d4, &b_TIP3dFourth);
   fChain->SetBranchAddress("IP3dError4", IP3dError4, &b_TIP3dErrorFourth);
   fChain->SetBranchAddress("IP3dDecayLength4", IP3dDecayLength4, &b_TIP3dDecayLengthFourth);
   fChain->SetBranchAddress("IP3dTransverseMomentum4", IP3dTransverseMomentum4, &b_TIP3dTransverseMomentumFourth);
   fChain->SetBranchAddress("IP3dEta4", IP3dEta4, &b_TIP3dEtaFourth);
   fChain->SetBranchAddress("IP3dPhi4", IP3dPhi4, &b_TIP3dPhiFourth);
   fChain->SetBranchAddress("tche", tche, &b_Ttche);
   fChain->SetBranchAddress("tchp", tchp, &b_Ttchp);
   fChain->SetBranchAddress("ssvhe", ssvhe, &b_Tssvhe);
   fChain->SetBranchAddress("ssvhp", ssvhp, &b_Tssvhp);
   fChain->SetBranchAddress("csv", csv, &b_Tcsv);
   fChain->SetBranchAddress("jp", jp, &b_Tjp);
   fChain->SetBranchAddress("jbp", jbp, &b_Tjbp);
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
// Int_t JetByJetComparison::Cut(Long64_t entry)
// {
// // This function may be called from Loop.
// // returns  1 if entry is accepted.
// // returns -1 otherwise.
//    return 1;
// }
#endif // #ifdef JetByJetComparison_cxx
