#include "TSystem.h"
#include <Riostream.h>
#include <TROOT.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TList.h>
#include <TTree.h>
#include <TF1.h>
#include <TFile.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TKey.h>
#include <TMath.h>
#include <Zbbstruct4JEC.h>
#include <Math/VectorUtil.h>
#include <fstream>
#include <iostream>

using namespace std;

// graphics
void setStyle();
Double_t getMaximum(TH1 *h1,TH1 *h2);   
Double_t getMaximum(TH1 *h,THStack *hs);
TLegend* MakeTLegend(TH1 *h,THStack *hs,TString *LegLabels);
void MakeNiceHistoStyle(TH1F *hist, Int_t color=-999);
void MakeNiceCanvas(TCanvas *canvas);
TCanvas* DrawCanvasWithRatio(TCanvas* canvas, bool emptyBinsNoUnc = true);

void CMSPrel(Double_t ZpT_min,Double_t ZpT_max, Double_t Lumi,Double_t alpha=-999);

Double_t* RatioExtrap(Bool_t MPF = false, Int_t ratioOpt = 1,Double_t ZpT_min = 10., Double_t ZpT_max = 1000., Bool_t plotspectra = false, Bool_t plottrend = false);

Bool_t isGoodForJECStudy(Z_info Z, bjet_info bjet, jet2_info jet2, Double_t ZpT_min = 10., Double_t ZpT_max = 1000.);

void ScanTheTree4JEC_bis(Bool_t MPF = true, Bool_t integrZpT = false, Bool_t plotspectra = false, Bool_t plottrend = false, Double_t ZpT_min = 10., Double_t ZpT_max = 1000.){
  
  const int n_ZpT_bin(4);
  Double_t ZpT_bin_bound[n_ZpT_bin+1] = {ZpT_min,50,100,200,ZpT_max};

  if (integrZpT){
    Double_t* ratio_int = new Double_t[2];
    ratio_int = RatioExtrap(MPF, 1, ZpT_min, ZpT_max, plotspectra, plottrend);
    cout<<"Integrating over all the range of Z pT"<<endl;
    cout<<"DATA/MC (extrapolated) = "<<ratio_int[0]<<" +- "<<ratio_int[1]<<endl;
  }
  else {
    Double_t ratio[n_ZpT_bin]; Double_t ratio_err[n_ZpT_bin];
    Double_t ratio_lowstat[n_ZpT_bin]; Double_t ratio_lowstat_err[n_ZpT_bin];
    Double_t ratio_a0a02[n_ZpT_bin]; Double_t ratio_a0a02_err[n_ZpT_bin];
    for (Int_t i=0; i<n_ZpT_bin; ++i){
      Double_t* ratioPluserror_1 = RatioExtrap(MPF, 1, ZpT_bin_bound[i], ZpT_bin_bound[i+1], false, plottrend);
      ratio[i]=ratioPluserror_1[0];
      ratio_err[i]=ratioPluserror_1[1];
      Double_t* ratioPluserror_2 = RatioExtrap(MPF, 2, ZpT_bin_bound[i], ZpT_bin_bound[i+1], false, plottrend);
      ratio_lowstat[i]=ratioPluserror_2[0];
      ratio_lowstat_err[i]=ratioPluserror_2[1];
      Double_t* ratioPluserror_3 = RatioExtrap(MPF, 3, ZpT_bin_bound[i], ZpT_bin_bound[i+1], false, plottrend);
      ratio_a0a02[i]=ratioPluserror_3[0];
      ratio_a0a02_err[i]=ratioPluserror_3[1];
    }

    TF1 *constant_1 = new TF1("constant_1","[0]");
    TF1 *constant_2 = new TF1("constant_2","[0]");
    TF1 *constant_3 = new TF1("constant_3","[0]");
    constant_1->SetLineWidth(1.2);
    constant_2->SetLineWidth(1.2);
    constant_3->SetLineWidth(1.2);
    Double_t ZpT[n_ZpT_bin];
    Double_t ZpT_err[n_ZpT_bin]; 
    for (Int_t i=0; i<n_ZpT_bin; ++i) ZpT_err[i]=0;
    for (Int_t j=0; j<n_ZpT_bin; ++j) ZpT[j] = ZpT_bin_bound[j] + (ZpT_bin_bound[j+1]-ZpT_bin_bound[j])/2;
    TGraphErrors *DATAMC_1 = new TGraphErrors(n_ZpT_bin, ZpT, ratio, ZpT_err, ratio_err);
    TGraphErrors *DATAMC_2 = new TGraphErrors(n_ZpT_bin, ZpT, ratio_lowstat, ZpT_err, ratio_lowstat_err);
    TGraphErrors *DATAMC_3 = new TGraphErrors(n_ZpT_bin, ZpT, ratio_a0a02, ZpT_err, ratio_a0a02_err);
    //DATAMC_1->SetTitle("Jet Energy Response DATA/MC ratio");
    DATAMC_1->GetYaxis()->SetRangeUser(0.7,1.3);
    DATAMC_1->GetYaxis()->SetTitle("DATA/MC (#alpha->0)");
    DATAMC_1->GetXaxis()->SetTitle("p_{T}^{Z} (GeV)");
    DATAMC_1->SetMarkerStyle(20);
    DATAMC_1->SetMarkerSize(1.2);
    DATAMC_1->SetMarkerColor(2);
    DATAMC_1->Fit(constant_1, "QC");
    DATAMC_2->GetYaxis()->SetRangeUser(0.7,1.3);
    DATAMC_2->GetYaxis()->SetTitle("DATA/MC (#alpha=0.3)");
    DATAMC_2->GetXaxis()->SetTitle("p_{T}^{Z} (GeV)");
    DATAMC_2->SetMarkerStyle(20);
    DATAMC_2->SetMarkerSize(1.2);
    DATAMC_2->SetMarkerColor(2);
    DATAMC_2->Fit(constant_2, "QC");
    DATAMC_3->GetYaxis()->SetRangeUser(0.7,1.3);
    DATAMC_3->GetYaxis()->SetTitle("DATA/MC (#alpha->0/#alpha=0.2)");
    DATAMC_3->GetXaxis()->SetTitle("p_{T}^{Z} (GeV)");
    DATAMC_3->SetMarkerStyle(20);
    DATAMC_3->SetMarkerSize(1.2);
    DATAMC_3->SetMarkerColor(2);
    DATAMC_3->Fit(constant_3, "QC");
    //TLegend* leg2 = new TLegend(0.65,0.7,0.9,0.87);
    //leg2->SetFillColor(0);
    //leg2->SetShadowColor(10);
    //leg2->SetTextFont(42);
    //leg2->AddEntry(constant_1,"Mean R_{DATA/MC}","l");
    new TCanvas;
    DATAMC_1->Draw("AP");
    new TCanvas;
    DATAMC_2->Draw("AP");
    new TCanvas;
    DATAMC_3->Draw("AP");
    //leg2->Draw("same");
    cout<<"Constant fit to the values per-bin of ZpT (standard extrapolation alpha->0)"<<endl;
    cout<<"DATA/MC (extrapolated) = "<<constant_1->GetParameter(0)<<" +- "<<constant_1->GetParError(0)<<endl;
    cout<<"Cabs = "<<1/constant_1->GetParameter(0)<<" +- "<<constant_1->GetParError(0)/((constant_1->GetParameter(0))*(constant_1->GetParameter(0)))<<endl;
    cout<<"Constant fit to the values per-bin of ZpT (face value at alpha=0.3)"<<endl;
    cout<<"DATA/MC (extrapolated) = "<<constant_2->GetParameter(0)<<" +- "<<constant_2->GetParError(0)<<endl;
    cout<<"Cabs = "<<1/constant_2->GetParameter(0)<<" +- "<<constant_2->GetParError(0)/((constant_2->GetParameter(0))*(constant_2->GetParameter(0)))<<endl;
    cout<<"Constant fit to the values per-bin of ZpT (alpha->0/alpha=0.2 ratio)"<<endl;
    cout<<"DATA/MC (extrapolated) = "<<constant_3->GetParameter(0)<<" +- "<<constant_3->GetParError(0)<<endl;
    cout<<"Cabs = "<<1/constant_3->GetParameter(0)<<" +- "<<constant_3->GetParError(0)/((constant_3->GetParameter(0))*(constant_3->GetParameter(0)))<<endl;
  }

}

Double_t* RatioExtrap(Bool_t MPF, Int_t ratioOpt, Double_t ZpT_min, Double_t ZpT_max, Bool_t plotspectra, Bool_t plottrend){

  setStyle();
  TH1F::SetDefaultSumw2(kTRUE);

  //Define constants
  const int n_MCSamples(4);
  const int  n_alphaValue(4);
 
  //Declare (and initialize) objects and arrays of objects  
  Float_t alpha[n_alphaValue] = {0.1,0.15,0.2,0.3};
  Float_t alpha_err[n_alphaValue] = {0,0,0,0};
  TString alpha_value[n_alphaValue] = {"alpha 0.1", "alpha 0.15", "alpha 0.2", "alpha 0.3"};
  Double_t nEvts_DATA[n_alphaValue] = {0.,0.,0.,0.};
  Double_t nEvts_MC[n_alphaValue] = {0.,0.,0.,0.};
  TFile *MCfile[n_MCSamples];
  TFile *datafile;
  TTree *ntupleDATA;
  TTree *ntupleMC[n_MCSamples];
  //if (ZpT_bin < 0 || ZpT_bin > 4) cout<<"Invalid number of ZpT bin: select a bin from 0 to 4"<<endl;
  //Int_t ZpT_bin_iterator = ZpT_bin;

  //Declare structures
  Event_info *Event = new Event_info[n_MCSamples+1];
  Z_info     *Z     = new Z_info[n_MCSamples+1];
  bjet_info  *bjet  = new bjet_info[n_MCSamples+1];
  jet2_info  *jet2  = new jet2_info[n_MCSamples+1];
  MET_info   *MET   = new MET_info[n_MCSamples+1];

  //Declare and book control histograms
  TH1F* h_pTjb_DATA = new TH1F("h_pTjb_DATA","b-jet p_{T};leading b-jet p_{T} [GeV];events", 50, 0., 200);
  TH1F* h_pTj2_DATA = new TH1F("h_pTj2_DATA","jet2 p_{T};subleading jet p_{T} [GeV];events", 50, 0., 200);
  TH1F* h_pTZ_DATA = new TH1F("h_pTZ_DATA","Z p_{T};Z p_{T} [GeV]; events", 50, 0., 200);
  TH1F* h_pTj2pTZratio_DATA = new TH1F("h_pTj2pTZratio_DATA","jet2 p_{T}/Z p_{T}; p^{j2}_{T}/p^{Z}_{T}; events", 80, 0., 2);
  TH1F* h_MET_DATA = new TH1F("h_MET_DATA","#slash{E}_{T};total #slash{E}_{T} [GeV]; events", 50, 0., 200);
  TH1F* h_METangle_DATA = new TH1F("h_METangle_DATA","cos#theta(#slash{E}_{T},#hat{p^{Z}}_{T});cos#theta(#slash{E}_{T},#hat{p}^{Z}_{T});events",50,-1.,1.);

  TH1F* h_pTjb_MC = new TH1F("h_pTjb_MC", "b-jet p_{T};leading b-jet p_{T} [GeV];events", 50, 0., 200);
  TH1F* h_pTj2_MC = new TH1F("h_pTj2_MC", "jet2 p_{T};subleading jet p_{T} [GeV];events", 50, 0., 200);
  TH1F* h_pTZ_MC = new TH1F("h_pTZ_MC", "Z p_{T}; Z p_{T} [GeV]; events", 50, 0., 200);
  TH1F* h_pTj2pTZratio_MC = new TH1F("h_pTj2pTZratio_MC", "jet2 p_{T}/Z p_{T}; p^{j2}_{T}/p^{Z}_{T}; events", 80, 0., 2);
  TH1F* h_MET_MC = new TH1F("h_MET_MC","#slash{E}_{T};total #slash{E}_{T} [GeV]; events", 50, 0., 200);
  TH1F* h_METangle_MC = new TH1F("h_METangle_MC","cos#theta(#slash{E}_{T},#hat{p^{Z}}_{T});cos#theta(#slash{E}_{T},#hat{p}^{Z}_{T});events",50,-1.,1.);

  TH1F* h_pTjb_MC_forStack[n_MCSamples];
  TH1F* h_pTj2_MC_forStack[n_MCSamples];
  TH1F* h_pTZ_MC_forStack[n_MCSamples];
  TH1F* h_pTj2pTZratio_MC_forStack[n_MCSamples];
  TH1F* h_MET_MC_forStack[n_MCSamples];
  TH1F* h_METangle_MC_forStack[n_MCSamples];

  THStack* h_pTjb_MC_Stack = new THStack("h_pTjb_MC_stack", "h_pTjb_MC_stack;b-jet p_{T} (GeV); events");
  THStack* h_pTj2_MC_Stack= new THStack("h_pTj2_MC_stack", "h_pTj2_MC_stack;jet2 p_{T} (GeV); events");
  THStack* h_pTZ_MC_Stack= new THStack("h_pTZ_MC_stack", "h_pTZ_MC_stack;Z p_{T} (GeV); events"); 
  THStack* h_pTj2pTZratio_MC_Stack =new THStack("h_pTj2pTZratio_MC_stack", "h_pTj2pTZratio_MC_stack;jet2 p_{T}/Z p_{T}; events");
  THStack* h_MET_MC_Stack = new THStack("h_MET_MC_stack","#slash{E}_{T};total #slash{E}_{T} [GeV]; events");
  THStack* h_METangle_MC_Stack = new THStack("h_METangle_MC_Stack","cos#theta(#slash{E}_{T},#hat{p^{Z}}_{T});cos#theta(#slash{E}_{T},#hat{p}^{Z}_{T});events");

  for (Int_t j=0; j<n_MCSamples; ++j) {
    h_pTjb_MC_forStack[j] = new TH1F(Form("h_pTjb_MC_%i",j),"b-jet p_{T}; b-jet p_{T} (GeV); events", 50, 0., 200);
    h_pTj2_MC_forStack[j] = new TH1F(Form("h_pTj2_MC_%i",j),"jet2 p_{T}; (GeV); jet2 p_{T} events", 50, 0., 200);
    h_pTZ_MC_forStack[j] = new TH1F(Form("h_pTZ_MC_%i",j),"Z p_{T}; Z p_{T} (GeV); events", 50, 0., 200);
    h_pTj2pTZratio_MC_forStack[j] = new TH1F(Form("h_pTj2pTZratio_MC_%i",j),"jet2 p_{T}/Z p_{T}; jet2 p_{T}/Z p_{T};events", 80, 0., 2.);
    h_MET_MC_forStack[j] =  new TH1F(Form("h_MET_MC_%i",j),"#slash{E}_{T};total #slash{E}_{T} [GeV]; events", 50, 0., 200);
    h_METangle_MC_forStack[j] =  new TH1F(Form("h_METangle_MC_%i",j),"cos#theta(#slash{E}_{T},#hat{p^{Z}}_{T});cos#theta(#slash{E}_{T},#hat{p}^{Z}_{T});events",50,-1.,1.);
  }

  //-----------------------------------------------------------
  //Discriminating variables
  TH1F* h_pTbal_DATA[n_alphaValue], *h_Rmpf_DATA[n_alphaValue];
  TH1F* h_pTbal_MC[n_alphaValue], *h_Rmpf_MC[n_alphaValue];
  THStack* h_pTbal_MC_Stack[n_alphaValue], *h_Rmpf_MC_Stack[n_alphaValue];
  for (Int_t i=0; i<n_alphaValue; ++i) {
    h_pTbal_DATA[i] = new TH1F("h_pTbal_DATA_"+alpha_value[i],"p_{T} b-jet/p_{T} Z for balanced events alpha="+alpha_value[i]+";p_{T} b-jet/p_{T} Z;events",25,0,2.5);
    h_Rmpf_DATA[i] = new TH1F("h_Rmpf_DATA_"+alpha_value[i],"R_{MPF} spectrum for balanced events alpha="+alpha_value[i]+";R_{MPF};events",25,0,2.5);
    h_pTbal_MC[i] = new TH1F("h_pTbal_MC_"+alpha_value[i],"p_{T} b-jet/p_{T} Z for balanced events alpha="+alpha_value[i]+";p_{T} b-jet/p_{T} Z;events",25,0,2.5);
    h_Rmpf_MC[i] = new TH1F("h_Rmpf_MC_"+alpha_value[i],"R_{MPF} spectrum for balanced events alpha="+alpha_value[i]+";R_{MPF};events",25,0,2.5);
    h_pTbal_MC_Stack[i]= new THStack("h_pTbal_MC_Stack_"+alpha_value[i],"p_{T} b-jet/p_{T} Z for balanced events alpha="+alpha_value[i]);
    h_Rmpf_MC_Stack[i]= new THStack("h_Rmpf_MC_Stack_"+alpha_value[i],"R_{MPF} spectrum for balanced events alpha="+alpha_value[i]);
  }

  //Stack variables for MC
  TH1F* h_pTbal_MC_forStack[n_MCSamples][n_alphaValue], *h_Rmpf_MC_forStack[n_MCSamples][n_alphaValue];
  for (Int_t j=0; j<n_MCSamples; ++j) {
    for (Int_t i=0; i<n_alphaValue; ++i) {
      h_pTbal_MC_forStack[j][i] = new TH1F(Form("h_pTbal_MC_%i_",j)+alpha_value[i],"p_{T} b-jet/p_{T} Z for balanced events alpha="+alpha_value[i]+";p_{T} b-jet/p_{T} Z;events",25,0,2.5);
      h_Rmpf_MC_forStack[j][i]  = new TH1F(Form("h_Rmpf_MC_%i_",j)+alpha_value[i],"R_{MPF} spectrum for balanced events alpha="+alpha_value[i]+";R_{MPF};events",25,0,2.5);
    }
  }

  //Open files
  MCfile[3] = TFile::Open("analyzePAT_MC_Zb5fToLL_All.root","READ");
  MCfile[0] = TFile::Open("analyzePAT_MC_ZcToLL_All.root","READ");
  MCfile[1] = TFile::Open("analyzePAT_MC_ZlJets_All.root","READ");
  MCfile[2] = TFile::Open("analyzePAT_MC_TTJets_All.root","READ");
  datafile  = TFile::Open("analyzePAT_DATA2011_merged.root","READ");

  //**************************
  //SCAN NTUPLES
  //**************************
  //N-tuple
  ntupleDATA = (TTree*)datafile->Get("finaldistros_ssvhpt/JetBTag/ZbbNtuple4JEC");
  for (Int_t i=0; i<n_MCSamples; i++) ntupleMC[i] = (TTree*)MCfile[i]->Get("finaldistros_ssvhpt/JetBTag/ZbbNtuple4JEC");

  ntupleDATA->SetBranchAddress("Event",&Event[0]);
  ntupleDATA->SetBranchAddress("Z",&Z[0]);
  ntupleDATA->SetBranchAddress("bjet",&bjet[0]);
  ntupleDATA->SetBranchAddress("jet2",&jet2[0]);
  ntupleDATA->SetBranchAddress("MET",&MET[0]);
  for (Int_t i=0; i<n_MCSamples; i++){
    ntupleMC[i]->SetBranchAddress("Event",&Event[i+1]);
    ntupleMC[i]->SetBranchAddress("Z",&Z[i+1]);
    ntupleMC[i]->SetBranchAddress("bjet",&bjet[i+1]);
    ntupleMC[i]->SetBranchAddress("jet2",&jet2[i+1]);
    ntupleMC[i]->SetBranchAddress("MET",&MET[i+1]);
  }

  Int_t nEntries_DATA = ntupleDATA->GetEntries();
  Int_t nEntries_MC[n_MCSamples];
  for (Int_t i=0; i<n_MCSamples; i++) nEntries_MC[i] = ntupleMC[i]->GetEntries();

  //Compute weights for MC samples
  TString LegLabels[5]={"CMS DATA","MC Z+c","MC Z+l","MC t#bart","MC Z+b"};
  Double_t lumi = 2113;
  Double_t weights[n_MCSamples],xsections[n_MCSamples],nevts[n_MCSamples];
  for (Int_t i = 0; i < n_MCSamples; i++) nevts[i] = ((TH1F*)MCfile[i]->Get("analyzePat/Selevents/SelectedEvts"))->GetBinContent(1);
  xsections[3] = 3048;
  xsections[0] = 3048;
  xsections[1] = 3048;
  xsections[2] = 165;
  for (Int_t i = 0; i < n_MCSamples; i++) weights[i] = lumi*xsections[i]/nevts[i];

  //Scan DATA ntuple  
  for (Int_t iEntry=0; iEntry<nEntries_DATA; ++iEntry){ 
    ntupleDATA->GetEntry(iEntry);
    if ( isGoodForJECStudy(Z[0], bjet[0], jet2[0], ZpT_min, ZpT_max) ){    
      Float_t Rmpf = 1 + (MET[0].px_MET*Z[0].pX_Z + MET[0].py_MET*Z[0].pY_Z)/(Z[0].pT_Z*Z[0].pT_Z);
      h_pTjb_DATA->Fill(bjet[0].pT_bjet);
      h_pTj2_DATA->Fill(jet2[0].pT_jet2);
      h_pTZ_DATA->Fill(Z[0].pT_Z);
      h_pTj2pTZratio_DATA->Fill(jet2[0].pT_jet2/Z[0].pT_Z);
      h_MET_DATA->Fill(MET[0].pT_MET);
      h_METangle_DATA->Fill((MET[0].px_MET*Z[0].pX_Z + MET[0].py_MET*Z[0].pY_Z)/(Z[0].pT_Z*MET[0].pT_MET));

      for (Int_t j=0; j<n_alphaValue; ++j) {
	if ( jet2[0].pT_jet2 < alpha[j]*Z[0].pT_Z  ) { // alpha selection
	  nEvts_DATA[j]++;
	  h_pTbal_DATA[j]->Fill(bjet[0].pT_bjet/Z[0].pT_Z); 
	  h_Rmpf_DATA[j]->Fill(Rmpf);
	}
      } //loop on alpha
    }//cut on Z+b candidates      			   
  }//loop on tree entries
 
  //Int_t theStackColors[n_MCSamples]={100,95,90,80};
  //Int_t theStackColors[n_MCSamples]={90,96,99,46};
  Int_t theStackColors[n_MCSamples]={409,856,393,628};
 
  //Scan MC ntuple
  for (Int_t k=0; k<n_MCSamples; k++){ //loop over MC samples
    for (Int_t iEntry=0; iEntry<nEntries_MC[k]; ++iEntry){
      ntupleMC[k]->GetEntry(iEntry);
      Float_t w = Event[k+1].weight_Event*weights[k];
      if ( isGoodForJECStudy(Z[k+1], bjet[k+1], jet2[k+1], ZpT_min, ZpT_max) ){
	Float_t Rmpf = 1 + (MET[k+1].px_MET*Z[k+1].pX_Z + MET[k+1].py_MET*Z[k+1].pY_Z)/(Z[k+1].pT_Z*Z[k+1].pT_Z);
	h_pTjb_MC->Fill(bjet[k+1].pT_bjet,w);
	h_pTj2_MC->Fill(jet2[k+1].pT_jet2,w);
	h_pTZ_MC->Fill(Z[k+1].pT_Z,w);
	h_pTj2pTZratio_MC->Fill(jet2[k+1].pT_jet2/Z[k+1].pT_Z,w);
	h_MET_MC->Fill(MET[k+1].pT_MET,w);
	h_METangle_MC->Fill((MET[k+1].px_MET*Z[k+1].pX_Z + MET[k+1].py_MET*Z[k+1].pY_Z)/(Z[k+1].pT_Z*MET[k+1].pT_MET),w);
	
	h_pTjb_MC_forStack[k]->Fill(bjet[k+1].pT_bjet,w);
	h_pTj2_MC_forStack[k]->Fill(jet2[k+1].pT_jet2,w);
	h_pTZ_MC_forStack[k]->Fill(Z[k+1].pT_Z,w);
	h_pTj2pTZratio_MC_forStack[k]->Fill(jet2[k+1].pT_jet2/Z[k+1].pT_Z,w);
	h_MET_MC_forStack[k]->Fill(MET[k+1].pT_MET,w);
	h_METangle_MC_forStack[k]->Fill((MET[k+1].px_MET*Z[k+1].pX_Z + MET[k+1].py_MET*Z[k+1].pY_Z)/(Z[k+1].pT_Z*MET[k+1].pT_MET),w);

	for (Int_t j=0; j<n_alphaValue; ++j) {
	  if ( jet2[k+1].pT_jet2 < alpha[j]*Z[k+1].pT_Z  ) { // alpha selection
	    nEvts_MC[j]+=w;

	    h_pTbal_MC[j]->Fill(bjet[k+1].pT_bjet/Z[k+1].pT_Z,w); 
	    h_Rmpf_MC[j]->Fill(Rmpf,w);
	 
	    h_pTbal_MC_forStack[k][j]->Fill(bjet[k+1].pT_bjet/Z[k+1].pT_Z,w); 
	    h_Rmpf_MC_forStack[k][j]->Fill(Rmpf,w);
	  }
	} //loop on alpha
      } //cut on Z+b candidates      			   
    } //loop on tree entries

    h_pTjb_MC_Stack->Add(h_pTjb_MC_forStack[k]);
    h_pTj2_MC_Stack->Add(h_pTj2_MC_forStack[k]);
    h_pTZ_MC_Stack->Add(h_pTZ_MC_forStack[k]);
    h_pTj2pTZratio_MC_Stack->Add(h_pTj2pTZratio_MC_forStack[k]); 
    h_MET_MC_Stack->Add(h_MET_MC_forStack[k]);   
    h_METangle_MC_Stack->Add(h_METangle_MC_forStack[k]);

    for (Int_t j=0; j<n_alphaValue; ++j) {
      h_pTbal_MC_Stack[j]->Add(h_pTbal_MC_forStack[k][j]);
      h_Rmpf_MC_Stack[j]->Add(h_Rmpf_MC_forStack[k][j]);
    }
    
  }//loop on MC samples

  TLegend* theStackedLeg = MakeTLegend(h_pTjb_DATA,h_pTjb_MC_Stack,LegLabels);

  //**************************
  //PLOTS
  //**************************
  Float_t maxYvalue[n_alphaValue], maxYvalueJ2[n_alphaValue];
  for (Int_t i=0; i<n_alphaValue; ++i){
    if (h_pTbal_DATA[i]->GetMaximum()>h_pTbal_MC[i]->GetMaximum()) maxYvalueJ2[i]=h_pTbal_DATA[i]->GetMaximum();
    else maxYvalueJ2[i]=h_pTbal_MC[i]->GetMaximum();
    if ( MPF ) {
      if (h_Rmpf_DATA[i]->GetMaximum()>h_Rmpf_MC[i]->GetMaximum()) maxYvalue[i]=h_Rmpf_DATA[i]->GetMaximum();
      else maxYvalue[i]=h_Rmpf_MC[i]->GetMaximum();
    } else {
      if (h_pTbal_DATA[i]->GetMaximum()>h_pTbal_MC[i]->GetMaximum()) maxYvalue[i]=h_pTbal_DATA[i]->GetMaximum();
      else maxYvalue[i]=h_pTbal_MC[i]->GetMaximum();
    }
  }


  if (plotspectra){

    h_pTjb_MC->SetFillColor(46);
    h_pTjb_DATA->GetYaxis()->SetRangeUser(0.1,getMaximum(h_pTjb_DATA,h_pTjb_MC)*1.10);

    h_pTj2_MC->SetFillColor(46);
    h_pTj2_DATA->GetYaxis()->SetRangeUser(0.1,getMaximum(h_pTj2_DATA,h_pTj2_MC)*1.10);
   
    h_pTZ_MC->SetFillColor(46);
    h_pTZ_DATA->GetYaxis()->SetRangeUser(0.1,getMaximum(h_pTZ_DATA,h_pTZ_MC)*1.10);
 
    h_pTj2pTZratio_MC->SetFillColor(46);
    h_pTj2pTZratio_DATA->GetYaxis()->SetRangeUser(0.1,getMaximum(h_pTj2pTZratio_DATA,h_pTj2pTZratio_MC)*1.10);

    h_MET_MC->SetFillColor(46);
    h_MET_DATA->GetYaxis()->SetRangeUser(0.1,getMaximum(h_MET_DATA,h_MET_MC)*1.10);
    
    h_METangle_MC->SetFillColor(46);
    h_METangle_DATA->GetYaxis()->SetRangeUser(0.1,getMaximum(h_METangle_DATA,h_METangle_MC)*1.10);

    for (Int_t k=0; k<n_MCSamples; k++){ //loop over MC samples
      MakeNiceHistoStyle(h_pTjb_MC_forStack[k],theStackColors[k]);
      MakeNiceHistoStyle(h_pTj2_MC_forStack[k],theStackColors[k]); 
      MakeNiceHistoStyle(h_pTZ_MC_forStack[k],theStackColors[k]);  
      MakeNiceHistoStyle(h_pTj2pTZratio_MC_forStack[k],theStackColors[k]);
      MakeNiceHistoStyle(h_MET_MC_forStack[k],theStackColors[k]);
      MakeNiceHistoStyle(h_METangle_MC_forStack[k],theStackColors[k]);
    }
 
    MakeNiceHistoStyle(h_pTjb_DATA);
    MakeNiceHistoStyle(h_pTj2_DATA); 
    MakeNiceHistoStyle(h_pTZ_DATA);  
    MakeNiceHistoStyle(h_pTj2pTZratio_DATA);
    MakeNiceHistoStyle(h_MET_DATA); 
    MakeNiceHistoStyle(h_METangle_DATA);  

    TLegend* leg_pTjb = new TLegend(0.65,0.7,0.9,0.87);
    leg_pTjb->SetFillColor(0);
    leg_pTjb->SetShadowColor(10);
    leg_pTjb->SetTextFont(42);
    leg_pTjb->AddEntry(h_pTjb_DATA,"CMS DATA","lep");
    //leg_pTjb->AddEntry(h_pTjb_MC,"MC (Madgraph)","f");

    TLegend* leg_pTj2 = new TLegend(0.65,0.7,0.9,0.87);
    leg_pTj2->SetFillColor(0);
    leg_pTj2->SetShadowColor(10);
    leg_pTj2->SetTextFont(42);
    leg_pTj2->AddEntry(h_pTj2_DATA,"CMS DATA","lep");
    //leg_pTj2->AddEntry(h_pTj2_MC,"MC (Madgraph)","f");

    TLegend* leg_pTZ = new TLegend(0.65,0.7,0.9,0.87);
    leg_pTZ->SetFillColor(0);
    leg_pTZ->SetShadowColor(10);
    leg_pTZ->SetTextFont(42);
    leg_pTZ->AddEntry(h_pTZ_DATA,"CMS DATA","lep");
    //leg_pTZ->AddEntry(h_pTZ_MC,"MC (Madgraph)","f");

    TLegend* leg_pTj2pTZratio = new TLegend(0.65,0.7,0.9,0.87);
    leg_pTj2pTZratio->SetFillColor(0);
    leg_pTj2pTZratio->SetShadowColor(10);
    leg_pTj2pTZratio->SetTextFont(42);
    leg_pTj2pTZratio->AddEntry(h_pTj2pTZratio_DATA,"CMS DATA","lep");
    //leg_pTj2pTZratio->AddEntry(h_pTj2pTZratio_MC,"MC (Madgraph)","f");

    TCanvas *c_pTjb =new TCanvas("c_pTjb","c_pTjb",800,600);
    MakeNiceCanvas(c_pTjb);
    c_pTjb->cd();
    //h_pTjb_MC->Draw();
    h_pTjb_DATA->Draw("Pe");
    h_pTjb_MC_Stack->Draw("HISTsame"); 
    h_pTjb_DATA->Draw("Pesames");  
    theStackedLeg->Draw("same");
    //leg_pTjb->Draw("same");
    CMSPrel(ZpT_min,ZpT_max,2.1);
    DrawCanvasWithRatio(c_pTjb);

    TCanvas *c_pTj2 =new TCanvas("c_pTj2","c_pTj2",800,600);
    MakeNiceCanvas(c_pTj2);
    c_pTj2->cd();
    //h_pTj2_MC->Draw();
    h_pTj2_DATA->Draw("Pe");
    h_pTj2_MC_Stack->Draw("HISTsame");
    h_pTj2_DATA->Draw("Pesames");
    //leg_pTj2->Draw("same");
    theStackedLeg->Draw("same");
    CMSPrel(ZpT_min,ZpT_max,2.1);
    DrawCanvasWithRatio(c_pTj2);

    TCanvas *c_pTZ =new TCanvas("c_pTZ","c_pTjZ",800,600);
    MakeNiceCanvas(c_pTZ);
    c_pTZ->cd();
    //h_pTZ_MC->Draw();
    h_pTZ_DATA->Draw("Pe");
    h_pTZ_MC_Stack->Draw("HISTsame");
    h_pTZ_DATA->Draw("Pesames");
    //leg_pTZ->Draw("same");
    theStackedLeg->Draw("same");
    CMSPrel(ZpT_min,ZpT_max,2.1);  
    DrawCanvasWithRatio(c_pTZ);
      
    TCanvas *c_pTratio =new TCanvas("c_pTratio","c_pTratio",800,600);
    MakeNiceCanvas(c_pTratio);
    c_pTratio->cd();
    //h_pTj2pTZratio_MC->Draw();
    h_pTj2pTZratio_DATA->Draw("Pe");
    h_pTj2pTZratio_MC_Stack->Draw("HISTsame");
    h_pTj2pTZratio_DATA->Draw("Pesames");
    //leg_pTj2pTZratio->Draw("same");
    theStackedLeg->Draw("same");
    CMSPrel(ZpT_min,ZpT_max,2.1);
    DrawCanvasWithRatio(c_pTratio);
 
    TCanvas *c_MET =new TCanvas("c_MET","c_MET",800,600);
    MakeNiceCanvas(c_MET);
    c_MET->cd();
    //h_MET_MC->Draw();
    h_MET_DATA->Draw("Pe");
    h_MET_MC_Stack->Draw("HISTsame");
    h_MET_DATA->Draw("Pesames");
    //leg_MET->Draw("same");
    theStackedLeg->Draw("same");
    CMSPrel(ZpT_min,ZpT_max,2.1);
    DrawCanvasWithRatio(c_MET);
    
    TCanvas *c_METangle =new TCanvas("c_METangle","c_METangle",800,600);
    MakeNiceCanvas(c_METangle);
    c_METangle->cd();
    //h_METangle_MC->Draw();
    h_METangle_DATA->Draw("Pe");
    h_METangle_MC_Stack->Draw("HISTsame");
    h_METangle_DATA->Draw("Pesames");
    //leg_METangle->Draw("same");
    theStackedLeg->Draw("same");
    CMSPrel(ZpT_min,ZpT_max,2.1);
    DrawCanvasWithRatio(c_METangle);
    
  
    // Loop on variables
    for (Int_t i=0; i<n_alphaValue; ++i){

      h_Rmpf_DATA[i]->GetYaxis()->SetRangeUser(0.1,getMaximum(h_Rmpf_DATA[i],h_Rmpf_MC[i])*1.10);
      h_pTbal_DATA[i]->GetYaxis()->SetRangeUser(0.1,getMaximum(h_pTbal_DATA[i],h_pTbal_MC[i])*1.10);

      TLegend* leg_Rmpf = new TLegend(0.65,0.7,0.9,0.87);
      leg_Rmpf->SetFillColor(0);
      leg_Rmpf->SetShadowColor(10);
      leg_Rmpf->SetTextFont(42);
      leg_Rmpf->AddEntry(h_Rmpf_DATA[i],"CMS DATA","lep");
      leg_Rmpf->AddEntry(h_Rmpf_MC[i],"MC (Madgraph)","f");

      TLegend* leg_pTbal = new TLegend(0.65,0.7,0.9,0.87);
      leg_pTbal->SetFillColor(0);
      leg_pTbal->SetShadowColor(10);
      leg_pTbal->SetTextFont(42);
      leg_pTbal->AddEntry(h_pTbal_DATA[i],"CMS DATA","lep");
      leg_pTbal->AddEntry(h_pTbal_MC[i],"MC (Madgraph)","f");

      new TCanvas;
      h_Rmpf_DATA[i]->Draw("Pe");
      h_Rmpf_MC_Stack[i]->Draw("HISTsame");
      h_Rmpf_DATA[i]->Draw("Pesames");
      //leg_Rmpf->Draw("same");
      theStackedLeg->Draw("same");
      CMSPrel(ZpT_min,ZpT_max,2.1,alpha[i]);
  
      new TCanvas;
      h_pTbal_DATA[i]->Draw("Pe");
      h_pTbal_MC_Stack[i]->Draw("HISTsame");
      h_pTbal_DATA[i]->Draw("Pesames");
      //leg_pTbal->Draw("same");
      theStackedLeg->Draw("same");     
      CMSPrel(ZpT_min,ZpT_max,2.1,alpha[i]);
 
      MakeNiceHistoStyle(h_Rmpf_DATA[i]);
      for (Int_t j=0; j<n_MCSamples; ++j) {
	MakeNiceHistoStyle(h_Rmpf_MC_forStack[j][i],theStackColors[j]);
      }

      

      MakeNiceHistoStyle(h_pTbal_DATA[i]);
      for (Int_t j=0; j<n_MCSamples; ++j) {
	MakeNiceHistoStyle(h_pTbal_MC_forStack[j][i],theStackColors[j]);
      }
    }
  }
  
  //**************************
  //DATA/MC ratio
  //**************************
  //Ratio
  Float_t Ratio[n_alphaValue], Ratio_err[n_alphaValue];
  Double_t xmin_fit = 0.;
  TH1F* h_DATA[n_alphaValue]; TH1F* h_MC[n_alphaValue];

  for (Int_t i=0; i<n_alphaValue; ++i){
    //cout<<alpha_value[i]<<endl;
    if ( MPF ){
      h_DATA[i] = (TH1F*)h_Rmpf_DATA[i]->Clone("h_DATA");
      h_MC[i] = (TH1F*)h_Rmpf_MC[i]->Clone("h_MC");
    }
      
    else {
      h_DATA[i] = (TH1F*)h_pTbal_DATA[i]->Clone("h_DATA");
      h_MC[i] = (TH1F*)h_pTbal_MC[i]->Clone("h_MC");
    }
    
    Float_t mDATA, mMC, mDATA_err, mMC_err;
    if (h_DATA[i]->GetEntries()==0 || h_MC[i]->GetEntries()==0) {
      xmin_fit = alpha[i];
      Ratio[i]=0; 
      Ratio_err[i]=0;
    }
    else {
      mDATA = h_DATA[i]->GetMean();
      mMC = h_MC[i]->GetMean();
      mDATA_err = h_DATA[i]->GetMeanError();
      //cout<<"MeanError (DATA) = "<<mDATA_err<<endl;
      mMC_err = h_MC[i]->GetRMS()/(TMath::Sqrt(nEvts_MC[i]));
      //cout<<"MeanError (MC) = "<<mMC_err<<endl<<endl;
      Ratio[i] = mDATA/mMC;
      Ratio_err[i] = Ratio[i]*TMath::Sqrt(TMath::Power(mDATA_err/mDATA,2)+TMath::Power(mMC_err/mMC,2)); //with MC statistical errors included
      //Ratio_err[i] = Ratio[i]*TMath::Sqrt(TMath::Power(mDATA_err/mDATA,2)); //no MC statistical error included
    }

  }
  //for (Int_t i=0; i<n_alphaValue; ++i) cout<<endl<<alpha_value[i]<<":  DATA/MC ratio = "<<Ratio[i]<<" +- "<<Ratio_err[i]<<endl; 

  //Fit
  TF1 *sline = new TF1("sline","[0]+[1]*x");
  sline->SetLineWidth(1);
  TGraphErrors *DATAMC_ratio = new TGraphErrors(n_alphaValue, alpha, Ratio, alpha_err, Ratio_err);
  DATAMC_ratio->SetTitle("Jet Energy Response DATA/MC ratio");
  DATAMC_ratio->GetYaxis()->SetRangeUser(0.6,1.4);
  DATAMC_ratio->GetYaxis()->SetTitle("DATA/MC");
  DATAMC_ratio->GetXaxis()->SetTitle("#frac{p_{T}^{jet2}}{p_{T}^{Z}}");
  DATAMC_ratio->SetMarkerStyle(20);
  DATAMC_ratio->SetMarkerSize(1.2);
  DATAMC_ratio->SetMarkerColor(2);
  DATAMC_ratio->Fit(sline,"QC","",xmin_fit,alpha[3]);
  if (plottrend){
    new TCanvas;
    DATAMC_ratio->Draw("AP");
    CMSPrel(ZpT_min,ZpT_max,2.1);
  }

  Double_t DATAMC_ratio_extr = sline->GetParameter(0);
  Double_t DATAMC_ratio_extr_err = sline->GetParError(0);
  Double_t* ratio = new Double_t[2];

  //Standard extrapolation alpha->0
  if (ratioOpt == 1){
    ratio[0]=DATAMC_ratio_extr;
    ratio[1]=DATAMC_ratio_extr_err;
  }

  else if (ratioOpt == 2) {
    //Option for low statistics (our case)
    ratio[0] = h_DATA[3]->GetMean()/h_MC[3]->GetMean();
    ratio[1] = ratio[0]*TMath::Sqrt(pow(h_DATA[3]->GetMeanError()/h_DATA[3]->GetMean(),2) + pow(h_MC[3]->GetMeanError()/h_MC[3]->GetMean(),2));
  }

  else if (ratioOpt == 3) {
    //alpha->0/alpha=0.2 option
    ratio[0]=DATAMC_ratio_extr/(h_DATA[2]->GetMean()/h_MC[2]->GetMean());
    ratio[1]=ratio[0]*TMath::Sqrt((TMath::Power(DATAMC_ratio_extr_err,2) + TMath::Power(h_DATA[2]->GetMeanError(),2) + TMath::Power(h_MC[2]->GetMeanError(),2)));
  }

  if ( Event !=0 ) delete[] Event;
  if ( Z !=0 ) delete[] Z;
  if ( bjet != 0 ) delete[] bjet;
  if ( jet2 != 0 ) delete[] jet2;
  if ( MET != 0 ) delete[] MET;

  return ratio;
  
}

////////////////////////////////////////////////////////////////////
Bool_t isGoodForJECStudy(Z_info Z, bjet_info bjet, jet2_info jet2, Double_t ZpT_min, Double_t ZpT_max){

  Bool_t decision(false);

  Float_t deltaPhiZbjet = TMath::Abs(Z.phi_Z - bjet.phi_bjet);
  Float_t deltaRl1jet2 = TMath::Sqrt( (jet2.eta_jet2-Z.eta_d1)*(jet2.eta_jet2-Z.eta_d1) + (jet2.phi_jet2-Z.phi_d1)*(jet2.phi_jet2-Z.phi_d1) );
  Float_t deltaRl2jet2 = TMath::Sqrt( (jet2.eta_jet2-Z.eta_d2)*(jet2.eta_jet2-Z.eta_d2) + (jet2.phi_jet2-Z.phi_d2)*(jet2.phi_jet2-Z.phi_d2) );
  
  if( bjet.pT_bjet > jet2.pT_jet2 
      && TMath::Abs(bjet.eta_bjet) < 1.3 
      && deltaPhiZbjet > 2.7    
      &&  Z.pT_Z > ZpT_min 
      &&  Z.pT_Z < ZpT_max 
      ) decision = true;

  return decision;  

}

////////////////////////////////////////////////////////////////////
void MakeNiceHistoStyle(TH1F *hist, Int_t color){
  if(color!=-999){
    hist->SetFillColor(color);
    hist->SetLineColor(color);
  }  else {
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.0);
  }
  //hist->SetTitleSize(0.09); 
  hist->SetMinimum(0.01);
  hist->SetTitleSize(0.09);
  hist->SetTitleFont(42);  
  hist->SetLineWidth(2);
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetXaxis()->SetTitleOffset(1);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetTitleFont(42); 
  hist->GetYaxis()->SetTitleFont(42);  
  hist->GetXaxis()->SetLabelSize(0.06); 
  hist->GetYaxis()->SetLabelSize(0.06);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
}

////////////////////////////////////////////////////////////////////
TLegend* MakeTLegend(TH1 *h,THStack *hs,TString *LegLabels){
  TLegend* leg = new TLegend(0.65,0.75,0.85,0.93);
  //TLegend* leg = new TLegend(0.6,0.75,0.8,0.90);
  leg->SetFillColor(10);
  leg->SetLineColor(10);
  leg->SetShadowColor(10);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);  
  leg->AddEntry(h,LegLabels[0],"P");
  TIter next(hs->GetHists());
  TH1 *h1;
  int i=1;
  while ((h1 = (TH1*)next())) {
    leg->AddEntry(h1,LegLabels[i],"Lf");
    i++;
  }
  return leg;
}

////////////////////////////////////////////////////////////////////
Double_t getMaximum(TH1 *h,THStack *hs){
 
  TObjArray *array = new TObjArray();
  array->Expand(2);
  array->Add(h);
  array->Add(hs->GetHistogram());

  // TIter next(hs->GetHists());
  // TH1 *h1;
  // int i=1;
  // while ((h1 = (TH1*)next())) {
  //   leg->AddEntry(h1,LegLabels[i],"Lf");
  //   i++;
  // }

  Double_t theMaximum = (static_cast<TH1*>(array->At(0)))->GetMaximum();
  for(Int_t i = 0; i< array->GetSize(); i++){
    if( (static_cast<TH1*>(array->At(i)))->GetMaximum() > theMaximum){
      theMaximum = (static_cast<TH1*>(array->At(i)))->GetMaximum();
      //cout<<"i= "<<i<<" theMaximum="<<theMaximum<<endl;
    }
  }
  return theMaximum;
}

////////////////////////////////////////////////////////////////////
Double_t getMaximum(TH1 *h1,TH1 *h2){
 
  TObjArray *array = new TObjArray();
  array->Expand(2);
  array->Add(h1);
  array->Add(h2);

  Double_t theMaximum = (static_cast<TH1*>(array->At(0)))->GetMaximum();
  for(Int_t i = 0; i< array->GetSize(); i++){
    if( (static_cast<TH1*>(array->At(i)))->GetMaximum() > theMaximum){
      theMaximum = (static_cast<TH1*>(array->At(i)))->GetMaximum();
      //cout<<"i= "<<i<<" theMaximum="<<theMaximum<<endl;
    }
  }
  return theMaximum;
}

////////////////////////////////////////////////////////////////////
void setStyle(){ 

  //Histogram style
  new TStyle;
  //gROOT->SetStyle("Pub");
  gStyle->SetOptTitle(0); //1 if you want titles to be shown
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetPadTickX(1); 
  gStyle->SetPadTickY(1);
  gStyle->SetPadBorderMode(0); 
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetTitleOffset(1.2, "X");
  //gStyle->SetCanvasDefH(600);
  //gStyle->SetCanvasDefW(600);
}

////////////////////////////////////////////////////////////////////
void CMSPrel(Double_t ZpT_min,Double_t ZpT_max, Double_t Lumi,Double_t alpha){
  
  TLatex *latexLabel = new TLatex();
  latexLabel->SetTextSize(0.03);
  latexLabel->SetTextFont(42);
  latexLabel->SetLineWidth(2);
  latexLabel->SetNDC();
  
  Int_t type(0);

  if(alpha!=-999){
    latexLabel->DrawLatex(0.18, 0.7, "#frac{p_{T}^{jet2}}{p_{T}^{Z}} < 0."+(TString)Form("%d",(Int_t)(100*alpha)));
    type=1;
  }

  Double_t x[2]  = {0.65,0.18};
  Double_t y1[2] = {0.63,0.9};
  Double_t y2[2] = {0.57,0.84};
  Double_t y3[2] = {0.51,0.78};

  latexLabel->DrawLatex(x[type],y1[type],"anti-k_{T} (R = 0.5) PF Jets ");
  latexLabel->DrawLatex(x[type],y2[type],(TString)Form("#sqrt{s} = 7 TeV  #int Ldt = %.1f fb^{-1}",Lumi));  
  latexLabel->DrawLatex(x[type],y3[type],(TString)Form("%.f",ZpT_min)+" GeV < p_{T}^{Z} < "+(TString)Form("%.f",ZpT_max)+" GeV");  

}

////////////////////////////////////////////////////////////////////
void MakeNiceCanvas(TCanvas *canvas){
  canvas->SetFillColor(0); 
}

////////////////////////////////////////////////////////////////////
TCanvas* DrawCanvasWithRatio(TCanvas* canvas, bool emptyBinsNoUnc)
{
  // get data and total MC from the canvas
  TIter next(canvas->GetListOfPrimitives());
  TIter* mc = NULL;
  TH1F* data = NULL;
  TObject* obj = NULL;
  while ((obj = next())) {
    if(obj->InheritsFrom("TH1")) {
      data = (TH1F*)obj;
    }
    if(obj->InheritsFrom("THStack")) {
      mc = new TIter(((THStack*)obj)->GetHists());
    }
  }
  TH1F* histo_ratio = (TH1F*)(data ? data->Clone() : NULL);
  TH1F* totmc = NULL;
  if(mc) {
    while ((obj = mc->Next())) {
      if(totmc) {
        totmc->Add((TH1*)obj);
      } else {
        totmc = (TH1F*)((TH1*)obj)->Clone();
      }
    }
  }
  // if data or MC is missing, simply return the input
  if(totmc == NULL || histo_ratio == NULL) {
    return (TCanvas*) canvas->DrawClone();
  }
  // create the ratio histogram
  histo_ratio->SetName("histo_ratio");
  histo_ratio->SetTitle("");
  histo_ratio->Sumw2();
  histo_ratio->Divide(totmc);
  // set ratio to 1 for empty bins
  //for(int i=0;i<=histo_ratio->GetNbinsX();++i) {
  //  if(histo_ratio->GetBinContent(i)==0) {
  //    histo_ratio->SetBinContent(i,1);
  //    histo_ratio->SetBinError(i,0);
  //  }
  //}
  // create the uncertainty histogram
  TH1F* mc_uncertainty = (TH1F*)totmc->Clone();
  //for(unsigned bin = 0; bin<=mc_uncertainty->GetNbinsx(); ++bin) mc_uncertainty->SetBinContent(mc_uncertainty->GetBinError());
  mc_uncertainty->Divide(totmc);
  // set uncertainty to 0 or 1 for empty bins
  for(int i=0;i<=mc_uncertainty->GetNbinsX();++i) {
    if(mc_uncertainty->GetBinContent(i)==0) {
      mc_uncertainty->SetBinContent(i,1);
      mc_uncertainty->SetBinError(i,emptyBinsNoUnc? 0. : 1.);
    }
  }
  // set the color
  mc_uncertainty->SetFillColor(kYellow);
  // create a new canvas with two pads
  TCanvas* c = new TCanvas(Form("%s_withRatio",canvas->GetName()),Form("%s with ratio",canvas->GetTitle()),800,640);
  TPad *canvas_1 = new TPad("canvas_1", canvas->GetTitle(),0,0.22,1.0,1.0);
  canvas_1->Draw();
  TPad *canvas_2 = new TPad("canvas_2", Form("%s ratio",canvas->GetTitle()),0,0.,1.0,0.22);
  canvas_2->Draw();
  // in pad 1, put a copy of the input
  canvas_1->cd();
  canvas->DrawClonePad();
  // in pad 2, put the ratio plot and the relative uncertainty from MC
  canvas_2->cd();
  gPad->SetBottomMargin(0.375);
  gPad->SetGridy();
  gPad->SetGridx();
  mc_uncertainty->Draw("E3");
  mc_uncertainty->GetYaxis()->SetTitle("Data/MC");
  mc_uncertainty->GetYaxis()->SetTitleFont(42);
  mc_uncertainty->GetYaxis()->SetTitleOffset( 0.4 );
  mc_uncertainty->GetYaxis()->SetTitleSize( 0.17 );
  mc_uncertainty->GetYaxis()->SetLabelFont(42);
  mc_uncertainty->GetYaxis()->SetLabelSize(0.16);
  mc_uncertainty->GetYaxis()->SetNdivisions( 505 );
  mc_uncertainty->GetXaxis()->SetTitle(data->GetXaxis()->GetTitle());
  mc_uncertainty->GetXaxis()->SetTitleFont(42);
  mc_uncertainty->GetXaxis()->SetTitleSize( 0.17 );
  mc_uncertainty->GetXaxis()->SetLabelSize(0.16);
  mc_uncertainty->GetXaxis()->SetLabelFont(42);
  mc_uncertainty->GetXaxis()->SetRange(data->GetXaxis()->GetFirst(), data->GetXaxis()->GetLast());
  mc_uncertainty->SetMinimum(0.);
  mc_uncertainty->SetMaximum(2.);
  histo_ratio->SetMarkerStyle(20);
  histo_ratio->SetMarkerSize(0.7);
  histo_ratio->Draw("E1X0 same");
  mc_uncertainty->Draw("AXIG same");
  // return the new canvas
  return c;
}


