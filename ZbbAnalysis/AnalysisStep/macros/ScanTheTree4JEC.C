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
#include <../interface/Zbbstruct4JEC.h>
#include <Math/VectorUtil.h>
#include <fstream>
#include <iostream>

using namespace std;

Double_t* RatioExtrap(Bool_t MPF = false, Double_t ZpT_min = 10., Double_t ZpT_max = 1000., Bool_t plotspectra = false, Bool_t plottrend = false);

Bool_t isGoodForJECStudy(Z_info Z, bjet_info bjet, jet2_info jet2, Double_t ZpT_min = 10., Double_t ZpT_max = 1000.);

void ScanTheTree4JEC(Bool_t MPF = false, Bool_t integrZpT = true, Bool_t plotspectra = false, Bool_t plottrend = false, Double_t ZpT_min = 10., Double_t ZpT_max = 1000.){

  const int n_ZpT_bin(4);
  Double_t ZpT_bin_bound[n_ZpT_bin+1] = {ZpT_min,50,100,200,ZpT_max};

  if (integrZpT){
    Double_t* ratio_int = new Double_t[2];
    ratio_int = RatioExtrap(MPF, ZpT_min, ZpT_max, plotspectra, plottrend);
    cout<<"Integrating over all the range of Z pT"<<endl;
    cout<<"DATA/MC (extrapolated) = "<<ratio_int[0]<<" +- "<<ratio_int[1]<<endl;
  }
  else {
    Double_t ratio[n_ZpT_bin]; Double_t ratio_err[n_ZpT_bin];
    for (Int_t i=0; i<n_ZpT_bin; ++i){
      Double_t* ratioPluserror = RatioExtrap(MPF, ZpT_bin_bound[i], ZpT_bin_bound[i+1], false, plottrend);
      ratio[i]=ratioPluserror[0];
      ratio_err[i]=ratioPluserror[1];
    }

      TF1 *constant = new TF1("constant","[0]");
      constant->SetLineWidth(1);
      Double_t ZpT[n_ZpT_bin];
      Double_t ZpT_err[n_ZpT_bin]; 
      for (Int_t i=0; i<n_ZpT_bin; ++i) ZpT_err[i]=0;
      for (Int_t j=0; j<n_ZpT_bin; ++j) ZpT[j] = ZpT_bin_bound[j] + (ZpT_bin_bound[j+1]-ZpT_bin_bound[j])/2;
      TGraphErrors *DATAMC = new TGraphErrors(n_ZpT_bin, ZpT, ratio, ZpT_err, ratio_err);
      //DATAMC->SetTitle("Jet Energy Response DATA/MC ratio");
      DATAMC->GetYaxis()->SetRangeUser(0.6,1.4);
      DATAMC->GetYaxis()->SetTitle("DATA/MC");
      DATAMC->GetXaxis()->SetTitle("p_{T}^{Z} (GeV)");
      DATAMC->SetMarkerStyle(20);
      DATAMC->SetMarkerSize(1.2);
      DATAMC->SetMarkerColor(2);
      DATAMC->Fit(constant, "QC");
      TLegend* leg2 = new TLegend(0.65,0.7,0.9,0.87);
      leg2->SetFillColor(0);
      leg2->SetShadowColor(10);
      leg2->SetTextFont(42);
      leg2->AddEntry(constant,"Mean R_{DATA/MC}","l");
      new TCanvas;
      DATAMC->Draw("AP");
      leg2->Draw("same");
      cout<<"Constant fit to the values per-bin of ZpT"<<endl;
      cout<<"DATA/MC (extrapolated) = "<<constant->GetParameter(0)<<" +- "<<constant->GetParError(0)<<endl;

  }

}

Double_t* RatioExtrap(Bool_t MPF, Double_t ZpT_min, Double_t ZpT_max, Bool_t plotspectra, Bool_t plottrend){

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

  //Declare and book histograms
  TH1F* h_pTjb_DATA = new TH1F("h_pTjb_DATA", "b-jet p_{T}", 50, 0., 500);
  TH1F* h_pTj2_DATA = new TH1F("h_pTj2_DATA", "jet2 p_{T}", 50, 0., 500);
  TH1F* h_pTZ_DATA = new TH1F("h_pTZ_DATA", "Z p_{T}", 50, 0., 500);
  TH1F* h_pTj2pTZratio_DATA = new TH1F("h_pTj2pTZratio_DATA", "jet2 p_{T}/Z p_{T}", 80, 0., 4);
  TH1F* h_pTjb_MC = new TH1F("h_pTjb_MC", "b-jet p_{T}", 50, 0., 500);
  TH1F* h_pTj2_MC = new TH1F("h_pTj2_MC", "jet2 p_{T}", 50, 0., 500);
  TH1F* h_pTZ_MC = new TH1F("h_pTZ_MC", "Z p_{T}", 50, 0., 500);
  TH1F* h_pTj2pTZratio_MC = new TH1F("h_pTj2pTZratio_MC", "jet2 p_{T}/Z p_{T}", 80, 0., 4);

  TH1F* h_pTbal_DATA[n_alphaValue], *h_Rmpf_DATA[n_alphaValue];
  TH1F* h_pTbal_MC[n_alphaValue], *h_Rmpf_MC[n_alphaValue];
  for (Int_t i=0; i<n_alphaValue; ++i) {
    h_pTbal_DATA[i] = new TH1F("h_pTbal_DATA_"+alpha_value[i],"p_{T} b-jet/p_{T} Z for balanced events alpha="+alpha_value[i], 20, 0, 2);
    h_Rmpf_DATA[i] = new TH1F("h_Rmpf_DATA_"+alpha_value[i],"R_{MPF} spectrum for balanced events alpha="+alpha_value[i], 20, 0, 2);
    h_pTbal_MC[i] = new TH1F("h_pTbal_MC_"+alpha_value[i],"p_{T} b-jet/p_{T} Z for balanced events alpha="+alpha_value[i], 20, 0, 2);
    h_Rmpf_MC[i] = new TH1F("h_Rmpf_MC_"+alpha_value[i],"R_{MPF} spectrum for balanced events alpha="+alpha_value[i], 20, 0, 2);
  }

  //Open files
  MCfile[3] = TFile::Open("/tmp/scasasso/scripts/analyzePAT_MC_Zb5fToLL_All.root","READ");
  MCfile[0] = TFile::Open("/tmp/scasasso/scripts/analyzePAT_MC_ZcToLL_All.root","READ");
  MCfile[1] = TFile::Open("/tmp/scasasso/scripts/analyzePAT_MC_ZlJets_All.root","READ");
  MCfile[2] = TFile::Open("/tmp/scasasso/scripts/analyzePAT_MC_TTJets_All.root","READ");
  datafile = TFile::Open("/tmp/scasasso/scripts/analyzePAT_DATA2011_merged.root","READ");

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
  Double_t lumi = 2110;
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
      for (Int_t j=0; j<n_alphaValue; ++j) {
	if ( jet2[0].pT_jet2 < alpha[j]*Z[0].pT_Z  ) { // alpha selection
	  nEvts_DATA[j]++;
	  h_pTbal_DATA[j]->Fill(bjet[0].pT_bjet/Z[0].pT_Z); 
	  h_Rmpf_DATA[j]->Fill(Rmpf);
	}
      } //loop on alpha
    }//cut on Z+b candidates      			   
  }//loop on tree entries

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
	for (Int_t j=0; j<n_alphaValue; ++j) {
	  if ( jet2[k+1].pT_jet2 < alpha[j]*Z[k+1].pT_Z  ) { // alpha selection
	    nEvts_MC[j]+=w;
	    h_pTbal_MC[j]->Fill(bjet[k+1].pT_bjet/Z[k+1].pT_Z,w); 
	    h_Rmpf_MC[j]->Fill(Rmpf,w);
	  }
	} //loop on alpha
      } //cut on Z+b candidates      			   
    } //loop on tree entries
    
  }//loop on MC samples

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
  
  //Histogram style
  new TStyle;
  gROOT->SetStyle("Plain");
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

  TLatex *latexLabel = new TLatex();
  latexLabel->SetTextSize(0.03);
  latexLabel->SetTextFont(42);
  latexLabel->SetLineWidth(2);
  latexLabel->SetNDC();

  if (plotspectra){

    h_pTjb_DATA->GetXaxis()->SetTitle("p_{T} b-jet");
    h_pTjb_DATA->GetYaxis()->SetTitle("Events");
    h_pTjb_MC->GetXaxis()->SetTitle("p_{T} b-jet");
    h_pTjb_MC->GetYaxis()->SetTitle("Events");
    h_pTjb_DATA->SetMarkerStyle(20);
    h_pTjb_DATA->SetMarkerSize(1.0);
    h_pTjb_MC->SetFillColor(46);
    //h_pTjb_MC->GetYaxis()->SetRangeUser(0,maxYvalue+50);
    h_pTj2_DATA->GetXaxis()->SetTitle("p_{T} jet2");
    h_pTj2_DATA->GetYaxis()->SetTitle("Events");
    h_pTj2_MC->GetXaxis()->SetTitle("p_{T} jet2");
    h_pTj2_MC->GetYaxis()->SetTitle("Events");
    h_pTj2_DATA->SetMarkerStyle(20);
    h_pTj2_DATA->SetMarkerSize(1.0);
    h_pTj2_MC->SetFillColor(46);
    //h_pTj2_MC->GetYaxis()->SetRangeUser(0,maxYvalueJ2+50);
    h_pTZ_DATA->GetXaxis()->SetTitle("p_{T} Z");
    h_pTZ_DATA->GetYaxis()->SetTitle("Events");
    h_pTZ_MC->GetXaxis()->SetTitle("p_{T} Z");
    h_pTZ_MC->GetYaxis()->SetTitle("Events");
    h_pTZ_DATA->SetMarkerStyle(20);
    h_pTZ_DATA->SetMarkerSize(1.0);
    h_pTZ_MC->SetFillColor(46);
    //h_pTZ_MC->GetYaxis()->SetRangeUser(0,maxYvalue[i]+50);
    h_pTj2pTZratio_DATA->GetXaxis()->SetTitle("p_{T} jet2/p_{T} Z");
    h_pTj2pTZratio_DATA->GetYaxis()->SetTitle("Events");
    h_pTj2pTZratio_MC->GetXaxis()->SetTitle("p_{T} jet2/p_{T} Z");
    h_pTj2pTZratio_MC->GetYaxis()->SetTitle("Events");
    h_pTj2pTZratio_DATA->SetMarkerStyle(20);
    h_pTj2pTZratio_DATA->SetMarkerSize(1.0);
    h_pTj2pTZratio_MC->SetFillColor(46);
    //h_pTj2pTZratio_MC[i]->GetYaxis()->SetRangeUser(0,maxYvalueJ2[i]+50);
    TLegend* leg_pTjb = new TLegend(0.65,0.7,0.9,0.87);
    leg_pTjb->SetFillColor(0);
    leg_pTjb->SetShadowColor(10);
    leg_pTjb->SetTextFont(42);
    leg_pTjb->AddEntry(h_pTjb_DATA,"CMS DATA","lep");
    leg_pTjb->AddEntry(h_pTjb_MC,"MC (Madgraph)","f");
    TLegend* leg_pTj2 = new TLegend(0.65,0.7,0.9,0.87);
    leg_pTj2->SetFillColor(0);
    leg_pTj2->SetShadowColor(10);
    leg_pTj2->SetTextFont(42);
    leg_pTj2->AddEntry(h_pTj2_DATA,"CMS DATA","lep");
    leg_pTj2->AddEntry(h_pTj2_MC,"MC (Madgraph)","f");
    TLegend* leg_pTZ = new TLegend(0.65,0.7,0.9,0.87);
    leg_pTZ->SetFillColor(0);
    leg_pTZ->SetShadowColor(10);
    leg_pTZ->SetTextFont(42);
    leg_pTZ->AddEntry(h_pTZ_DATA,"CMS DATA","lep");
    leg_pTZ->AddEntry(h_pTZ_MC,"MC (Madgraph)","f");
    TLegend* leg_pTj2pTZratio = new TLegend(0.65,0.7,0.9,0.87);
    leg_pTj2pTZratio->SetFillColor(0);
    leg_pTj2pTZratio->SetShadowColor(10);
    leg_pTj2pTZratio->SetTextFont(42);
    leg_pTj2pTZratio->AddEntry(h_pTj2pTZratio_DATA,"CMS DATA","lep");
    leg_pTj2pTZratio->AddEntry(h_pTj2pTZratio_MC,"MC (Madgraph)","f");
    new TCanvas;
    h_pTjb_MC->Draw();
    h_pTjb_DATA->Draw("Pesames");
    leg_pTjb->Draw("same");
    latexLabel->DrawLatex(0.18, 0.9, "anti-k_{T} (R = 0.5) PF Jets ");
    latexLabel->DrawLatex(0.18, 0.84, "#sqrt{s} = 7 TeV  #int Ldt = 2.1 fb^{-1}");  
    latexLabel->DrawLatex(0.18, 0.78, (TString)Form("%d",(Int_t)ZpT_min)+" GeV < p_{T}^{Z} < "+(TString)Form("%d",(Int_t)ZpT_max)+" GeV");  
    new TCanvas;
    h_pTj2_MC->Draw();
    h_pTj2_DATA->Draw("Pesames");
    leg_pTj2->Draw("same");
    latexLabel->DrawLatex(0.18, 0.9, "anti-k_{T} (R = 0.5) PF Jets ");
    latexLabel->DrawLatex(0.18, 0.84, "#sqrt{s} = 7 TeV  #int Ldt = 2.1 fb^{-1}");  
    latexLabel->DrawLatex(0.18, 0.78, (TString)Form("%d",(Int_t)ZpT_min)+" GeV < p_{T}^{Z} < "+(TString)Form("%d",(Int_t)ZpT_max)+" GeV");  
    new TCanvas;
    h_pTZ_MC->Draw();
    h_pTZ_DATA->Draw("Pesames");
    leg_pTZ->Draw("same");
    latexLabel->DrawLatex(0.18, 0.9, "anti-k_{T} (R = 0.5) PF Jets ");
    latexLabel->DrawLatex(0.18, 0.84, "#sqrt{s} = 7 TeV  #int Ldt = 2.1 fb^{-1}");  
    latexLabel->DrawLatex(0.18, 0.78, (TString)Form("%d",(Int_t)ZpT_min)+" GeV < p_{T}^{Z} < "+(TString)Form("%d",(Int_t)ZpT_max)+" GeV");  
    new TCanvas;
    h_pTj2pTZratio_MC->Draw();
    h_pTj2pTZratio_DATA->Draw("Pesames");
    leg_pTj2pTZratio->Draw("same");
    latexLabel->DrawLatex(0.18, 0.9, "anti-k_{T} (R = 0.5) PF Jets ");
    latexLabel->DrawLatex(0.18, 0.84, "#sqrt{s} = 7 TeV  #int Ldt = 2.1 fb^{-1}");  
    latexLabel->DrawLatex(0.18, 0.78, (TString)Form("%d",(Int_t)ZpT_min)+" GeV < p_{T}^{Z} < "+(TString)Form("%d",(Int_t)ZpT_max)+" GeV");  
        
    for (Int_t i=0; i<n_alphaValue; ++i){
      h_Rmpf_DATA[i]->GetXaxis()->SetTitle("Z+b MPF");
      h_Rmpf_DATA[i]->GetYaxis()->SetTitle("Events");
      h_Rmpf_MC[i]->GetXaxis()->SetTitle("Z+b MPF");
      h_Rmpf_MC[i]->GetYaxis()->SetTitle("Events");
      h_Rmpf_DATA[i]->SetMarkerStyle(20);
      h_Rmpf_DATA[i]->SetMarkerSize(1.0);
      h_Rmpf_MC[i]->SetFillColor(46);
      h_Rmpf_MC[i]->GetYaxis()->SetRangeUser(0,maxYvalue[i]+50);
      h_pTbal_DATA[i]->GetXaxis()->SetTitle("p_{T} b-jet/p_{T} Z");
      h_pTbal_DATA[i]->GetYaxis()->SetTitle("Events");
      h_pTbal_MC[i]->GetXaxis()->SetTitle("p_{T} b-jet/ p_{T} Z");
      h_pTbal_MC[i]->GetYaxis()->SetTitle("Events");
      h_pTbal_DATA[i]->SetMarkerStyle(20);
      h_pTbal_DATA[i]->SetMarkerSize(1.0);
      h_pTbal_MC[i]->SetFillColor(46);
      h_pTbal_MC[i]->GetYaxis()->SetRangeUser(0,maxYvalue[i]+50);
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
      h_Rmpf_MC[i]->Draw();
      h_Rmpf_DATA[i]->Draw("Pesames");
      leg_Rmpf->Draw("same");
      latexLabel->DrawLatex(0.18, 0.9, "anti-k_{T} (R = 0.5) PF Jets ");
      latexLabel->DrawLatex(0.18, 0.84, "#sqrt{s} = 7 TeV  #int Ldt = 2.1 fb^{-1}");  
      latexLabel->DrawLatex(0.18, 0.78, (TString)Form("%d",(Int_t)ZpT_min)+" GeV < p_{T}^{Z} < "+(TString)Form("%d",(Int_t)ZpT_max)+" GeV");  
      latexLabel->DrawLatex(0.18, 0.7, "#frac{p_{T}^{jet2}}{p_{T}^{Z}} < 0."+(TString)Form("%d",(Int_t)(100*alpha[i])));
      new TCanvas;
      h_pTbal_MC[i]->Draw();
      h_pTbal_DATA[i]->Draw("Pesames");
      leg_pTbal->Draw("same");
      latexLabel->DrawLatex(0.18, 0.9, "anti-k_{T} (R = 0.5) PF Jets ");
      latexLabel->DrawLatex(0.18, 0.84, "#sqrt{s} = 7 TeV  #int Ldt = 2.1 fb^{-1}");  
      latexLabel->DrawLatex(0.18, 0.78, (TString)Form("%d",(Int_t)ZpT_min)+" GeV < p_{T}^{Z} < "+(TString)Form("%d",(Int_t)ZpT_max)+" GeV");  
      latexLabel->DrawLatex(0.18, 0.7, "#frac{p_{T}^{jet2}}{p_{T}^{Z}} < 0."+(TString)Form("%d",(Int_t)(100*alpha[i])));
    }
  }
  
  //**************************
  //DATA/MC ratio
  //**************************
  //Ratio
  Float_t Ratio[n_alphaValue], Ratio_err[n_alphaValue];
  Double_t xmin_fit = 0.;
  for (Int_t i=0; i<n_alphaValue; ++i){
    //cout<<alpha_value[i]<<endl;
    TH1F* h_DATA[n_alphaValue]; TH1F* h_MC[n_alphaValue];
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
    latexLabel->DrawLatex(0.18, 0.9, "anti-k_{T} (R = 0.5) PF Jets ");
    latexLabel->DrawLatex(0.18, 0.84, "#sqrt{s} = 7 TeV  #int Ldt = 2.1 fb^{-1}");  
    latexLabel->DrawLatex(0.18, 0.78, (TString)Form("%d",(Int_t)ZpT_min)+" GeV < p_{T}^{Z} < "+(TString)Form("%d",(Int_t)ZpT_max)+" GeV");  
  }

  Double_t DATAMC_ratio_extr = sline->GetParameter(0);
  Double_t DATAMC_ratio_extr_err = sline->GetParError(0);
  /*
  cout<<"******Fit summary*********"<<endl;
  cout<<"R_{data/mc} = "<<DATAMC_ratio_extr<<" +- "<<DATAMC_ratio_extr_err<<endl;
  cout<<"**************************"<<endl;
  */
  Double_t* ratio = new Double_t[2];
  ratio[0]=DATAMC_ratio_extr;
  ratio[1]=DATAMC_ratio_extr_err;


  if ( Event !=0 ) delete[] Event;
  if ( Z !=0 ) delete[] Z;
  if ( bjet != 0 ) delete[] bjet;
  if ( jet2 != 0 ) delete[] jet2;
  if ( MET != 0 ) delete[] MET;

  return ratio;
  
}

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
      //&& jet2.frac2b_jet2 > 0.1
      ) decision = true;

  return decision;  

}

