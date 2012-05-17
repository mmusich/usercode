//
//  rooPurityEstimator.C
//
//  Created by Alberto Traverso on 15/05/12.
//  Copyright (c) 2012 Università degli Studi di Torino. All rights reserved.
//
//  Usage: 
//  root [0] .L rooPurityEstimator.C++g
//  root [1] BFractioFitter("SVmass","","ssvhpt")
//
//  arguments = "variable=SVmass,JPdisc,TCHEdisc","channel=_ee,_mm","tagger=ssvhpt,ssvhem"
//


#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TString.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooAbsData.h"
#include "RooFitResult.h" 
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"

#include <iostream>
#include "TSystem.h"
#include <Riostream.h>
#include <TROOT.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TList.h>
#include <TF1.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TVectorF.h>
#include <TKey.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <Math/VectorUtil.h>
#include <fstream>

using namespace RooFit;
using namespace std;

// graphics
void CMSPrel(Double_t Lumi);
void setStyle();
void MakeNiceHistoStyle(TH1F *hist, Int_t color=-999,Bool_t fillInside=true);
void MakeNiceGraphStyle(TGraph *hist, Int_t color);
TLegend* MakeSimpleTLegend(TObjArray *array,TString *LegLabels);
Double_t getMaximum(TObjArray *array);
TH1F* takeTheRatio(TH1F* h1,TH1F* h2);
void makeWaterMark(Double_t xNDC,Double_t yNDC, TString mytext);

// per jet purity calculator
std::pair<Float_t,Float_t> BFractionFitter(TString tagger_,TString channel_,TString variable_, TFile* fileout);

// per event purity calculator
pair<Float_t,Float_t> EventPurityCalculator(TH1F *nbjetshisto,Float_t perJetPurity,Float_t perJetPurityerr);
// function to get the calculate step by step event purity
Float_t purFunction(Float_t P,Int_t njets);

// helper variables
ofstream textfile_;
Float_t ch[3],ech[3],f_b[3],f_c[3],f_l[3],ef_b[3],ef_c[3],ef_l[3];
Int_t theStackColors[3]={628,409,856};
TString binlabels[3]={"ll","#mu#mu","ee"};

// main method
void BPurityFitter(TString tagger_="ssvhpt",TString variable_="SVmass"){

  TString fileoutname = "fits_"+tagger_+variable_+".root";
  TString textoutname = "fits_"+tagger_+variable_+".txt";
  
  TFile* fileout = new TFile(fileoutname,"RECREATE");

  textfile_.open(textoutname);
  textfile_.precision(2);

  pair<Float_t,Float_t> llJetPur = BFractionFitter(tagger_,"",variable_,fileout);
  pair<Float_t,Float_t> mmJetPur = BFractionFitter(tagger_,"_mm",variable_,fileout);
  pair<Float_t,Float_t> eeJetPur = BFractionFitter(tagger_,"_ee",variable_,fileout);

  for (UInt_t i=0; i<3; i++){
    ch[i]=i+1;
    ech[i]=0;
  }

  TGraphErrors *gr_b = new TGraphErrors(3,ch,f_b,ech,ef_b);
  gr_b->SetName("gr_b");
  TGraphErrors *gr_c = new TGraphErrors(3,ch,f_c,ech,ef_c);
  gr_c->SetName("gr_c");
  TGraphErrors *gr_l = new TGraphErrors(3,ch,f_l,ech,ef_l);
  gr_l->SetName("gr_l");

  TCanvas *c4 =new TCanvas("c4","c4",800,600);
  c4->cd();
  
  MakeNiceGraphStyle(gr_b,theStackColors[0]);
  MakeNiceGraphStyle(gr_c,theStackColors[1]);
  MakeNiceGraphStyle(gr_l,theStackColors[2]);
  
  TLegend* myleg = new TLegend(0.18,0.72,0.38,0.93);
  myleg->SetFillColor(10);
  myleg->AddEntry(gr_b,"b-jet fraction","P");
  myleg->AddEntry(gr_c,"c-jet fraction","P");
  myleg->AddEntry(gr_l,"l-jet fraction","P");

  gr_b->GetYaxis()->SetRangeUser(0.,gr_b->GetMaximum()*1.35);
  gr_b->GetYaxis()->SetTitle("flavour fraction");
  gr_b->GetXaxis()->SetTitle("channel");

  gr_b->GetHistogram()->GetXaxis()->Set(3,0.5,3.5);
  
  for(Int_t k=0;k<3;k++){
    gr_b->GetHistogram()->GetXaxis()->SetBinLabel(k+1,binlabels[k]);   
  } 
  
  c4->cd();
  gr_b->Draw("AP");
  gr_c->Draw("Psame");
  gr_l->Draw("Psame");

  myleg->Draw("same");
  CMSPrel(2.2);

  c4->Write();
  fileout->Close();

  // calculation of the event-purity

  TString filename =" templates_"+tagger_+".root";
  TFile* file = TFile::Open(filename);
  
  TH1F *nbjets_all = (TH1F*)file->Get("nbjetsdata");
  TH1F *nbjets_mm  = (TH1F*)file->Get("nbjetsdata_mm");
  TH1F *nbjets_ee  = (TH1F*)file->Get("nbjetsdata_ee");
  
  pair<Float_t,Float_t> llEventPur = EventPurityCalculator(nbjets_all,llJetPur.first,llJetPur.second);
  pair<Float_t,Float_t> mmEventPur = EventPurityCalculator(nbjets_mm ,mmJetPur.first,mmJetPur.second);
  pair<Float_t,Float_t> eeEventPur = EventPurityCalculator(nbjets_ee ,eeJetPur.first,eeJetPur.second);

  cout<<"=============== final report =============="<<endl;
  cout<<"Jet-purity ll = "<<llJetPur.first<<"+/-"<<llJetPur.second<<" Event-purity = "<< llEventPur.first <<"+/-"<<llEventPur.second<<endl;
  cout<<"Jet-purity mm = "<<mmJetPur.first<<"+/-"<<mmJetPur.second<<" Event-purity = "<< mmEventPur.first <<"+/-"<<mmEventPur.second<<endl;
  cout<<"Jet-purity ee = "<<eeJetPur.first<<"+/-"<<eeJetPur.second<<" Event-purity = "<< eeEventPur.first <<"+/-"<<eeEventPur.second<<endl;
  cout<<"==========================================="<<endl;

  textfile_<<"=============== final report =============="<<endl;
  textfile_<<"Jet-purity ll = "<<llJetPur.first<<"+/-"<<llJetPur.second<<" Event-purity = "<< llEventPur.first <<"+/-"<<llEventPur.second<<endl;
  textfile_<<"Jet-purity mm = "<<mmJetPur.first<<"+/-"<<mmJetPur.second<<" Event-purity = "<< mmEventPur.first <<"+/-"<<mmEventPur.second<<endl;
  textfile_<<"Jet-purity ee = "<<eeJetPur.first<<"+/-"<<eeJetPur.second<<" Event-purity = "<< eeEventPur.first <<"+/-"<<eeEventPur.second<<endl;
  textfile_<<"==========================================="<<endl;

  textfile_.close();

}

// Getting the file and its histograms
pair<Float_t,Float_t> BFractionFitter(TString tagger_,TString channel_,TString variable_,TFile* fileout){

  setStyle();
  
  TString n_ch="b+Z#rightarrow ll";
  if (channel_=="_mm"){
    n_ch = "b+Z#rightarrow #mu#mu";
  } else if (channel_=="_ee"){
    n_ch = "b+Z#rightarrow ee";
  }
    
  TString filename =" templates_"+tagger_+".root";
  TFile* file = TFile::Open(filename);
  TH1F *data = (TH1F*)file->Get(variable_+"data"+channel_);  // data histogram
  TH1F *mcB  = (TH1F*)file->Get(variable_+"_mc_b"+channel_); // first MC histogram
  TH1F *mcC  = (TH1F*)file->Get(variable_+"_mc_c"+channel_); // second MC histogram
  TH1F *mcL  = (TH1F*)file->Get(variable_+"_mc_l"+channel_); // third MC histogram
  
  MakeNiceHistoStyle(data,kBlack,false);
  mcB->SetFillColor(theStackColors[0]);
  mcC->SetFillColor(theStackColors[1]);
  mcL->SetFillColor(theStackColors[2]);

  //Normalized data
  //data->Sumw2();
  cout << data->Integral() << endl;
  
  RooRealVar  SVmass   = RooRealVar("SVmass","SVmass",0.,5.);
  RooDataHist DataHist = RooDataHist("DataHist","DataHist",RooArgList(SVmass),data); 
  
  //Creating signal and template PDF
  // RooDataHist SignalHist= RooDataHist("SignalHist","SignalHist",RooArgList(SVmass),data); 
  // RooHistPdf SignalPDF  = RooHistPdf("SignalPDF","SignalPDF",RooArgSet(SVmass),SignalHist); 

  RooDataHist mcBTemplateHist = RooDataHist("MCbtemplate","MCbtemplate",RooArgList(SVmass),mcB); 
  RooHistPdf MCbtemplatePDF   = RooHistPdf("MCbtemplatePDF","MCbtemplatePDF",RooArgSet(SVmass),mcBTemplateHist);
  
  RooDataHist mcCTemplateHist = RooDataHist("MCctemplate","MCctemplate",RooArgList(SVmass),mcC); 
  RooHistPdf MCctemplatePDF   = RooHistPdf("MCctemplatePDF","McctemplatePDF",RooArgSet(SVmass),mcCTemplateHist);
  
  RooDataHist mcLTemplateHist = RooDataHist("MCltemplate","MCltemplate",RooArgList(SVmass),mcL); 
  RooHistPdf MCltemplatePDF   = RooHistPdf("MCltemplatePDF","McltemplatePDF",RooArgSet(SVmass),mcLTemplateHist);
  
  RooRealVar totalB   = RooRealVar("f_{b}","fraction of bottom",0.,1.);
  RooRealVar totalC   = RooRealVar("f_{c}","fraction of charm",0.,1.);
  RooRealVar totalL   = RooRealVar("f_{l}","fraction of ligth",0.,1.);

  RooAddPdf model = RooAddPdf("model","model",RooArgList(MCbtemplatePDF,MCctemplatePDF,MCltemplatePDF),RooArgList(totalB,totalC));
  
  RooFitResult *fitres  = model.fitTo(DataHist,SumW2Error(1),RooFit::Minos(1),RooFit::Save(1),RooFit::Extended(0));
  TMatrixDSym covMatrix = fitres->covarianceMatrix();
  
  cout<<"cov[0][0]= "<<covMatrix[0][0]<<" cov[0][1]= "<<covMatrix[0][1]<<" cov[1][0]= "<<covMatrix[1][0]<<" cov[1][1]= "<<covMatrix[1][1]<<endl;

  RooPlot *xframe = SVmass.frame() ;

  xframe->SetTitle("Secondary vertex invariant mass"); 
  xframe->GetXaxis()->SetTitleSize(0.05);
  xframe->GetYaxis()->SetTitleSize(0.05);
  xframe->GetXaxis()->SetLabelSize(0.05);
  xframe->GetYaxis()->SetLabelSize(0.05);
  xframe->GetXaxis()->SetTitleOffset(0.9);
  xframe->GetYaxis()->SetTitleOffset(1.);
  xframe->GetXaxis()->SetTitleFont(42);
  xframe->GetXaxis()->SetLabelFont(42);
  xframe->GetYaxis()->SetTitleFont(42);
  xframe->GetYaxis()->SetLabelFont(42);
  xframe->GetYaxis()->SetRangeUser(1.,data->GetMaximum()*1.2);
  xframe->GetXaxis()->SetTitle("m_{SV} (GeV/c^{2})");

  DataHist.plotOn(xframe,MarkerColor(1),LineWidth(2),MarkerSize(0.9),MarkerStyle(20));
  
  //DataHist.statOn(xframe,Layout(0.63,0.93,0.7)); 
  //model.paramOn(xframe,Layout(0.63,0.93,0.5));

  //model.plotOn(xframe,Components("MCbtemplatePDF"),DrawOption("F"),LineColor(kRed),FillColor(kRed));
  //model.plotOn(xframe,Components("MCctemplatePDF"),DrawOption("F"),LineColor(kGreen),FillColor(kGreen));
  //model.plotOn(xframe,Components("MCltemplatePDF"),DrawOption("F"),LineColor(kBlue),FillColor(kBlue));

  model.plotOn(xframe,Components(RooArgSet(MCbtemplatePDF,MCctemplatePDF,MCltemplatePDF)),DrawOption("F"),LineColor(theStackColors[2]),FillColor(theStackColors[2])); 
  model.plotOn(xframe,Components(RooArgSet(MCbtemplatePDF,MCctemplatePDF)),DrawOption("F"),LineColor(theStackColors[1]),FillColor(theStackColors[1]));
  model.plotOn(xframe,Components(RooArgSet(MCbtemplatePDF)),DrawOption("F"),LineColor(theStackColors[0]),FillColor(theStackColors[0]));

  //model.plotOn(xframe,LineColor(1));
  DataHist.plotOn(xframe,MarkerColor(1),LineWidth(2),MarkerSize(0.9),MarkerStyle(20));

  TString output = Form("#chi^{2}=%5.3f",xframe->chiSquare());
  
  TCanvas* c = new TCanvas("finalfit"+tagger_+variable_+channel_,"finalfit"+tagger_+variable_+channel_,800,600) ;
  c->cd(); 
  gPad->SetLeftMargin(0.12);

  xframe->Draw();
  CMSPrel(2.2);
  makeWaterMark(0.14,0.91,n_ch);
  makeWaterMark(0.14,0.86,"after fit");
 
  Float_t fb    = totalB.getVal();
  Float_t fberr = totalB.getError();
  
  Float_t fc    = totalC.getVal();
  Float_t fcerr = totalC.getError();

  Float_t fl    = 1-fb-fc;
  Float_t flerr = TMath::Sqrt(TMath::Power(fberr,2)+TMath::Power(fcerr,2)); //+2*covMatrix[1][0]);

  TLegend *leg = new TLegend(0.65,0.47,0.92,0.74);
  leg->SetTextFont(42);
  leg->AddEntry(data,"CMS data","p");
  leg->AddEntry(mcB,Form("f_{b} = %5.3f #pm %5.3f",fb,fberr),"f");
  leg->AddEntry(mcC,Form("f_{c} = %5.3f #pm %5.3f",fc,fcerr),"f");
  leg->AddEntry(mcL,Form("f_{l} = %5.3f #pm %5.3f",fl,flerr),"f");
  leg->SetFillColor(10);
  leg->SetLineColor(10);
  leg->Draw(); 

  textfile_<<"=================== "<<n_ch<<" ================="<<endl;
  textfile_<<Form("f_{b} = %5.3f +/- %5.3f",fb,fberr)<<endl;
  textfile_<<Form("f_{c} = %5.3f +/- %5.3f",fc,fberr)<<endl;
  textfile_<<Form("f_{l} = %5.3f +/- %5.3f",fl,flerr)<<endl;
  textfile_<<endl;
  textfile_<<output<<endl;
  
  if(channel_=="") {
    f_b[0]=fb;
    f_c[0]=fc; 
    f_l[0]=fl;
    ef_b[0]=fberr;
    ef_c[0]=fcerr; 
    ef_l[0]=flerr;
  } else if (channel_=="_mm"){
    f_b[1]=fb;
    f_c[1]=fc; 
    f_l[1]=fl;
    ef_b[1]=fberr;
    ef_c[1]=fcerr; 
    ef_l[1]=flerr;
  } else {
    f_b[2]=fb;
    f_c[2]=fc; 
    f_l[2]=fl;
    ef_b[2]=fberr;
    ef_c[2]=fcerr; 
    ef_l[2]=flerr;
  }

  TLatex *ll = new TLatex();
  ll->SetTextSize(0.04);
  ll->SetTextFont(42);
  ll->SetLineWidth(2);
  ll->SetNDC();
  ll->DrawLatex(0.65,0.40,output);
  
  c->SaveAs("finalfit"+tagger_+variable_+channel_+".png");
  c->SaveAs("finalfit"+tagger_+variable_+channel_+".pdf");

  totalB.Print();
  totalC.Print();
  //totalL.Print();

  c->Draw(); 

  fileout->cd();
  c->Write();

  pair<Float_t,Float_t> result;
  result = make_pair(totalB.getVal(),totalB.getError());
  
  return result;
  
}

////////////////////////////////////////////////////////////////////
void setStyle(){ 

  //Histogram style
  new TStyle;
  //gROOT->SetStyle("Pub");
  gStyle->SetOptTitle(0); //1 if you want titles to be shown
  gStyle->SetTitleFont(42);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetCanvasColor(kWhite);
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
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.03);
  //gStyle->SetTitleOffset(1.2, "X");
  //gStyle->SetCanvasDefH(600);
  //gStyle->SetCanvasDefW(600);
}

////////////////////////////////////////////////////////////////////
void MakeNiceHistoStyle(TH1F *hist, Int_t color,Bool_t fillInside){
  if(color!=-999){
    if(fillInside){
      hist->SetFillColor(color);
    } else {
      hist->SetLineColor(color);
      hist->SetLineWidth(2);
      hist->SetMarkerStyle(20);
      hist->SetMarkerSize(0.9);
    }
  }  
  //hist->SetTitleSize(0.09); 
  hist->SetMinimum(0.01);
  hist->SetTitleSize(0.09);
  hist->SetTitleFont(42);  
  hist->SetLineWidth(2);
  hist->GetYaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleOffset(1.1);
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
void MakeNiceGraphStyle(TGraph *graph, Int_t color){

  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerColor(color);

  graph->SetMinimum(0.01);
  graph->SetFillColor(color);
  graph->SetLineColor(color);
  graph->SetLineWidth(2);
  graph->GetYaxis()->SetTitleOffset(1.3);
  graph->GetXaxis()->SetTitleOffset(1.1);
  graph->GetYaxis()->SetTitleSize(0.06);
  graph->GetXaxis()->SetTitleSize(0.06);
  graph->GetXaxis()->SetTitleFont(42); 
  graph->GetYaxis()->SetTitleFont(42);  
  graph->GetXaxis()->SetLabelSize(0.07); 
  graph->GetYaxis()->SetLabelSize(0.06);
  graph->GetXaxis()->SetLabelFont(42);
  graph->GetYaxis()->SetLabelFont(42);
}


////////////////////////////////////////////////////////////////////
void CMSPrel(Double_t Lumi){
  
  TLatex *latexLabel = new TLatex();
  latexLabel->SetTextSize(0.035);
  latexLabel->SetTextFont(42);
  latexLabel->SetLineWidth(2);
  latexLabel->SetNDC();

  latexLabel->DrawLatex(0.66,0.9,"CMS VHF Preliminary");
  latexLabel->DrawLatex(0.66,0.84,"anti-k_{T} (R = 0.5) PF Jets ");
  latexLabel->DrawLatex(0.66,0.78,(TString)Form("#sqrt{s} = 7 TeV  #int Ldt = %.1f fb^{-1}",Lumi));  

}

////////////////////////////////////////////////////////////////////
TLegend* MakeSimpleTLegend(TObjArray *array,TString *LegLabels){
  TLegend* leg = new TLegend(0.60,0.53,0.80,0.74);
  leg->SetFillColor(10);
  leg->SetLineColor(10);
  leg->SetShadowColor(10);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);  
  for(Int_t j=0; j < array->GetSize(); j++){
    if(j==0){
      leg->AddEntry(static_cast<TH1*>(array->At(j)),LegLabels[j],"P");
    } else {
      if(static_cast<TH1*>(array->At(j))->GetEntries()>0)
	leg->AddEntry(static_cast<TH1*>(array->At(j)),LegLabels[j],"f");
    }
  }
  return leg;
}

////////////////////////////////////////////////////////////////////
Double_t getMaximum(TObjArray *array){

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
TH1F* takeTheRatio(TH1F* h1,TH1F* h2){

  TH1F* histo_ratio = (TH1F*)h2 ? (TH1F*)h1->Clone() : NULL;
  histo_ratio->SetName("histo_ratio");
  histo_ratio->Divide(h2);
  histo_ratio->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
  histo_ratio->GetYaxis()->SetTitle("Ratio");
  histo_ratio->SetMinimum(0.5);
  histo_ratio->SetMaximum(1.5);

  histo_ratio->GetYaxis()->SetTitleOffset(0.5);
  histo_ratio->GetYaxis()->SetTitleFont(42);
  histo_ratio->GetYaxis()->SetTitleSize(0.15);
  histo_ratio->GetYaxis()->SetLabelFont(42);
  histo_ratio->GetYaxis()->SetLabelSize(0.15);
  histo_ratio->GetYaxis()->SetNdivisions( 505 );

  histo_ratio->GetXaxis()->SetTitleFont(42);
  histo_ratio->GetXaxis()->SetTitleSize(0.15);
  histo_ratio->GetXaxis()->SetLabelFont(42);
  histo_ratio->GetXaxis()->SetLabelSize(0.15);
  
  histo_ratio->SetLineWidth(1.);
  histo_ratio->SetMarkerStyle(20);
  histo_ratio->SetMarkerColor(kRed);
  histo_ratio->SetMarkerSize(1.);

  return histo_ratio;

}

////////////////////////////////////////////////////////////////////
void makeWaterMark(Double_t xNDC,Double_t yNDC, TString mytext){
  
  TLatex *ll = new TLatex();
  ll->SetTextSize(0.04);
  ll->SetTextFont(42);
  ll->SetLineWidth(2);
  //ll->SetTextColor(kRed);
  ll->SetNDC();
  ll->DrawLatex(xNDC,yNDC,mytext);

}

////////////////////////////////////////////////////////////////////
pair<Float_t,Float_t> EventPurityCalculator(TH1F *nbjetshisto,Float_t perJetPurity,Float_t perJetPurityerr){

  cout.precision(3);

  Float_t P_b = perJetPurity;
  vector<pair<Int_t,Float_t> > mult_and_frac;
  
  for(Int_t n=0; n<nbjetshisto->GetNbinsX(); n++){
    if(nbjetshisto->GetBinContent(n)>0.){ 
      pair<Int_t,Float_t> bin_and_frac;
      bin_and_frac.first = n-1;
      bin_and_frac.second = nbjetshisto->GetBinContent(n)/nbjetshisto->GetSumOfWeights();
      mult_and_frac.push_back(bin_and_frac);
    }
  }

  Float_t perEventPurity(0.);
  Float_t diff(0);
  Float_t perEventPurityError(0.);

  for(UInt_t i=0; i<mult_and_frac.size(); i++){
    perEventPurity+=(mult_and_frac[i].second)*purFunction(P_b,mult_and_frac[i].first);
    diff+=(mult_and_frac[i].second)*(mult_and_frac[i].first)*TMath::Power((1.-P_b),(mult_and_frac[i].first-1));
    cout<<"bin= "<<fixed<< mult_and_frac[i].first<<" f("<<mult_and_frac[i].first<<"b)="<<fixed<<mult_and_frac[i].second<<
      " P("<<mult_and_frac[i].first<<"b)="<<fixed<<purFunction(P_b,mult_and_frac[i].first)<<" Pur="<<fixed<<perEventPurity<<endl;
  }

  cout<<"........."<<endl;

  perEventPurityError=TMath::Sqrt(TMath::Power(diff,2)*TMath::Power(perJetPurityerr,2));

  pair<Float_t,Float_t> result;
  result = make_pair(perEventPurity,perEventPurityError);

  return result;

}

////////////////////////////////////////////////////////////////////
Float_t purFunction(Float_t P,Int_t njets){
  
  Float_t P_e = 1.-TMath::Power((1.-P),njets);
  return P_e;
}


 











