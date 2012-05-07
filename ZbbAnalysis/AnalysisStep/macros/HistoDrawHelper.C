#include "Riostream.h"
#include "TCanvas.h"
#include "TEnv.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>
#include <sstream>

void setStyle(){
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //  Graphic settings
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  //gStyle->SetOptStat(00000000);
  //gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFontSize(0.06);
  //gStyle->SetOptStat("emr");
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadBorderMode(0);
  //gStyle->SetNdivisions(303);
  //@@ gStyle->SetNdivisions(505);  
  //TH1::StatOverflows(kTRUE);
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.05);///---> gStyle->SetStatFontSize(0.025);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
  gStyle->SetPadBorderMode(0); 
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadBottomMargin(0.13);
  //@@ gStyle->SetPadBottomMargin(0.20);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetCanvasDefH(600); //Height of canvas
  gStyle->SetCanvasDefW(600); //Width of canvas
  // end of style

}

void MakeNiceHistoStyle(TH1 *hist, Int_t color,Bool_t StackIt, Bool_t doMCDATAratio=false){
  
  TString histName =hist->GetName();

  //hist->SetTitleSize(0.09); 
  //hist->SetMinimum(0.01);
  hist->SetTitleSize(0.09);
  hist->SetTitleFont(42); 
  //hist->Smooth();
  if(StackIt){
    hist->SetLineWidth(2);
    hist->SetLineColor(color);
    hist->SetFillColor(color);
  } else {
    hist->SetLineWidth(3);
    hist->SetLineColor(color);
  }

  if (!doMCDATAratio) {
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetTitleOffset(1);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetTitleFont(42); 
    hist->GetYaxis()->SetTitleFont(42); 

    if(histName.Contains("eventCategory")){
      hist->GetXaxis()->SetLabelSize(0.04);
    } else {
      hist->GetXaxis()->SetLabelSize(0.05); 
    }

    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
  } else {
    
    hist->GetYaxis()->SetTitleSize(25);
    hist->GetXaxis()->SetTitleSize(25);
    hist->GetXaxis()->SetTitleFont(63); 
    hist->GetYaxis()->SetTitleFont(63);  
    hist->GetXaxis()->SetLabelSize(25); 
    hist->GetYaxis()->SetLabelSize(25);
    hist->GetXaxis()->SetLabelFont(63);
    hist->GetYaxis()->SetLabelFont(63);  
    hist->GetYaxis()->SetTitleOffset(2);
    hist->GetXaxis()->SetTitleOffset(1);
   
    if(histName.Contains("eventCategory")){
      hist->GetXaxis()->SetLabelSize(20);
      hist->GetXaxis()->LabelsOption("u");
      //  Float_t x, y;
      // y = 0.1; // gPad->GetUymin() - 0.2*hist->GetYaxis()->GetBinWidth(1);
      // TText t;
      // t.SetTextAngle(60);
      // t.SetTextSize(12);
      // t.SetTextAlign(33);
      // for (Int_t i=0;i<hist->GetNbinsX();i++) {
      //	x = hist->GetXaxis()->GetBinCenter(i+1);
      //	t.DrawText(x,y,hist->GetXaxis()->GetBinLabel(i+1));
      // }
    }   
  }
}

void MakeNiceCanvas(TCanvas *canvas,Bool_t StackIt){

  canvas->SetFillColor(0); 
  TString objName =canvas->GetName();
  if(StackIt){
    //    if(!objName.Contains("Presel")){
    //       canvas->cd()->SetLogy();
    //     }
  }
  else {
    if(!objName.Contains("Tag")){
      // doMCDATAratio      if(objName.Contains("SSV")||objName.Contains("pt"))
      canvas->cd()->SetLogy();
    }
  }
  // canvas->cd()->SetGridx();
  // canvas->cd()->SetGridy(); 

}

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


TLegend* MakeTLegend(TObjArray *array,TString *LegLabels){
  TLegend* leg = new TLegend(0.67,0.65,0.87,0.95);
  //  TLegend* leg = new TLegend(0.6,0.75,0.8,0.90);
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

TString getPitchString(TH1 *h, int prec)  {
  float pitch = h->GetBinWidth(1);
  stringstream ss;
  ss << setprecision(prec);
  ss << pitch;
  TString buffer;
  ss >> buffer;
  return buffer;
}

void setAxisTitle(TH1* s, TString yTitle) {
  if (yTitle!="Arbitraty units") {
    TString pitch = getPitchString(s,2);
    yTitle = "Events/"+pitch+" "+yTitle;
  }
  s->GetYaxis()->SetTitle(yTitle.Data());
}

void setAxisTitle(THStack* s, TString yTitle) {
  TListIter hist(s->GetHists());
  TH1F* h = (TH1F*) hist.Next();
  if (yTitle!="Arbitraty units") {
    TString pitch = getPitchString(h,2);
    yTitle = "Events/"+pitch+" "+yTitle;
  }
  s->GetYaxis()->SetTitle(yTitle.Data());
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
  //histo_ratio->Sumw2();
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
  TCanvas* c = new TCanvas(Form("%s_withRatio",canvas->GetName()),Form("%s with ratio",canvas->GetTitle()),800,800);
  TPad *canvas_1 = new TPad("canvas_1", canvas->GetTitle(),0,0.19,1.0,1.0);
  canvas_1->SetFillColor(10);
  canvas_1->Draw();
  TPad *canvas_2 = new TPad("canvas_2", Form("%s ratio",canvas->GetTitle()),0,0.,1.0,0.25);
  canvas_2->SetFillColor(10);
  canvas_2->Draw();
  // in pad 1, put a copy of the input
  canvas_1->cd();
  canvas->DrawClonePad();
  // in pad 2, put the ratio plot and the relative uncertainty from MC
  canvas_2->cd();
  gPad->SetBottomMargin(0.375);
  gPad->SetGridy();
  gPad->SetGridx();
  //mc_uncertainty->SetFillStyle(3013);
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
  histo_ratio->SetMarkerStyle(33);
  //histo_ratio->SetMarkerColor(kBlue);
  histo_ratio->SetLineWidth(2);
  histo_ratio->SetMarkerSize(1.3);
  histo_ratio->Draw("E1X0 same"); 
  mc_uncertainty->Draw("AXIG same");;
  // return the new canvas
  return c;
}

void testKolmogorovAndChi2(TH1* h1, TH1* h2)
{

  TLatex *latexLabel = new TLatex();
  latexLabel->SetTextSize(0.04);
  latexLabel->SetTextFont(42);
  latexLabel->SetLineWidth(2);
  latexLabel->SetNDC();
  
  Double_t km1 =  h1->KolmogorovTest(h2,"");
  Double_t chi2 = h1->Chi2Test(h2,"UW");

  // Kolmogorov test
  //std::cout << "output of the Kolmogorov test: " <<km1<< std::endl;
  //std::cout << "output of the Kolmogorov test including normalization: " <<km2<< std::endl;
 
  latexLabel->DrawLatex(0.18,0.79,Form("k-test =%.2f #chi^{2}-test =%.2f",km1,chi2));

}
