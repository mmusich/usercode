#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TKey.h"
#include "Riostream.h"
#include <vector>
#include <sstream>
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TPaveStats.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TVectorT.h"
#include "TEnv.h"
#include "TGraphErrors.h"

Double_t getMaximum(TObjArray *array);
void MakeNiceHistoStyle(TH1 *hist, Int_t color, Bool_t StackIt);
void MakeNiceCanvas(TCanvas *canvas, Bool_t StackIt);
TLegend* MakeTLegend(TObjArray *array,TString *LegLabels);

void FastStackAndSuperImposeHistos(TString namesandlabels,Int_t nOfFiles, Bool_t StackIt){

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //  Graphic settings
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  //gStyle->SetOptStat(00000000);
  //gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(1);
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
  gStyle->SetPadTopMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetCanvasDefH(600); //Height of canvas
  gStyle->SetCanvasDefW(600); //Width of canvas

  Int_t LineColors[6]={1,6,2,3,46};
  // Normalization for: DATA/ Zbb (MC) / Zcc (MC) / TTjets (MC) / DYjtets (MC)
  //Double_t scalings[6]={1,1,1,1,1,1};
  Double_t scalings[6]={1,0.0306,0.0362,0.0278,0.0208,1};
  //Double_t scalings[6]={1,0.0295,0.03502,0.0268,0.0202,1};   // alternate scalings
  Int_t FillColors[6]={1,5,2,3,10};

  TList *FileList = new TList();
  TList *LabelList = new TList();
  
  TObjArray *nameandlabelpairs = namesandlabels.Tokenize(",");
  for (Int_t i = 0; i < nameandlabelpairs->GetEntries(); ++i) {
    TObjArray *aFileLegPair = TString(nameandlabelpairs->At(i)->GetName()).Tokenize("=");
    
    if(aFileLegPair->GetEntries() == 2) {
      FileList->Add( TFile::Open(aFileLegPair->At(0)->GetName())  );  // 2
      LabelList->Add( aFileLegPair->At(1) );
    }
    else {
      std::cout << "Please give file name and legend entry in the following form:\n" 
		<< " filename1=legendentry1,filename2=legendentry2\n";   
    }    
  }
  
  Int_t NOfFiles =  FileList->GetSize();
  if(NOfFiles!=nOfFiles){
    cout<<"&MSG-e: NOfFiles = "<<nOfFiles<<endl;
    return;
  }  

  TString LegLabels[nOfFiles];  
  
  for(Int_t j=0; j < nOfFiles; j++) {   
    TObjString* legend = (TObjString*)LabelList->At(j);
    LegLabels[j] = legend->String();
    //cout<<"LegLabels["<<j<<"]"<<LegLabels[j]<<endl;
  }
  
  TFile *file[nOfFiles];
  TObject *statObj[nOfFiles];
  TPaveStats *stats[nOfFiles];
 
  for(Int_t j=0; j < nOfFiles; j++) {  
    file[j] = (TFile*)FileList->At(j);
  }
  
  //file[0]->cd("Zplots");
  file[0]->cd("analyzeBasicPat");

  TDirectory *current_sourcedir = gDirectory;
  //gain time, do not add the objects in the list in memory
  //Bool_t status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  //loop over all keys in this directory
  //TIter nextkey( current_sourcedir->GetListOfKeys() );
  //TKey *key;
  //while ( (key = (TKey*)nextkey())) {
  
  TObject *obj;
  //TIter next(file[0]->GetListOfKeys());
  TIter next(current_sourcedir->GetListOfKeys());

  while ((obj=next())) {
    TString objName =obj->GetName();
    std::cout << "Reading: " << obj->GetName() <<std::endl;
    if (!objName.Contains("analyzeBasicPat")&& !objName.Contains("L2L3")) {

      TH1* histo[nOfFiles];
      TCanvas *theCanvas=new TCanvas();
      theCanvas->SetName(objName);
      MakeNiceCanvas(theCanvas,StackIt);
      TString theName = objName+".png";
      
      TObjArray *arrayHistos = new TObjArray();
      arrayHistos->Expand(nOfFiles);
      THStack *hs = new THStack(objName,objName);
      
      for(Int_t j=0; j < nOfFiles; j++) { 
	histo[j] = (TH1F*)file[j]->Get("analyzeBasicPat/"+objName);
	histo[j]->Sumw2();
	if(StackIt){
	  histo[j]->Scale(scalings[j]);
	  MakeNiceHistoStyle(histo[j],FillColors[j],StackIt);
	} else {
	  histo[j]->Scale(histo[0]->GetSumOfWeights()/histo[j]->GetSumOfWeights());
	  MakeNiceHistoStyle(histo[j],LineColors[j],StackIt); 
	}
	
	arrayHistos->Add(histo[j]);

	if(j!=0){
	  hs->Add(histo[j]);
	}
      }
      
      theCanvas->cd();
      histo[0]->SetMarkerStyle(20);
      histo[0]->SetMarkerSize(1.2);
      histo[0]->Draw("e");
      if(StackIt){
	hs->Draw("hits,same");
      } else {
	hs->Draw("nostack,hits,same");
      }
      histo[0]->Draw("esame");

      Double_t theMaximum = getMaximum(arrayHistos);
      histo[0]->GetYaxis()->SetRangeUser(0.001,theMaximum*1.30);
      
      theCanvas->Draw();
        
      //    for(Int_t j=0; j < nOfFiles; j++) { 
      // 	statObj[j] = histo[j]->GetListOfFunctions()->FindObject("stats");
      // 	stats[j]= static_cast<TPaveStats*>(statObj[j]);
      // 	stats[j]->SetFillColor(10);
      // 	stats[j]->SetLineWidth(1);
      // 	stats[j]->SetShadowColor(0);
      // 	stats[j]->SetTextFont(42);
      // 	stats[j]->SetLineColor(LineColors[j]);
      // 	stats[j]->SetTextColor(LineColors[j]);
      // 	stats[j]->SetX1NDC(0.76);
      // 	stats[j]->SetY1NDC(0.87-j*0.12);
      // 	stats[j]->SetX2NDC(0.97);
      // 	stats[j]->SetY2NDC(0.98-j*0.12);
      // 	stats[j]->Draw("same"); 
      //      }
      
      TLegend *legend = MakeTLegend(arrayHistos,LegLabels);
      legend->Draw("same");
      
      theCanvas->SaveAs(theName);
    }
  }
}

void MakeNiceHistoStyle(TH1 *hist, Int_t color,Bool_t StackIt){
  
  // hist->SetTitleSize(0.09); 
  //hist->SetMinimum(0.01);
  hist->SetTitleSize(0.09);
  hist->SetTitleFont(42); 
  if(StackIt){
    hist->SetLineWidth(2);
    hist->SetLineColor(1);
    hist->SetFillColor(color);
  } else {
    hist->SetLineWidth(3);
    hist->SetLineColor(color);
  }
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetXaxis()->SetTitleOffset(1);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetTitleFont(42); 
  hist->GetYaxis()->SetTitleFont(42);  
  hist->GetXaxis()->SetLabelSize(0.05); 
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
}

void MakeNiceCanvas(TCanvas *canvas,Bool_t StackIt){

  canvas->SetFillColor(0); 
  TString objName =canvas->GetName();
  if(StackIt){
    if(!objName.Contains("Presel")){
      canvas->cd()->SetLogy();
    }
  }
  else {
    if(objName.Contains("SSV")||objName.Contains("pt")){
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
  TLegend* leg = new TLegend(0.7,0.7,0.99,0.85);
  leg->SetFillColor(10);
  leg->SetShadowColor(10);
  for(Int_t j=0; j < array->GetSize(); j++){
    if(j==0){
      leg->AddEntry(static_cast<TH1*>(array->At(j)),LegLabels[j],"P");
    } else {
      leg->AddEntry(static_cast<TH1*>(array->At(j)),LegLabels[j],"f");
    }
  }
  return leg;
}
