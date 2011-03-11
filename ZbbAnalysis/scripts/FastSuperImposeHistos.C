////////////////////////////////////////////
//
// Simple macro to superimpose histograms
// from files having the same structure
//
// Original author: M. Musich INFN Torino
//
////////////////////////////////////////////

#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
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
void MakeNiceHistoStyle(TH1 *hist, Int_t color);
void MakeNiceCanvas(TCanvas *canvas);
TLegend* MakeTLegend(TObjArray *array,TString *LegLabels);

void FastSuperImposeHistos(TString namesandlabels,Int_t nOfFiles){

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
  gStyle->SetOptStat("emr");
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
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.04);

  Int_t LineColors[6]={1,2,4,6,7,3};
  // Int_t FillColors[6]={1,2,4,6,7,3};

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
  
  file[0]->cd("ZPlots");

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
    if (!objName.Contains("ZPlots")) {

      TH1* histo[nOfFiles];
      TCanvas *theCanvas=new TCanvas();
      theCanvas->SetName(objName);
      MakeNiceCanvas(theCanvas);
      TString theName = objName+".eps";
      
      TObjArray *arrayHistos = new TObjArray();
      arrayHistos->Expand(nOfFiles);
      for(Int_t j=0; j < nOfFiles; j++) { 
	histo[j] = (TH1F*)file[j]->Get("ZPlots/"+objName);
	arrayHistos->Add(histo[j]);
	MakeNiceHistoStyle(histo[j],LineColors[j]);
	theCanvas->cd();
	if(j==0){
	  histo[j]->Draw();
	} else{
	  histo[j]->Draw("sames");
	}
      }
      Double_t theMaximum = getMaximum(arrayHistos);
      histo[0]->GetYaxis()->SetRangeUser(0.01,theMaximum*1.10);
      
      theCanvas->Draw();
        
      for(Int_t j=0; j < nOfFiles; j++) { 
	statObj[j] = histo[j]->GetListOfFunctions()->FindObject("stats");
	stats[j]= static_cast<TPaveStats*>(statObj[j]);
	stats[j]->SetFillColor(10);
	stats[j]->SetLineWidth(1);
	stats[j]->SetShadowColor(0);
	stats[j]->SetTextFont(42);
	stats[j]->SetLineColor(LineColors[j]);
	stats[j]->SetTextColor(LineColors[j]);
	stats[j]->SetX1NDC(0.76);
	stats[j]->SetY1NDC(0.87-j*0.12);
	stats[j]->SetX2NDC(0.97);
	stats[j]->SetY2NDC(0.98-j*0.12);
	stats[j]->Draw("same"); 
      }
      
      TLegend *legend = MakeTLegend(arrayHistos,LegLabels);
      legend->Draw("same");

      theCanvas->SaveAs(theName);
    }
  }
}

void MakeNiceHistoStyle(TH1 *hist, Int_t color){
  
  // hist->SetTitleSize(0.09); 
  hist->SetMinimum(0.01);
  hist->SetTitleSize(0.09);
  hist->SetTitleFont(42);  
  hist->SetLineWidth(2);
  hist->SetLineColor(color);
  //hist->SetFillColor(color);
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

void MakeNiceCanvas(TCanvas *canvas){

  canvas->SetFillColor(0); 
  TString objName =canvas->GetName();
  if(objName.Contains("SSV")){
     canvas->cd()->SetLogy();
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
  TLegend *leg = new TLegend(0.54,0.90,0.75,0.99);
  leg->SetFillColor(10);
  leg->SetShadowColor(10);
  for(Int_t j=0; j < array->GetSize(); j++){
    leg->AddEntry(static_cast<TH1*>(array->At(j)),LegLabels[j]);
  }
  return leg;
}
