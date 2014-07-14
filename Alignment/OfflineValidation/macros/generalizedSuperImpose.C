#include "TKey.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TList.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TObject.h"
#include "TROOT.h"
#include "TStyle.h"
#include <Riostream.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>
#include <TText.h>
#include <iomanip>
#include <iostream>
#include <stdio.h>

void makeNiceCanv(TCanvas *c);
std::pair<Double_t,Double_t> getExtrema(TObjArray *array);
void makeNicePlotStyleAndColor(TH1F *hist,int color, std::pair<Double_t,Double_t> the_extrema);

void generalizedSuperImpose(TString namesandlabels,Int_t nOfFiles){

  gStyle->SetOptTitle(1);
  gStyle->SetTitleFont(62);
  //gStyle->SetTitleColor(4);
  gStyle->SetTitleTextColor(4);
  gStyle->SetTitleFillColor(10);

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
  
  const Int_t NOfFiles = FileList->GetSize();
  if(NOfFiles!=nOfFiles){
    cout<<"&MSG-e: NOfFiles = "<<nOfFiles<<endl;
    return;
  }  
  
  TString legs[10];  
  
  for(Int_t j=0; j < NOfFiles; j++) {   
    TObjString* legend = (TObjString*)LabelList->At(j);
    legs[j] = legend->String();
    cout<<"LegLabels["<<j<<"]: "<<legs[j]<<endl;
  }
  
  TFile *file[NOfFiles];
  
  for(Int_t j=0; j < NOfFiles; j++) { 
    file[j] = (TFile*)FileList->At(j);
  }
 
  const char* folders[2] = {"Pulls","rmsDMR"};
  const char* vars[3]    = {"Eta","Phi","Z"};  

  TCanvas* cBPixX[2][3];
  TCanvas* cFPixX[2][3];
  	 	   
  TCanvas* cBPixY[2][3];
  TCanvas* cFPixY[2][3];
  	 	   
  TCanvas* cTIB[2][3];  
  TCanvas* cTID[2][3];  
  TCanvas* cTOB[2][3]; 
  TCanvas* cTEC[2][3];  

  for (Int_t kF=0; kF < 2; kF++) { 
    for (Int_t kV=0; kV < 3; kV++) { 

      cBPixX[kF][kV] = new TCanvas(Form("cBPixX_%s_Vs_%s",folders[kF],vars[kV]),Form("BPix-x %s vs %s",folders[kF],vars[kV]),700,700);
      cFPixX[kF][kV] = new TCanvas(Form("cFPixX_%s_Vs_%s",folders[kF],vars[kV]),Form("FPix-x %s vs %s",folders[kF],vars[kV]),700,700);
      
      cBPixY[kF][kV] = new TCanvas(Form("cBPixY_%s_Vs_%s",folders[kF],vars[kV]),Form("BPix-y %s vs %s",folders[kF],vars[kV]),700,700);
      cFPixY[kF][kV] = new TCanvas(Form("cFPixY_%s_Vs_%s",folders[kF],vars[kV]),Form("FPix-y %s vs %s",folders[kF],vars[kV]),700,700);
      
      cTIB[kF][kV]   = new TCanvas(Form("cTIB_%s_Vs_%s",folders[kF],vars[kV]),Form("TIB %s vs %s",folders[kF],vars[kV]),700,700);
      cTID[kF][kV]   = new TCanvas(Form("cTID_%s_Vs_%s",folders[kF],vars[kV]),Form("TID %s vs %s",folders[kF],vars[kV]),700,700);
      cTOB[kF][kV]   = new TCanvas(Form("cTOB_%s_Vs_%s",folders[kF],vars[kV]),Form("TOB %s vs %s",folders[kF],vars[kV]),700,700);
      cTEC[kF][kV]   = new TCanvas(Form("cTEC_%s_Vs_%s",folders[kF],vars[kV]),Form("TEC %s vs %s",folders[kF],vars[kV]),700,700);
    
      makeNiceCanv(cBPixX[kF][kV]);
      makeNiceCanv(cFPixX[kF][kV]);
      makeNiceCanv(cBPixY[kF][kV]);
      makeNiceCanv(cFPixY[kF][kV]);
      
      makeNiceCanv(cTIB[kF][kV]);
      makeNiceCanv(cTID[kF][kV]);
      makeNiceCanv(cTOB[kF][kV]);
      makeNiceCanv(cTEC[kF][kV]);
    
    }
  }

  TH1F* h1xBPix[nOfFiles][2][3];
  TH1F* h1yBPix[nOfFiles][2][3];
  TH1F* h1xFPix[nOfFiles][2][3];
  TH1F* h1yFPix[nOfFiles][2][3];
  TH1F* h1TIB[nOfFiles][2][3];
  TH1F* h1TID[nOfFiles][2][3];
  TH1F* h1TOB[nOfFiles][2][3];
  TH1F* h1TEC[nOfFiles][2][3]; 

  for (Int_t kF=0; kF < 2; kF++) { 
    
    file[0]->cd(folders[kF]);
    
    TDirectory *currentDir = gDirectory;
    TKey *key;
    TIter nextkey(currentDir->GetListOfKeys());
    
    TObject *obj;

    while ((key = (TKey*)nextkey())) {
      obj = key->ReadObj();
      //std::cout<<"object name is:"<<TString((obj->GetName()))<<std::endl;
      TString objName = TString((obj->GetName()));
    
      Int_t _tI=0;
      if(objName.Contains("Eta")){
	_tI = 0;
      } else if(objName.Contains("Phi")){
	_tI = 1;
      } else if(objName.Contains("Z")){
	_tI = 2;
      } else {
	_tI = -1;
	std::cout<<objName<<"unknown variable. This is gonna break!"<<std::endl;
      }

      for(Int_t j=0; j < NOfFiles; j++) {     
	
	TCanvas *c1 =(TCanvas*)file[j]->Get(TString(folders[kF])+"/"+objName);
	
	TObject *oo1;
	TIter next1(c1->GetListOfPrimitives());
	
	while ((oo1=next1())) {
	  if(oo1->InheritsFrom("TH1")){
	    //std::cout<<"h1 name:"<<oo1->GetName()<<std::endl;
	    TString oo1Name = oo1->GetName();
	    if(oo1Name.Contains("BPix-x")){
	      h1xBPix[j][kF][_tI]=(TH1F*)oo1;
	    } else if(oo1Name.Contains("BPix-y")){
	      h1yBPix[j][kF][_tI]=(TH1F*)oo1;
	    } else if(oo1Name.Contains("FPix-x")){
	      h1xFPix[j][kF][_tI]=(TH1F*)oo1;
	    } else if(oo1Name.Contains("FPix-y")){
	      h1yFPix[j][kF][_tI]=(TH1F*)oo1;
	    } else if(oo1Name.Contains("TIB")){
	      h1TIB[j][kF][_tI]=(TH1F*)oo1;
	    } else if(oo1Name.Contains("TID")){
	      h1TID[j][kF][_tI]=(TH1F*)oo1;
	    } else if(oo1Name.Contains("TOB")){
	      h1TOB[j][kF][_tI]=(TH1F*)oo1;
	    } else if(oo1Name.Contains("TEC")){
	      h1TEC[j][kF][_tI]=(TH1F*)oo1;
	    }
	  }
	} 
	//delete oo1;
	//delete c1;
      }
    }
  }

  TLegend *leg = new TLegend(0.15,0.75,0.50,0.90);
  //leg->SetNColumns(2);
  leg->SetLineColor(10);
  leg->SetFillColor(10);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  
  // Get all the extrema
  
  std::pair<Double_t,Double_t> the_extremaxBPix[2][3];
  std::pair<Double_t,Double_t> the_extremayBPix[2][3];
  std::pair<Double_t,Double_t> the_extremaxFPix[2][3];
  std::pair<Double_t,Double_t> the_extremayFPix[2][3];
  std::pair<Double_t,Double_t> the_extremaTIB[2][3];
  std::pair<Double_t,Double_t> the_extremaTOB[2][3];
  std::pair<Double_t,Double_t> the_extremaTID[2][3];
  std::pair<Double_t,Double_t> the_extremaTEC[2][3];   

  for (Int_t kF=0; kF < 2; kF++) { 
    for (Int_t kV=0; kV < 3; kV++) {

      TObjArray *arrxBPix = new TObjArray();
      TObjArray *arryBPix = new TObjArray();
      TObjArray *arrxFPix = new TObjArray();
      TObjArray *arryFPix = new TObjArray();
      TObjArray *arrTIB   = new TObjArray();
      TObjArray *arrTOB   = new TObjArray();
      TObjArray *arrTID   = new TObjArray();
      TObjArray *arrTEC   = new TObjArray();

      arrxBPix->Expand(nOfFiles);
      arryBPix->Expand(nOfFiles);
      arrxFPix->Expand(nOfFiles);
      arryFPix->Expand(nOfFiles);
      arrTIB->Expand(nOfFiles);  
      arrTOB->Expand(nOfFiles);  
      arrTID->Expand(nOfFiles);  
      arrTEC->Expand(nOfFiles);  

      for(Int_t i =0;i<nOfFiles;i++){
	arrxBPix->Add(h1xBPix[i][kF][kV]);
	arryBPix->Add(h1yBPix[i][kF][kV]);
	arrxFPix->Add(h1xFPix[i][kF][kV]);
	arryFPix->Add(h1yFPix[i][kF][kV]);
	arrTIB->Add(h1TIB[i][kF][kV]);  
	arrTOB->Add(h1TOB[i][kF][kV]);  
	arrTID->Add(h1TID[i][kF][kV]);  
	arrTEC->Add(h1TEC[i][kF][kV]);  
      }

      the_extremaxBPix[kF][kV] = getExtrema(arrxBPix);
      the_extremayBPix[kF][kV] = getExtrema(arryBPix);
      the_extremaxFPix[kF][kV] = getExtrema(arrxFPix);
      the_extremayFPix[kF][kV] = getExtrema(arryFPix);
      the_extremaTIB[kF][kV]   = getExtrema(arrTIB);
      the_extremaTOB[kF][kV]   = getExtrema(arrTOB);
      the_extremaTID[kF][kV]   = getExtrema(arrTID);
      the_extremaTEC[kF][kV]   = getExtrema(arrTEC);
    }
  }

  
  for (Int_t kF=0; kF < 2; kF++) { 
    for (Int_t kV=0; kV < 3; kV++) { 

      //%%%%%%%%%%%%%%%%%%%%%%%%%%
      cBPixX[kF][kV]->cd(); 
      for(Int_t i =0;i<nOfFiles;i++){
	//  std::cout<<"i:"<<i<<std::endl;
	makeNicePlotStyleAndColor(h1xBPix[i][kF][kV],i,the_extremaxBPix[kF][kV]);
	//std::cout<<"i:"<<i<<std::endl;
	if(i==0) h1xBPix[i][kF][kV]->Draw("");
	else h1xBPix[i][kF][kV]->Draw("same");
	if(kF==0 && kV==0) leg->AddEntry(h1xBPix[i][kF][kV],legs[i]);
      }
      leg->Draw("same");
      cBPixX[kF][kV]->SaveAs(Form("cBPixX_%s_Vs_%s.pdf",folders[kF],vars[kV]));
      
      //%%%%%%%%%%%%%%%%%%%%%%%%%%
      cFPixX[kF][kV]->cd();
      for(Int_t i =0;i<nOfFiles;i++){
	makeNicePlotStyleAndColor(h1xFPix[i][kF][kV],i,the_extremaxFPix[kF][kV]);
	if(i==0)h1xFPix[i][kF][kV]->Draw("");
	else h1xFPix[i][kF][kV]->Draw("same");
      }
      leg->Draw("same");
      cFPixX[kF][kV]->SaveAs(Form("cFPixX_%s_Vs_%s.pdf",folders[kF],vars[kV]));
      
      //%%%%%%%%%%%%%%%%%%%%%%%%%%
      cBPixY[kF][kV]->cd(); 
      for(Int_t i =0;i<nOfFiles;i++){
	makeNicePlotStyleAndColor(h1yBPix[i][kF][kV],i,the_extremayBPix[kF][kV]);
	if(i==0) h1yBPix[i][kF][kV]->Draw("");
	else h1yBPix[i][kF][kV]->Draw("same");
      }
      leg->Draw("same");
      cBPixY[kF][kV]->SaveAs(Form("cBPixY_%s_Vs_%s.pdf",folders[kF],vars[kV]));
      
      //%%%%%%%%%%%%%%%%%%%%%%%%%%
      cFPixY[kF][kV]->cd();
      for(Int_t i =0;i<nOfFiles;i++){
	makeNicePlotStyleAndColor(h1yFPix[i][kF][kV],i,the_extremayFPix[kF][kV]);
	if(i==0) h1yFPix[i][kF][kV]->Draw("");
	else h1yFPix[i][kF][kV]->Draw("same");
      }
      leg->Draw("same");
      cFPixY[kF][kV]->SaveAs(Form("cFPixY_%s_Vs_%s.pdf",folders[kF],vars[kV]));
      
      //%%%%%%%%%%%%%%%%%%%%%%%%%%
      cTIB[kF][kV]->cd(); 
      for(Int_t i =0;i<nOfFiles;i++){
	makeNicePlotStyleAndColor(h1TIB[i][kF][kV],i,the_extremaTIB[kF][kV]);
	if(i==0) h1TIB[i][kF][kV]->Draw("");
	else h1TIB[i][kF][kV]->Draw("same");
      }
      leg->Draw("same");
      cTIB[kF][kV]->SaveAs(Form("cTIB_%s_Vs_%s.pdf",folders[kF],vars[kV]));
      
      //%%%%%%%%%%%%%%%%%%%%%%%%%%
      cTID[kF][kV]->cd(); 
      for(Int_t i =0;i<nOfFiles;i++){
	makeNicePlotStyleAndColor(h1TID[i][kF][kV],i,the_extremaTID[kF][kV]);
	if(i==0) h1TID[i][kF][kV]->Draw("");
	else h1TID[i][kF][kV]->Draw("same");
      }
      leg->Draw("same");
      cTID[kF][kV]->SaveAs(Form("cTID_%s_Vs_%s.pdf",folders[kF],vars[kV]));
      
      //%%%%%%%%%%%%%%%%%%%%%%%%%%
      cTOB[kF][kV]->cd();  
      for(Int_t i =0;i<nOfFiles;i++){
	makeNicePlotStyleAndColor(h1TOB[i][kF][kV],i,the_extremaTOB[kF][kV]);
	if(i==0) h1TOB[i][kF][kV]->Draw("");
	else h1TOB[i][kF][kV]->Draw("same");
      }
      leg->Draw("same");
      cTOB[kF][kV]->SaveAs(Form("cTOB_%s_Vs_%s.pdf",folders[kF],vars[kV]));
      
      //%%%%%%%%%%%%%%%%%%%%%%%%%%
      cTEC[kF][kV]->cd(); 
      for(Int_t i =0;i<nOfFiles;i++){
	makeNicePlotStyleAndColor(h1TEC[i][kF][kV],i,the_extremaTEC[kF][kV]);
	if(i==0) h1TEC[i][kF][kV]->Draw("");
	else h1TEC[i][kF][kV]->Draw("same");
      }
      leg->Draw("same");
      cTEC[kF][kV]->SaveAs(Form("cTEC_%s_Vs_%s.pdf",folders[kF],vars[kV]));
    }
  }
}

/*--------------------------------------------------------------------*/
void makeNicePlotStyleAndColor(TH1F *hist,int color,std::pair<Double_t,Double_t> the_extrema){
/*--------------------------------------------------------------------*/

  //Int_t styles[4] = {kFullCircle,kFullSquare,kFullTriangleUp,kFullTriangleDown};
  Int_t styles[4] = {kFullCircle,kFullSquare,kFullCircle,kFullSquare};
  Int_t colors[4] = {kRed,kBlue,kMagenta,kMagenta-2};

  hist->SetStats(false);
  //hist->GetYaxis()->SetRangeUser(the_extrema.first,the_extrema.second*1.5);
  //hist->GetYaxis()->SetRangeUser(0.1,the_extrema.second*1.5);
  //std::cout<<"the_extrema.second: "<<the_extrema.second<<std::endl;
  hist->GetYaxis()->SetRangeUser(0.01,the_extrema.second*1.6);
  //hist->GetYaxis()->UnZoom();
  hist->GetXaxis()->SetNdivisions(505);
  hist->SetLineWidth(2);
  hist->SetLineColor(colors[color]);
  hist->SetMarkerColor(colors[color]);
  hist->SetMarkerStyle(styles[color]);
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  hist->GetXaxis()->SetTitleFont(62); 
  hist->GetYaxis()->SetTitleFont(62);  
  hist->GetXaxis()->SetTitleSize(0.055);
  hist->GetYaxis()->SetTitleSize(0.055);
  hist->GetXaxis()->SetTitleOffset(1.0);
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(.045);
  hist->GetXaxis()->SetLabelSize(.045);
  hist->SetMarkerSize(1.8);
}

/*--------------------------------------------------------------------*/
void makeNiceCanv(TCanvas *c){
/*--------------------------------------------------------------------*/  

  //c->SetGridx();
  c->SetGridy();
  c->cd()->SetBottomMargin(0.12);
  c->cd()->SetLeftMargin(0.14);
  c->cd()->SetRightMargin(0.07);
  c->cd()->SetTopMargin(0.08);  

}

/*--------------------------------------------------------------------*/
std::pair<Double_t,Double_t> getExtrema(TObjArray *array)
/*--------------------------------------------------------------------*/
{
  TIter histoIter(array);
  TH1F *h = NULL;
  Double_t the_max(0.);
  Double_t the_min(999.);

  while((h=(TH1F*)histoIter.Next())){
    Double_t this_max = h->GetMaximum();
    Double_t this_min = h->GetMinimum();
    if(this_max>the_max) the_max = this_max;
    if(this_min<the_min) the_min = this_min;
  }

  std::pair <Double_t,Double_t> result;
  result = std::make_pair(the_min,the_max);

  //std::cout<<"Minimum: "<<the_min<<", Maximum: "<<the_max<<std::endl;

  return result;
}
