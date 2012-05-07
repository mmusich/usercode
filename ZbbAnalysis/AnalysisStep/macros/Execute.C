//////////////////////////////////////////////////////////
//
// Simple script to run ROOT analysis job:
//
// Usage:
//
// root -b
// [0] .L Execute.C+g                                       
// [1] Execute("yourSample.xml",true)
//
// (true) for DATA/MC ratio  ----  (false) for no ratio
// Original Author: M. Musich INFN Torino
//
//////////////////////////////////////////////////////////

#include "TSystem.h"
#include <Riostream.h>
#include <TROOT.h>
#include <TString.h>
#include <TList.h>
#include <TFile.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TKey.h>
#include <fstream>
#include <iostream>
#include <TStopwatch.h>

using namespace std;

// --- helper class and functions --- //
#include "Sample.C" 
#include "SampleList.C"
#include "HistoDrawHelper.C"
#include "diow.C"

// -- produce table of event yields --//
void produceEvtYiedlsTable(SampleList allSamples, const Int_t nOfFiles){
  
  Sample *aSample = 0;  
  TIter nextsample(allSamples.listOfSample);  
  
  TFile   *file[nOfFiles];  

  Int_t jj(0);
  while (aSample = (Sample*)nextsample() ) {
    file[jj] = TFile::Open(aSample->inputrootfile_);    
    jj++;
  }

  Sample *bSample = 0;  

  ofstream outfile1_;
  outfile1_.open ("EventYields_ssvhem.txt");
  outfile1_.precision(2);
  
  ofstream outfile2_;
  outfile2_.open ("EventYields_ssvhpt.txt");
  outfile2_.precision(2);

  ofstream outfile3_;
  outfile3_.open ("EventYields_csvm.txt");
  outfile3_.precision(2);

  ofstream outfile4_;
  outfile4_.open ("EventYields_csvt.txt");
  outfile4_.precision(2);

  TH1* histoMuSSVHPT[nOfFiles];
  TH1* histoMuSSVHEM[nOfFiles];
  TH1* histoMuCSVT[nOfFiles];
  TH1* histoMuCSVM[nOfFiles];
  THStack *hsMuHE   = new THStack("eventCategoryMuSSVHEM","eventCategoryMuSSVHEM");
  THStack *hsMuHP   = new THStack("eventCategoryMuSSVHPT","eventCategoryMuSSVHPT");
  THStack *hsMuCSVM = new THStack("eventCategoryMuCSVM","eventCategoryMuCSVM");
  THStack *hsMuCSVT = new THStack("eventCategoryMuCSVT","eventCategoryMuCSVT");
  
  TH1* histoEleSSVHEM[nOfFiles];
  TH1* histoEleSSVHPT[nOfFiles];
  TH1* histoEleCSVM[nOfFiles];
  TH1* histoEleCSVT[nOfFiles];
  THStack *hsEleHE = new THStack("eventCategoryEleSSVHEM","eventCategoryEleSSVHEM");
  THStack *hsEleHP = new THStack("eventCategoryEleSSVHPT","eventCategoryEleSSVHPT");
  THStack *hsEleCSVM = new THStack("eventCategoryEleCSVM","eventCategoryEleCSVM");
  THStack *hsEleCSVT = new THStack("eventCategoryEleCSVT","eventCategoryEleCSVT");

  TH1* histoAllSSVHEM[nOfFiles];
  TH1* histoAllSSVHPT[nOfFiles];
  TH1* histoAllCSVM[nOfFiles];
  TH1* histoAllCSVT[nOfFiles];

  outfile1_<<"-------------------------------------------------------------------------- Muon Channel --------------------------------------------------------------------------"<<endl;
  outfile2_<<"-------------------------------------------------------------------------- Muon Channel --------------------------------------------------------------------------"<<endl;
  outfile3_<<"-------------------------------------------------------------------------- Muon Channel --------------------------------------------------------------------------"<<endl;
  outfile4_<<"-------------------------------------------------------------------------- Muon Channel --------------------------------------------------------------------------"<<endl;

  for(Int_t j=0; j < nOfFiles; j++) { 

    bSample = (Sample*)allSamples.listOfSample->At(j);
    
    histoMuSSVHEM[j]  = (TH1F*)file[j]->Get("finaldistros_ssvhem/EventYields/eventCategoryMu");
    histoMuSSVHPT[j]  = (TH1F*)file[j]->Get("finaldistros_ssvhpt/EventYields/eventCategoryMu");
    histoEleSSVHEM[j] = (TH1F*)file[j]->Get("finaldistros_ssvhem/EventYields/eventCategoryEle");
    histoEleSSVHPT[j] = (TH1F*)file[j]->Get("finaldistros_ssvhpt/EventYields/eventCategoryEle");

    histoMuCSVM[j]  = (TH1F*)file[j]->Get("finaldistros_csvm/EventYields/eventCategoryMu");
    histoMuCSVT[j]  = (TH1F*)file[j]->Get("finaldistros_csvt/EventYields/eventCategoryMu");
    histoEleCSVM[j] = (TH1F*)file[j]->Get("finaldistros_csvm/EventYields/eventCategoryEle");
    histoEleCSVT[j] = (TH1F*)file[j]->Get("finaldistros_csvt/EventYields/eventCategoryEle");

    histoMuSSVHEM[j]->Scale(bSample->weight_);
    histoMuSSVHPT[j]->Scale(bSample->weight_);                   
    histoEleSSVHEM[j]->Scale(bSample->weight_); 
    histoEleSSVHPT[j]->Scale(bSample->weight_); 
    
    histoMuCSVM[j]->Scale(bSample->weight_);  
    histoMuCSVT[j]->Scale(bSample->weight_);  
    histoEleCSVM[j]->Scale(bSample->weight_);
    histoEleCSVT[j]->Scale(bSample->weight_);

    //histoMuSSVHEM[j]->Sumw2(); 
    //histoMuSSVHPT[j]->Sumw2(); 
    //histoEleSSVHEM[j]->Sumw2(); 
    //histoEleSSVHPT[j]->Sumw2();  
   
    histoAllSSVHEM[j]=(TH1*)histoMuSSVHEM[j]->Clone((TString)histoMuSSVHEM[j]->GetName()+"_new");
    histoAllSSVHPT[j]=(TH1*)histoMuSSVHPT[j]->Clone((TString)histoMuSSVHPT[j]->GetName()+"_new");
    histoAllCSVM[j]=(TH1*)histoMuCSVM[j]->Clone((TString)histoMuCSVM[j]->GetName()+"_new");
    histoAllCSVT[j]=(TH1*)histoMuCSVT[j]->Clone((TString)histoMuCSVT[j]->GetName()+"_new");

    histoAllSSVHEM[j]->Add(histoEleSSVHEM[j]);
    histoAllSSVHPT[j]->Add(histoEleSSVHPT[j]);
    histoAllCSVM[j]->Add(histoEleCSVM[j]);
    histoAllCSVT[j]->Add(histoEleCSVT[j]);
 
    if(j!=0){
      hsMuHE->Add(histoMuSSVHEM[j]);
      hsMuHP->Add(histoMuSSVHPT[j]);
      hsEleHE->Add(histoEleSSVHEM[j]);
      hsEleHP->Add(histoEleSSVHPT[j]);
      
      hsMuCSVM->Add(histoMuCSVM[j]);
      hsMuCSVT->Add(histoMuCSVT[j]);
      hsEleCSVM->Add(histoEleCSVM[j]);
      hsEleCSVT->Add(histoEleCSVT[j]);
    }
    
    if(j==0){
      outfile1_<<setw(15)<<" Selection"<<"|";
      outfile2_<<setw(15)<<" Selection"<<"|";
      outfile3_<<setw(15)<<" Selection"<<"|";
      outfile4_<<setw(15)<<" Selection"<<"|";
    } 
    outfile1_<<setw(21)<<bSample->histolabel_<<"|";
    outfile2_<<setw(21)<<bSample->histolabel_<<"|";
    outfile3_<<setw(21)<<bSample->histolabel_<<"|";
    outfile4_<<setw(21)<<bSample->histolabel_<<"|";
  }
  
  outfile1_<<setw(21)<<"Sum MC"<<"|"<<endl;
  outfile2_<<setw(21)<<"Sum MC"<<"|"<<endl;
  outfile3_<<setw(21)<<"Sum MC"<<"|"<<endl;
  outfile4_<<setw(21)<<"Sum MC"<<"|"<<endl;
  
  TH1F *sumMuHE  = (TH1F*)hsMuHE->GetStack()->Last();
  TH1F *sumMuHP  = (TH1F*)hsMuHP->GetStack()->Last();
  TH1F *sumEleHE = (TH1F*)hsEleHE->GetStack()->Last();
  TH1F *sumEleHP = (TH1F*)hsEleHP->GetStack()->Last();
  
  TH1F *sumMuCSVM  = (TH1F*)hsMuCSVM->GetStack()->Last();
  TH1F *sumMuCSVT  = (TH1F*)hsMuCSVT->GetStack()->Last();
  TH1F *sumEleCSVM = (TH1F*)hsEleCSVM->GetStack()->Last();
  TH1F *sumEleCSVT = (TH1F*)hsEleCSVT->GetStack()->Last();

  TH1F *sumAllHE =  (TH1F*)sumMuHE->Clone((TString)sumMuHE->GetName()+"_new");
  TH1F *sumAllHP =  (TH1F*)sumMuHP->Clone((TString)sumMuHP->GetName()+"_new");
  
  sumAllHE->Add(sumEleHE);
  sumAllHP->Add(sumEleHP);

  TH1F *sumAllCSVM =  (TH1F*)sumMuCSVM->Clone((TString)sumMuCSVM->GetName()+"_new");
  TH1F *sumAllCSVT =  (TH1F*)sumMuCSVT->Clone((TString)sumMuCSVT->GetName()+"_new");
  
  sumAllCSVM->Add(sumEleCSVM);
  sumAllCSVT->Add(sumEleCSVT);

  for(Int_t bin=3;  bin<=histoMuSSVHEM[0]->GetNbinsX(); bin++){

    TString binlabel = histoMuSSVHEM[0]->GetXaxis()->GetBinLabel(bin);
    outfile1_<<setw(15)<<binlabel<<"|";
    outfile2_<<setw(15)<<binlabel<<"|";
    outfile3_<<setw(15)<<binlabel<<"|";
    outfile4_<<setw(15)<<binlabel<<"|";
    for(Int_t j=0; j < nOfFiles; j++) { 
      
      float bincontentMuHE  = histoMuSSVHEM[j]->GetBinContent(bin);
      float binerrorMuHE    = histoMuSSVHEM[j]->GetBinError(bin); 	
      float bincontentMuHP  = histoMuSSVHPT[j]->GetBinContent(bin);
      float binerrorMuHP    = histoMuSSVHPT[j]->GetBinError(bin);
      
      float bincontentMuCSVM  = histoMuCSVM[j]->GetBinContent(bin);
      float binerrorMuCSVM    = histoMuCSVM[j]->GetBinError(bin); 	
      float bincontentMuCSVT  = histoMuCSVT[j]->GetBinContent(bin);
      float binerrorMuCSVT    = histoMuCSVT[j]->GetBinError(bin);
	
      outfile1_<<fixed<<setw(12)<<bincontentMuHE<<"+/-"<<setw(6)<<binerrorMuHE<<"|";
      outfile2_<<fixed<<setw(12)<<bincontentMuHP<<"+/-"<<setw(6)<<binerrorMuHP<<"|";
      outfile3_<<fixed<<setw(12)<<bincontentMuCSVM<<"+/-"<<setw(6)<<binerrorMuCSVM<<"|";
      outfile4_<<fixed<<setw(12)<<bincontentMuCSVT<<"+/-"<<setw(6)<<binerrorMuCSVT<<"|";
    }
   
    outfile1_<<fixed<<setw(12)<<sumMuHE->GetBinContent(bin)<<"+/-"<<setw(6)<<sumMuHE->GetBinError(bin)<<"|";
    outfile2_<<fixed<<setw(12)<<sumMuHP->GetBinContent(bin)<<"+/-"<<setw(6)<<sumMuHP->GetBinError(bin)<<"|";
    outfile3_<<fixed<<setw(12)<<sumMuCSVM->GetBinContent(bin)<<"+/-"<<setw(6)<<sumMuCSVM->GetBinError(bin)<<"|";
    outfile4_<<fixed<<setw(12)<<sumMuCSVT->GetBinContent(bin)<<"+/-"<<setw(6)<<sumMuCSVT->GetBinError(bin)<<"|";
    outfile1_<<std::endl;
    outfile2_<<std::endl;
    outfile3_<<std::endl;
    outfile4_<<std::endl;
  }

  //outfile1_<<setw(80)<<setfill('-')<<endl;
  //outfile2_<<setw(80)<<setfill('-')<<endl;

  // now for the Electron Channel
  outfile1_<<"---------------------------------------------------------------------- Electron Channel --------------------------------------------------------------------------"<<endl;
  outfile2_<<"---------------------------------------------------------------------- Electron Channel --------------------------------------------------------------------------"<<endl;
  outfile3_<<"---------------------------------------------------------------------- Electron Channel --------------------------------------------------------------------------"<<endl;
  outfile4_<<"---------------------------------------------------------------------- Electron Channel --------------------------------------------------------------------------"<<endl;

  for(Int_t j=0; j < nOfFiles; j++) { 
    bSample = (Sample*)allSamples.listOfSample->At(j);
    if(j==0){
      outfile1_<<setw(15)<<" Selection"<<"|";
      outfile2_<<setw(15)<<" Selection"<<"|";
      outfile3_<<setw(15)<<" Selection"<<"|";
      outfile4_<<setw(15)<<" Selection"<<"|";
    } 
    outfile1_<<setw(21)<<bSample->histolabel_<<"|";
    outfile2_<<setw(21)<<bSample->histolabel_<<"|";
    outfile3_<<setw(21)<<bSample->histolabel_<<"|";
    outfile4_<<setw(21)<<bSample->histolabel_<<"|";
  }
  
  outfile1_<<setw(21)<<"Sum MC"<<"|"<<endl;
  outfile2_<<setw(21)<<"Sum MC"<<"|"<<endl;
  outfile3_<<setw(21)<<"Sum MC"<<"|"<<endl;
  outfile4_<<setw(21)<<"Sum MC"<<"|"<<endl;

  for(Int_t bin=3;  bin<=histoEleSSVHEM[0]->GetNbinsX(); bin++){
    
    TString binlabel = histoEleSSVHEM[0]->GetXaxis()->GetBinLabel(bin);
    outfile1_<<setw(15)<<binlabel<<"|";
    outfile2_<<setw(15)<<binlabel<<"|";
    outfile3_<<setw(15)<<binlabel<<"|";
    outfile4_<<setw(15)<<binlabel<<"|";
    for(Int_t j=0; j < nOfFiles; j++) { 
      float bincontentEleHE   = histoEleSSVHEM[j]->GetBinContent(bin);
      float binerrorEleHE     = histoEleSSVHEM[j]->GetBinError(bin); 	
      float bincontentEleHP   = histoEleSSVHPT[j]->GetBinContent(bin);
      float binerrorEleHP     = histoEleSSVHPT[j]->GetBinError(bin); 	
      float bincontentEleCSVM = histoEleCSVM[j]->GetBinContent(bin);
      float binerrorEleCSVM   = histoEleCSVM[j]->GetBinError(bin); 	
      float bincontentEleCSVT = histoEleCSVT[j]->GetBinContent(bin);
      float binerrorEleCSVT   = histoEleCSVT[j]->GetBinError(bin); 
      outfile1_<<fixed<<setw(12)<<bincontentEleHE<<"+/-"<<setw(6)<<binerrorEleHE<<"|";
      outfile2_<<fixed<<setw(12)<<bincontentEleHP<<"+/-"<<setw(6)<<binerrorEleHP<<"|";
      outfile3_<<fixed<<setw(12)<<bincontentEleCSVM<<"+/-"<<setw(6)<<binerrorEleCSVM<<"|";
      outfile4_<<fixed<<setw(12)<<bincontentEleCSVT<<"+/-"<<setw(6)<<binerrorEleCSVT<<"|";
    }
   
    outfile1_<<fixed<<setw(12)<<sumEleHE->GetBinContent(bin)<<"+/-"<<setw(6)<<sumEleHE->GetBinError(bin)<<"|";
    outfile2_<<fixed<<setw(12)<<sumEleHP->GetBinContent(bin)<<"+/-"<<setw(6)<<sumEleHP->GetBinError(bin)<<"|";
    outfile3_<<fixed<<setw(12)<<sumEleCSVM->GetBinContent(bin)<<"+/-"<<setw(6)<<sumEleCSVM->GetBinError(bin)<<"|";
    outfile4_<<fixed<<setw(12)<<sumEleCSVT->GetBinContent(bin)<<"+/-"<<setw(6)<<sumEleCSVT->GetBinError(bin)<<"|";
    outfile1_<<std::endl;
    outfile2_<<std::endl;
    outfile3_<<std::endl;
    outfile4_<<std::endl;
  }

  // now for the combining channels
  outfile1_<<"---------------------------------------------------------------------- Combined Channel --------------------------------------------------------------------------"<<endl;
  outfile2_<<"---------------------------------------------------------------------- Combined Channel --------------------------------------------------------------------------"<<endl;
  outfile3_<<"---------------------------------------------------------------------- Combined Channel --------------------------------------------------------------------------"<<endl;
  outfile4_<<"---------------------------------------------------------------------- Combined Channel --------------------------------------------------------------------------"<<endl;

  TString CombinedBinLabels[10]={"ll","ll tight","l-HLT match","Z(ll)","Z+j","Z+jvtx","Z+b","Z+1b","Z+b+l","Z+b+2l"};

  for(Int_t j=0; j < nOfFiles; j++) { 
    bSample = (Sample*)allSamples.listOfSample->At(j);
    if(j==0){
      outfile1_<<setw(15)<<" Selection"<<"|";
      outfile2_<<setw(15)<<" Selection"<<"|";
      outfile3_<<setw(15)<<" Selection"<<"|";
      outfile4_<<setw(15)<<" Selection"<<"|";
    } 
    outfile1_<<setw(21)<<bSample->histolabel_<<"|";
    outfile2_<<setw(21)<<bSample->histolabel_<<"|";
    outfile3_<<setw(21)<<bSample->histolabel_<<"|";
    outfile4_<<setw(21)<<bSample->histolabel_<<"|";
  }
  
  outfile1_<<setw(21)<<"Sum MC"<<"|"<<endl;
  outfile2_<<setw(21)<<"Sum MC"<<"|"<<endl;
  outfile3_<<setw(21)<<"Sum MC"<<"|"<<endl;
  outfile4_<<setw(21)<<"Sum MC"<<"|"<<endl;

  for(Int_t bin=3;  bin<=histoAllSSVHEM[0]->GetNbinsX(); bin++){
    
    TString binlabel = histoAllSSVHEM[0]->GetXaxis()->GetBinLabel(bin);
    outfile1_<<setw(15)<<CombinedBinLabels[bin-3]<<"|";
    outfile2_<<setw(15)<<CombinedBinLabels[bin-3]<<"|";
    outfile3_<<setw(15)<<CombinedBinLabels[bin-3]<<"|";
    outfile4_<<setw(15)<<CombinedBinLabels[bin-3]<<"|";

    for(Int_t j=0; j < nOfFiles; j++) { 
      float bincontentAllHE  = histoAllSSVHEM[j]->GetBinContent(bin);
      float binerrorAllHE    = histoAllSSVHEM[j]->GetBinError(bin); 	
      float bincontentAllHP  = histoAllSSVHPT[j]->GetBinContent(bin);
      float binerrorAllHP    = histoAllSSVHPT[j]->GetBinError(bin); 	
      float bincontentAllCSVM  = histoAllCSVM[j]->GetBinContent(bin);
      float binerrorAllCSVM    = histoAllCSVM[j]->GetBinError(bin); 	
      float bincontentAllCSVT  = histoAllCSVT[j]->GetBinContent(bin);
      float binerrorAllCSVT    = histoAllCSVT[j]->GetBinError(bin); 
      outfile1_<<fixed<<setw(12)<<bincontentAllHE<<"+/-"<<setw(6)<<binerrorAllHE<<"|";
      outfile2_<<fixed<<setw(12)<<bincontentAllHP<<"+/-"<<setw(6)<<binerrorAllHP<<"|";
      outfile3_<<fixed<<setw(12)<<bincontentAllCSVM<<"+/-"<<setw(6)<<binerrorAllCSVM<<"|";
      outfile4_<<fixed<<setw(12)<<bincontentAllCSVT<<"+/-"<<setw(6)<<binerrorAllCSVT<<"|";
    }
   
    outfile1_<<fixed<<setw(12)<<sumAllHE->GetBinContent(bin)<<"+/-"<<setw(6)<<sumAllHE->GetBinError(bin)<<"|";
    outfile2_<<fixed<<setw(12)<<sumAllHP->GetBinContent(bin)<<"+/-"<<setw(6)<<sumAllHP->GetBinError(bin)<<"|";
    outfile3_<<fixed<<setw(12)<<sumAllCSVM->GetBinContent(bin)<<"+/-"<<setw(6)<<sumAllCSVM->GetBinError(bin)<<"|";
    outfile4_<<fixed<<setw(12)<<sumAllCSVT->GetBinContent(bin)<<"+/-"<<setw(6)<<sumAllCSVT->GetBinError(bin)<<"|";
    
    outfile1_<<std::endl;
    outfile1_<<std::endl;
    outfile3_<<std::endl;
    outfile4_<<std::endl;
  }

  outfile1_.close();
  outfile2_.close();
  outfile3_.close();
  outfile4_.close();

}
 
//----- main plotting function --------//
void FastStackAndSuperImposeHistosAllDir(SampleList allSamples, const Int_t nOfFiles, Bool_t StackIt, TString initialfoldername, Bool_t doMCDATAratio,Bool_t doStatTest,Int_t nMaxPlots,TFile* ff){

  Int_t thePlotCounter_(0);

  Sample *aSample = 0;  
  TIter nextsample(allSamples.listOfSample);  
   
  TFile   *file[nOfFiles];  
  TString LegLabels[nOfFiles];  

  Float_t lumiProcessed(0.);
  Int_t jj(0);
  while (aSample = ((Sample*)nextsample()) ) {
    lumiProcessed=(aSample->lumi_)/1000.;
    file[jj] = TFile::Open(aSample->inputrootfile_);    
    TString strlegend = aSample->histolabel_; 
    LegLabels[jj] = strlegend;
    cout<<"LegLabels["<<jj<<"] "<<LegLabels[jj]<<endl;
    jj++;
  }

  cout<<"Lumi processed: "<<lumiProcessed<<"/fb"<<endl;
  
  file[0]->cd(initialfoldername);
  ff->mkdir(initialfoldername);
  TDirectory *stardir = gDirectory;

  TObject *thesourcedir;
  TIter nextdir(stardir->GetListOfKeys());

  while((thesourcedir=nextdir())){

    TString dirName = thesourcedir->GetName();
       
    stardir->cd(dirName);
    TDirectory *current_sourcedir = gDirectory;
    TH1::AddDirectory(kFALSE);

    std::cout << "Reading Directory: " << dirName <<std::endl;

    ff->cd(initialfoldername);
    TDirectory *savdir = gDirectory;
    TDirectory *adir = savdir->mkdir(dirName);

    TIter next(current_sourcedir->GetListOfKeys());
    TKey *key;
  
    while ( (key = (TKey*)next())) {
      TObject *obj = key->ReadObj();
      TString objName =obj->GetName();
      thePlotCounter_++;
      if(thePlotCounter_>nMaxPlots) return;
      if (obj->IsA()->InheritsFrom( "TH1" ) && !objName.Contains("2D") && !objName.Contains("vs") && !objName.Contains("Vs") 
	  && !objName.Contains("GENP") && !objName.Contains("LHE") 
	  && !objName.Contains("UnwEventCategory")
	  ){
	
	TH1* histo[nOfFiles];
	Double_t  ytmp = doMCDATAratio ? 800 : 700 ;
	TCanvas *theCanvas=new TCanvas(dirName+objName,objName,600,ytmp);
	theCanvas->SetName(objName);
	MakeNiceCanvas(theCanvas,StackIt);
	TString theName = objName+".png";
	
	TObjArray *arrayHistos = new TObjArray();
	arrayHistos->Expand(nOfFiles);
	THStack *hs = new THStack(objName,objName);
	
	Sample *bSample = 0;  

	//std::cout << "Reading object: " <<initialfoldername << "/" << dirName << "/"<<objName <<std::endl;

	for(Int_t j=0; j < nOfFiles; j++) { 
	  bSample = (Sample*)allSamples.listOfSample->At(j);
	  histo[j] = (TH1F*)file[j]->Get(initialfoldername+"/"+dirName+"/"+objName);
	  if(StackIt){
	    histo[j]->Scale(bSample->weight_);
	    MakeNiceHistoStyle(histo[j],bSample->histocolor_,StackIt,doMCDATAratio);
	    if(initialfoldername.Contains("Pat")){ 
	      histo[j]->Sumw2();
	    }
	  } else {
	    histo[j]->Scale(histo[0]->GetSumOfWeights()/histo[j]->GetSumOfWeights());
	    MakeNiceHistoStyle(histo[j],bSample->histocolor_,StackIt,doMCDATAratio); 
	  }
	  arrayHistos->Add(histo[j]);

	  if(j!=0){
	    hs->Add(histo[j]);
	  }
	}
	
	TLatex *latexLabel = new TLatex();
	latexLabel->SetTextSize(0.04);
	latexLabel->SetTextFont(42);
	latexLabel->SetLineWidth(2);
	latexLabel->SetNDC();

	theCanvas->cd();
	
	TPad top_pad("top_pad", "top_pad",0.0,0.3, 1.0, 1.0);
	top_pad.SetFillColor(0);

	TPad bottom_pad("bottom_pad", "bottom_pad",0., 0.05, 1.0, 0.3);
	bottom_pad.SetFillColor(0);
	  
	// draw top pat only if ratio is desired
	if ( doMCDATAratio ) {	  
	  top_pad.SetBottomMargin(0);
	  top_pad.Draw();
	  top_pad.cd();
	}

	TString theHistoTitle =  histo[0]->GetXaxis()->GetTitle();
	histo[0]->SetMarkerStyle(20);
	histo[0]->SetMarkerSize(1.2);
	if( theHistoTitle.Contains("GeV") ){
	  setAxisTitle(histo[0],"Gev");
	} else if( theHistoTitle.Contains("rad")){
	  setAxisTitle(histo[0],"rad");
	} else {
	  setAxisTitle(histo[0],"");
	}
	  
	histo[0]->Draw("Pe");
	if(StackIt){
	  hs->Draw("HISTSAME");
	} else {
	  hs->Draw("nostack,hits,same");
	}
	histo[0]->Draw("esame");

	Double_t  yoffset = doMCDATAratio ? 0.00: 0.02 ;
	latexLabel->DrawLatex(0.18, 0.88+yoffset, "CMS Z+b WG");
	latexLabel->DrawLatex(0.18, 0.83+yoffset,Form("#sqrt{s} = 7 TeV  L =%.2f fb^{-1}",lumiProcessed));
	
	if(doStatTest) testKolmogorovAndChi2(histo[0],(TH1*)hs->GetStack()->Last());

	Double_t theMaximum = getMaximum(arrayHistos);
	
	if ( !doMCDATAratio ) theCanvas->Draw();
	
	TLegend *legend = MakeTLegend(arrayHistos,LegLabels);
        legend->SetFillColor(0);
	legend->Draw("same");

	if ( !doMCDATAratio ){
	  if(!objName.Contains("JetBTag")){ 
	    theCanvas->cd()->SetLogy();
	    histo[0]->GetYaxis()->SetRangeUser(3,theMaximum*45);
	  } else {
	    histo[0]->GetYaxis()->SetRangeUser(0.,theMaximum*1.6);
	  }
	} else {
	  if(!objName.Contains("JetBTag") ){ 
	    top_pad.cd()->SetLogy();
	    histo[0]->GetYaxis()->SetRangeUser(3,theMaximum*45);
	  } else {
	    histo[0]->GetYaxis()->SetRangeUser(0.,theMaximum*1.6);
	  }
	}

	if ( doMCDATAratio ) {
	// ------------------------------------------------------------------------------------------

	  theCanvas->cd();  
	  bottom_pad.SetTopMargin(0);
	  bottom_pad.Draw();

	  TH1F *sum = (TH1F*)hs->GetStack()->Last();

	  //store sum of squares of weights (if not already done)
	  //if(initialfoldername.Contains("analyzePat")) histo[0]->Sumw2();
	  //if(initialfoldername.Contains("analyzePat")) sum->Sumw2();
	  
	  TH1F *rat  = (TH1F*)histo[0]->Clone(); 
	  TH1F *rat2 = (TH1F*)histo[0]->Clone();
	  TH1F *rat3 = (TH1F*)histo[0]->Clone();
	  TH1F *rat4 = (TH1F*)histo[0]->Clone();
	  rat->SetName("Ratio"); // Clone one of the histograms
	  rat->Divide(histo[0],sum,1.,1.,"B");
	  
	  for(Int_t nbin=1; nbin<=rat2->GetNbinsX(); nbin++){
	    rat2->SetBinContent(nbin,1);
	    rat2->SetBinError(nbin,0);
	    rat3->SetBinContent(nbin,1.5);
	    rat3->SetBinError(nbin,0);
	    rat4->SetBinContent(nbin,0.5);
	    rat4->SetBinError(nbin,0);
	  }	

	  rat->SetMarkerStyle(27);
	  bottom_pad.cd();
	  rat->Draw("e");
	  rat->GetYaxis()->SetRangeUser(0.1,1.9);
	  rat->GetYaxis()->SetNdivisions(303);
	  rat->GetYaxis()->SetTitleOffset(2);
	  rat->GetXaxis()->SetTitleOffset(5);
	  
	  rat2->SetFillColor(0);
	  rat2->SetLineColor(4);
	  rat2->SetLineWidth(2);
	  rat2->Draw("same");
	  rat3->SetFillColor(0);
	  rat3->SetLineColor(4);
	  rat3->SetLineWidth(2);
	  rat3->SetLineStyle(4);
	  rat3->Draw("same");
	  rat4->SetFillColor(0);
	  rat4->SetLineColor(4);
	  rat4->SetLineWidth(2);
	  rat4->Draw("same");	
	  rat4->SetLineStyle(4);
	  
	  rat->GetYaxis()->SetTitle("DATA/MC ratio");
	  rat2->GetYaxis()->SetRangeUser(0.1,1.9);
	} else {
	  TCanvas* canvasWithRatio = new TCanvas();
	  canvasWithRatio =  DrawCanvasWithRatio(theCanvas,false);
	  canvasWithRatio->SaveAs(objName+".png");
	  ff->cd(initialfoldername+"/"+dirName);
	  canvasWithRatio->Write();
	  canvasWithRatio->Close();
	  delete canvasWithRatio;
	}

	//theCanvas->SaveAs(theName);
	//ff->cd(initialfoldername+"/"+dirName);
	//theCanvas->Write();

      }
    } // loop on histograms 
  } //loop on TDirectories

  std::cout<<"total analyzed histograms: "<<thePlotCounter_<<std::endl;

}  

//-------------------------------------------------------- EXECUTE FUNCTION ------------------------------------------------------------------//
void Execute(TString sampleFileName_,Bool_t MCDRatio=true,Bool_t doStatTest=false,Int_t nMaxPlots=10000){

  TStopwatch timer; 	 
  timer.Start();

  gROOT->Reset();

  // style
  setStyle();

  SampleList allSamples;
  if ( allSamples.ParseFile(sampleFileName_)  == 0  ) {
    TString outdirnames[6] ={"analyzePat","analyzePatAfterCleaning","finaldistros_ssvhem","finaldistros_ssvhpt","finaldistros_csvm","finaldistros_csvt"};

    TFile* fileout = new TFile("histos.root","RECREATE");
    produceEvtYiedlsTable(allSamples,allSamples.GetSize());

    for(UInt_t i=0; i<6; i++){
      TString clear = Form(".!rm \t -fr %s",outdirnames[i].Data());
      //cout<<clear<<endl;
      gROOT->ProcessLine(clear.Data());
      gSystem->mkdir(outdirnames[i]);      
      FastStackAndSuperImposeHistosAllDir(allSamples,allSamples.GetSize(),true,outdirnames[i],MCDRatio,doStatTest,nMaxPlots,fileout);
      TString processline = Form(".! mv *.png %s",outdirnames[i].Data()) ;
      gROOT->ProcessLine(processline.Data());
      gSystem->Sleep(100);
    }

    fileout->Close();

    //gROOT->LoadMacro("./diow.C++");
    for(UInt_t i=0; i<6; i++){
      diow("./"+outdirnames[i],"index.html");
    }

    timer.Stop(); 	 
    timer.Print();
    
  }
}
