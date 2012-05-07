//////////////////////////////////////////////////////////
//
// Simple script to run ROOT analysis job
// and extract event yields
//
// Usage:
//
// root -b
// [0] .L extractSummaryTable.C+g                               
// [1] extactSummaryTable("yourSample.xml")
//
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
#include<fstream>
#include<sstream>
#include<cmath>
#include<cassert>
#include <vector>
#include <algorithm>
#include <iomanip>

using namespace std;

// --- helper class and functions --- //
#include "Sample.C" 
#include "SampleList.C"

struct Data {
  string cutName;
  double val;
  double err;
  Data():
    cutName(""),
    val(0),
    err(0)
  {;}

  void operator+=(const Data & rhs){
    this->val += rhs.val;
    this->err = sqrt(this->err*this->err + rhs.err*rhs.err);
  }
  
  void print() {
    std::cout << cutName << " " << val << " " << err << std::endl;
  }
};

void produceEvtYiedlsTable(SampleList allSamples, const Int_t nOfFiles,TString theCorrection);

//-- main method --//
void extractSummaryTable(TString sampleFileName_){
  
  SampleList allSamples;

  TString theCorrection="";
  if(sampleFileName_.Contains("All")){
    theCorrection="PU+T\\&P+b-eff";
  }else if(sampleFileName_.Contains("TNP")){
    theCorrection="PU+T\\&P";
  } else if(sampleFileName_.Contains("PU")){
    theCorrection="PU";
  } else if(sampleFileName_.Contains("OOB")){
    theCorrection="Out-of-the-box";
  } else {
    cout<<"unknown correction type"<<endl;
    return;
  }
  
  if ( allSamples.ParseFile(sampleFileName_)  == 0  ) { 
    produceEvtYiedlsTable(allSamples,allSamples.GetSize(),theCorrection);
  }
  gSystem->Sleep(50);
  TString processline = ".! pdflatex EventYields.tex";
  gROOT->ProcessLine(processline.Data());

}

// -- produce table of event yields --//
void produceEvtYiedlsTable(SampleList allSamples, const Int_t nOfFiles, TString theCorrection){
  
  // ------------------------------------- Initialization -----------------------------------------------------------------

  Sample *aSample = 0;  
  TIter nextsample(allSamples.listOfSample);  
 
  TFile   *file[nOfFiles];  
  const unsigned int  nMC= nOfFiles-1;

  Int_t jj(0);
  while (aSample = (Sample*)nextsample() ) {
    file[jj] = TFile::Open(aSample->inputrootfile_);    
    jj++;
  }

  Sample *bSample = 0;  

  ofstream outfile_;
  outfile_.open ("EventYields.tex");
  outfile_.precision(2);

  // ---------------------------------------  Start Retrieving data ----------------------------------------------------

  TH1* histoMuSSVHPT[nOfFiles];
  TH1* histoMuSSVHEM[nOfFiles]; 
  TH1* histoEleSSVHEM[nOfFiles];
  TH1* histoEleSSVHPT[nOfFiles];
  
  for(Int_t j=0; j < nOfFiles; j++) { 

    bSample = (Sample*)allSamples.listOfSample->At(j);
    
    histoMuSSVHEM[j] = (TH1F*)file[j]->Get("finaldistros_ssvhem/EventYields/eventCategoryMu");
    histoMuSSVHPT[j]  = (TH1F*)file[j]->Get("finaldistros_ssvhpt/EventYields/eventCategoryMu");
    histoEleSSVHEM[j] = (TH1F*)file[j]->Get("finaldistros_ssvhem/EventYields/eventCategoryEle");
    histoEleSSVHPT[j]  = (TH1F*)file[j]->Get("finaldistros_ssvhpt/EventYields/eventCategoryEle");

    histoMuSSVHEM[j]->Scale(bSample->weight_);
    histoMuSSVHPT[j]->Scale(bSample->weight_);                   
    histoEleSSVHEM[j]->Scale(bSample->weight_); 
    histoEleSSVHPT[j]->Scale(bSample->weight_); 

    histoMuSSVHEM[j]->Sumw2(); 
    histoMuSSVHPT[j]->Sumw2(); 
    histoEleSSVHEM[j]->Sumw2(); 
    histoEleSSVHPT[j]->Sumw2();  
  }
    
  std::vector<Data> lDataMu;
  std::vector<Data> lDataEle;
  std::vector<Data> lDataAll;

  std::vector<Data> lMCMu[nMC];
  std::vector<Data> lMCEle[nMC];
  std::vector<Data> lMCAll[nMC];

  // SSVHEM selection
  for(Int_t bin=3;  bin<=histoMuSSVHEM[0]->GetNbinsX(); bin++){
    
    // --- muons

    Data lTmpMu;
    lTmpMu.cutName = histoMuSSVHEM[0]->GetXaxis()->GetBinLabel(bin);
    lTmpMu.val  = histoMuSSVHEM[0]->GetBinContent(bin);
    lTmpMu.err = histoMuSSVHEM[0]->GetBinError(bin); 
	 
    if (lTmpMu.val==0 && lTmpMu.err==0) break;

    lDataMu.push_back(lTmpMu);
    // lTmpMu.print();

    // --- electrons

    Data lTmpEle;
    lTmpEle.cutName = histoEleSSVHEM[0]->GetXaxis()->GetBinLabel(bin);
    lTmpEle.val  = histoEleSSVHEM[0]->GetBinContent(bin);
    lTmpEle.err = histoEleSSVHEM[0]->GetBinError(bin); 
	 
    if (lTmpEle.val==0 && lTmpEle.err==0) break;

    lDataEle.push_back(lTmpEle);
    // lTmpEle.print();
 
    //--- combined    
    Data lTmpAll;
    lTmpAll += lTmpMu;
    lTmpAll += lTmpEle;
    if (lTmpAll.val==0 && lTmpAll.err==0) break;
    lDataAll.push_back(lTmpAll);
    // lTmpAll.print();

    // MC
    for (unsigned int iMC(1); iMC<=nMC; ++iMC){
      Data lTmpMCAll;

      lTmpMu.cutName = histoMuSSVHEM[iMC]->GetXaxis()->GetBinLabel(bin);
      lTmpMu.val  = histoMuSSVHEM[iMC]->GetBinContent(bin);
      lTmpMu.err = histoMuSSVHEM[iMC]->GetBinError(bin); 
      lMCMu[iMC-1].push_back(lTmpMu);
      // lTmpMu.print();
      
      lTmpEle.cutName = histoEleSSVHEM[iMC]->GetXaxis()->GetBinLabel(bin);
      lTmpEle.val  = histoEleSSVHEM[iMC]->GetBinContent(bin);
      lTmpEle.err = histoEleSSVHEM[iMC]->GetBinError(bin); 
      lMCEle[iMC-1].push_back(lTmpEle);
      //lTmpEle.print();

      lTmpMCAll += lTmpMu;
      lTmpMCAll += lTmpEle;
      lMCAll[iMC-1].push_back(lTmpMCAll);
    }
  }

  // SSVHPT selection
  for(Int_t bin=8;  bin<=histoMuSSVHPT[0]->GetNbinsX(); bin++){
    
    // --- muons
    Data lTmpMu;
    lTmpMu.cutName = histoMuSSVHPT[0]->GetXaxis()->GetBinLabel(bin);
    lTmpMu.val  = histoMuSSVHPT[0]->GetBinContent(bin);
    lTmpMu.err = histoMuSSVHPT[0]->GetBinError(bin); 
    if (lTmpMu.val==0 && lTmpMu.err==0) break;
    lDataMu.push_back(lTmpMu);
    // lTmpMu.print();

    // --- electrons
    Data lTmpEle;
    lTmpEle.cutName = histoEleSSVHPT[0]->GetXaxis()->GetBinLabel(bin);
    lTmpEle.val  = histoEleSSVHPT[0]->GetBinContent(bin);
    lTmpEle.err = histoEleSSVHPT[0]->GetBinError(bin); 
    if (lTmpEle.val==0 && lTmpEle.err==0) break;
    lDataEle.push_back(lTmpEle);
    // lTmpEle.print();
 
    //--- combined 
    Data lTmpAll;
    lTmpAll += lTmpMu;
    lTmpAll += lTmpEle;
    if (lTmpAll.val==0 && lTmpAll.err==0) break;
    lDataAll.push_back(lTmpAll);
    // lTmpAll.print();
    
    // MC
    for (unsigned int iMC(1); iMC<=nMC; ++iMC){
      Data lTmpMCAll;
      lTmpMu.cutName = histoMuSSVHPT[iMC]->GetXaxis()->GetBinLabel(bin);
      lTmpMu.val  = histoMuSSVHPT[iMC]->GetBinContent(bin);
      lTmpMu.err = histoMuSSVHPT[iMC]->GetBinError(bin); 
      lMCMu[iMC-1].push_back(lTmpMu);
      // lTmpMu.print();
      
      lTmpEle.cutName = histoEleSSVHPT[iMC]->GetXaxis()->GetBinLabel(bin);
      lTmpEle.val  = histoEleSSVHPT[iMC]->GetBinContent(bin);
      lTmpEle.err = histoEleSSVHPT[iMC]->GetBinError(bin); 
      lMCEle[iMC-1].push_back(lTmpEle);
      //lTmpEle.print();

      lTmpMCAll += lTmpMu;
      lTmpMCAll += lTmpEle;
      lMCAll[iMC-1].push_back(lTmpMCAll);
    }
  }

  // ---------------------------------------  Start Printing ----------------------------------------------------

  //---------------------------------------------------------------------------------
  // Muons
  //---------------------------------------------------------------------------------

  outfile_ << "\\documentclass{beamer}"<< std::endl
	   <<"\\setbeamertemplate{navigation symbols}{} "<<std::endl
	   <<"\\begin{document}" << std::endl
	   <<"\\begin{frame}" << std::endl
	   << "\\frametitle{Efficiency tables for $ \\ell \\ell$ +b(b) selection ("<<theCorrection<<")}"
	   << std::endl << std::endl
	   << "\\begin{table}" << std::endl
       	   << "{\\tiny" << std::endl
	   << "\\caption{$\\mu\\mu$ +b event selection} " << std::endl
   	   << "\\hspace*{-0.5cm}\\begin{tabular}{|l|c|c|" ;
  for (unsigned int i(0); i<nMC; i++){
    outfile_ << "c|";
  }
  outfile_ << "}" << std::endl
	   << "\\hline" << std::endl;
  outfile_ << "Cut / Sample & " ;
  if (nMC == 4) outfile_ << " tt & Z+b & Z+c & Z+l & " ;
  outfile_ << " SumMC & Data \\\\" 
	   << std::endl
	   << "\\hline" 
	   << std::endl;
  
  TString  lCutName[9]={"ll","ll tight","l-HLT match","Z(ll)","Z+j","Z+b (HE)","Z+b (HE)+MET","Z+b (HP)","Z+b (HP)+MET"};
  const unsigned int nCuts = 9;

  for (unsigned int iC(0); iC<nCuts; ++iC){
    outfile_ << lCutName[iC] << " & " ;
    Data lTotMCMu;
    
    for (unsigned int iMC(0); iMC<nMC; ++iMC){
      outfile_ << std::fixed << std::setprecision(0) 
    	       << " $ " << lMCMu[iMC][iC].val << " \\pm " << lMCMu[iMC][iC].err << "$ & ";
      lTotMCMu += lMCMu[iMC][iC];
    }
    
    outfile_ << std::fixed << std::setprecision(0)
    	     << " $ " << lTotMCMu.val << " \\pm " << lTotMCMu.err << "$ & "
    	     << " $ " << lDataMu[iC].val << " \\pm " << lDataMu[iC].err << "$ \\\\"
    	     << std::endl;
  }
  
  outfile_ << "\\hline" << std::endl
	   << "\\end{tabular}" << std::endl
	   << "}" << std::endl
	   << "\\end{table}" << std::endl << std::endl;
 
  //---------------------------------------------------------------------------------
  // Electrons
  //---------------------------------------------------------------------------------

  outfile_ << "\\begin{table}" << std::endl
       	   << "{\\tiny" << std::endl
	   << "\\caption{$ee$ +b event selection} " << std::endl
   	   << "\\hspace*{-0.5cm}\\begin{tabular}{|l|c|c|" ;
  for (unsigned int i(0); i<nMC; i++){
    outfile_ << "c|";
  }
  outfile_ << "}" << std::endl
	   << "\\hline" << std::endl;
  outfile_ << "Cut / Sample & " ;
  if (nMC == 4) outfile_ << " tt & Z+b & Z+c & Z+l & " ;
  outfile_ << " SumMC & Data \\\\" 
	   << std::endl
	   << "\\hline" 
	   << std::endl;
  
  for (unsigned int iC(0); iC<nCuts; ++iC){
    outfile_ << lCutName[iC] << " & " ;
    Data lTotMCEle;
    
    for (unsigned int iMC(0); iMC<nMC; ++iMC){
      outfile_ << std::fixed << std::setprecision(0) 
    	       << " $ " << lMCEle[iMC][iC].val << " \\pm " << lMCEle[iMC][iC].err << "$ & ";
      lTotMCEle += lMCEle[iMC][iC];
    }
    
    outfile_ << std::fixed << std::setprecision(0)
    	     << " $ " << lTotMCEle.val << " \\pm " << lTotMCEle.err << "$ & "
    	     << " $ " << lDataEle[iC].val << " \\pm " << lDataEle[iC].err << "$ \\\\"
    	     << std::endl;
  }
  
  outfile_ << "\\hline" << std::endl
	   << "\\end{tabular}" << std::endl
	   << "}" << std::endl
	   << "\\end{table}" << std::endl << std::endl
	   << "\\end{frame}"  << std::endl;
  // << "\\end{document}"  << std::endl;
  
  //---------------------------------------------------------------------------------
  // Combined
  //---------------------------------------------------------------------------------

  outfile_ <<"\\begin{frame}" << std::endl
	   << "\\frametitle{Efficiency tables for $ \\ell \\ell$+b(b) selection  ("<<theCorrection<<")}"
	   << std::endl << std::endl
	   << "\\begin{table}" << std::endl
       	   << "{\\tiny" << std::endl
  	   << "\\caption{Combined $\\ell \\ell$ +b event selection} " << std::endl
   	   << "\\hspace*{-0.5cm}\\begin{tabular}{|l|c|c|" ;
  for (unsigned int i(0); i<nMC; i++){
    outfile_ << "c|";
  }
  outfile_ << "}" << std::endl
  	   << "\\hline" << std::endl;
  outfile_ << "Cut / Sample & " ;
  if (nMC == 4) outfile_ << " tt & Z+b & Z+c & Z+l & " ;
  outfile_ << " SumMC & Data \\\\" 
  	   << std::endl
  	   << "\\hline" 
  	   << std::endl;
  
  for (unsigned int iC(0); iC<nCuts; ++iC){
    outfile_ << lCutName[iC] << " & " ;
    Data lTotMCAll;
    
    for (unsigned int iMC(0); iMC<nMC; ++iMC){
      
      outfile_ << std::fixed << std::setprecision(0) 
       	       << " $ " << lMCAll[iMC][iC].val << " \\pm " << lMCAll[iMC][iC].err << "$ & ";
      lTotMCAll += lMCAll[iMC][iC];
    }
    
    outfile_ << std::fixed << std::setprecision(0)
	     << " $ " << lTotMCAll.val << " \\pm " << lTotMCAll.err << "$ & "
     	     << " $ " << lDataAll[iC].val << " \\pm " << lDataAll[iC].err << "$ \\\\"
    	     << std::endl;
  }
  
  outfile_ << "\\hline" << std::endl
  	   << "\\end{tabular}" << std::endl
  	   << "}" << std::endl
  	   << "\\end{table}" << std::endl << std::endl
  	   << "\\end{frame}"  << std::endl
  	   << "\\end{document}"  << std::endl;

  outfile_.close();

}
