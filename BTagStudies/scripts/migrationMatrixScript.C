#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <vector>
#include <map>
#include <TStyle.h>
#include <TColor.h>
#include <TROOT.h>
#include <iostream>

void resetHistoName(TH2F *h,TString dirname){
  h->SetName(h->GetName()+dirname);
}

void SetHistoStyle(TH2F *h2){
  
  h2->SetMarkerSize(4);  
  h2->SetMarkerColor(kRed);
  //h2->GetYaxis()->LabelsOption("d");
  //h2->GetXaxis()->LabelsOption("u");
  h2->GetYaxis()->CenterTitle();
  h2->GetYaxis()->SetTitleOffset(1.9);
  h2->GetXaxis()->SetTitleOffset(1.8);
  h2->GetYaxis()->SetTitleSize(0.06);
  h2->GetXaxis()->SetTitleSize(0.06);
  h2->GetXaxis()->SetTitleFont(42); 
  h2->GetYaxis()->SetTitleFont(42);  
  h2->GetXaxis()->SetLabelSize(0.05); 
  h2->GetYaxis()->SetLabelSize(0.05);
  h2->GetXaxis()->SetLabelOffset(0.02);
  h2->GetXaxis()->SetLabelFont(42);
  h2->GetYaxis()->SetLabelFont(42);
  h2->LabelsDeflate("X");
  h2->LabelsDeflate("Y");
  h2->LabelsOption("d");
}

void migrationMatrixScript(TString filename="JetByJetComparisonPlots_DeltaZ500_vs_Ideal.root"){
  
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  Double_t stops[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
  Double_t red[5]   = {1.00, 0.84, 0.61, 0.34, 0.00};
  Double_t green[5] = {1.00, 0.84, 0.61, 0.34, 0.00};
  Double_t blue[5]  = {1.00, 1.00, 1.00, 1.00, 1.00};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  gROOT->SetStyle("Pub");
  gStyle->SetPaintTextFormat("6.3g");
  //gStyle->SetPalette(1);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFontSize(0.06);
  //gStyle->SetOptStat("em");
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
  gStyle->SetStatBorderSize(1);
  gStyle->SetPadTickX(1);      // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
  gStyle->SetPadBorderMode(0); 
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadBottomMargin(0.25);
  gStyle->SetPadLeftMargin(0.25);
  gStyle->SetPadRightMargin(0.12);

  std::vector<TH2F*> h2vec;
  std::vector<std::string> dirnames;

  TFile *f = new TFile(filename);
  
  f->cd();
  TDirectory *stardir = gDirectory;
  TObject *thesourcedir;
  TIter nextdir(stardir->GetListOfKeys());
  UInt_t theDirNumber(0);

  while((thesourcedir=nextdir())){
    
    theDirNumber++;
    TString dirName = thesourcedir->GetName();
    std::string dirname_s=dirName.Data();
    dirnames.push_back(dirname_s);

    //stardir->cd(dirName);
    //TDirectory *current_sourcedir = gDirectory;
    //TH1::AddDirectory(kFALSE);

    std::cout << "*************************" <<std::endl;
    std::cout << "Reading Directory: " << dirName <<std::endl;
    
    TH2F *h2TCHE  = (TH2F*)f->Get(dirName+"/h2MigrationMatrixTCHE");
    TH2F *h2TCHP  = (TH2F*)f->Get(dirName+"/h2MigrationMatrixTCHP");
    TH2F *h2SSVHE = (TH2F*)f->Get(dirName+"/h2MigrationMatrixSSVHE");
    TH2F *h2SSVHP = (TH2F*)f->Get(dirName+"/h2MigrationMatrixSSVHP");
    TH2F *h2CSV   = (TH2F*)f->Get(dirName+"/h2MigrationMatrixCSV");
    TH2F *h2JP    = (TH2F*)f->Get(dirName+"/h2MigrationMatrixJP");
    TH2F *h2JBP   = (TH2F*)f->Get(dirName+"/h2MigrationMatrixJBP");
   
    resetHistoName(h2TCHE,dirName);
    resetHistoName(h2TCHP,dirName);
    resetHistoName(h2SSVHE,dirName);
    resetHistoName(h2SSVHP,dirName);
    resetHistoName(h2CSV,dirName);
    resetHistoName(h2JP,dirName);
    resetHistoName(h2JBP,dirName);   

    h2vec.push_back(h2TCHE);
    h2vec.push_back(h2TCHP);
    h2vec.push_back(h2SSVHE);
    h2vec.push_back(h2SSVHP);
    h2vec.push_back(h2CSV);
    h2vec.push_back(h2JP);
    h2vec.push_back(h2JBP);
    
  }
  
  for(UInt_t num=0; num<theDirNumber; num++){

    TCanvas *c1 = new TCanvas(Form("defCanvas_%d",num),Form("defCanvas_%d",num),1200,600);
    c1->Divide(4,2);
    
    Int_t histcount(0);
    for(UInt_t h=0+(num*7); h<7+(num*7); h++){
      c1->cd(histcount+1);
      SetHistoStyle(h2vec[h]);
      h2vec[h]->Scale(100./h2vec[h]->GetEntries());
      h2vec[h]->Draw("colz");
      h2vec[h]->Draw("TEXTsame");
      histcount++;
    }
    
    c1->Update(); 

    std::cout<<"dirnames["<<num<<"]= "<<dirnames[num]<<std::endl;
    c1->SaveAs(Form("matrices_%s.png",dirnames[num].c_str()));

  }
}




