/* 
 Macro to make PV Validation plots

 usage:
 $ root -l
 $ root [0] .L PlotPVValidation.C++
 $ root [1] PlotPVValidation("filename1.root=legendentry1,filename2.root=legendentry2", .... , number_of_files,"YourCut") 

 - YourCut is a TCut (default ="")
 
 // M. Musich - INFN Torino

*/

#include <string>
#include <sstream>
#include <vector>
#include <Riostream.h>
#include <iomanip>
#include "TFile.h"
#include "TPaveStats.h"
#include "TList.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLegendEntry.h"
#include <TPaveText.h>
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TVectorD.h"
#include <string>
#include <sstream>
#include <vector>
#include <Riostream.h>
#include <iomanip>
#include "TFile.h"
#include "TPaveStats.h"
#include "THStack.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

TList *FileList;
TList *LabelList;
TFile *Target;

void PlotPVValidation(TString namesandlabels,Int_t nOfFiles, TCut theCut, Bool_t useBS);
void fitResiduals(TH1 *hist);
void DoubleGausFitResiduals(TH1 *hist);

Double_t fULine(Double_t *x, Double_t *par);
Double_t fDLine(Double_t *x, Double_t *par);
void FitULine(TH1 *hist);
void FitDLine(TH1 *hist);

void makeNiceTF1Style(TF1 *f1,int color);
void MakeNiceHistoStyle(TH1 *hist);
void MakeNiceTrendPlotStyle(TH1 *hist);

void PlotPVValidation(TString namesandlabels,Int_t nOfFiles, TCut theCut,Bool_t useBS){

  TCut theDefaultCut= "hasRecVertex==1&&isGoodTrack==1";

  ofstream outfile("FittedDeltaZ.txt");

  TH1::StatOverflows(kTRUE);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat("e");
  gStyle->SetPadTopMargin(0.12);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadBorderMode(0);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.05);///---> gStyle->SetStatFontSize(0.025);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptFit(1);
  gStyle->SetNdivisions(510);

  using namespace std;
  using namespace TMath;

  FileList = new TList();
  LabelList = new TList();
  
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
    cout<<"&MSG-e: NOfFiles = "<<nOfFiles<<" != "<<NOfFiles<<" given as input"<<endl;
    return;
  }  
  
  TTree *trees[nOfFiles]; 
  TString LegLabels[nOfFiles];  

  for(Int_t j=0; j < nOfFiles; j++) {
    
    // Retrieve files

    TFile *fin = (TFile*)FileList->At(j);    
    trees[j] = (TTree*)fin->Get("tree");

    // Retrieve labels

    TObjString* legend = (TObjString*)LabelList->At(j);
    LegLabels[j] = legend->String();
    cout<<"LegLabels["<<j<<"]"<<LegLabels[j]<<endl;
    
  }

  TString label = LegLabels[nOfFiles-1];
  TFile fileout("histos"+label+".root", "new");

  //#################### Canvases #####################

  Int_t nplots = 12; 

  TCanvas *c1 = new TCanvas("c1","c1",800,550);
  c1->SetFillColor(10);  
  
  TCanvas *c2 = new TCanvas("c2","c2",1460,600);
  c2->SetFillColor(10);  
  c2->Divide(nplots/2,2);
  
  TCanvas *c3 = new TCanvas("c3","c3",1460,600);
  c3->SetFillColor(10);  
  c3->Divide(nplots/2,2);

  TCanvas *c2Eta = new TCanvas("c2Eta","c2Eta",1460,600);
  c2Eta->SetFillColor(10);  
  c2Eta->Divide(nplots/2,2);
  
  TCanvas *c3Eta = new TCanvas("c3Eta","c3Eta",1460,600);
  c3Eta->SetFillColor(10);  
  c3Eta->Divide(nplots/2,2);
   
  TCanvas *c4 = new TCanvas("c4","c4",1460,550);
  c4->SetFillColor(10);  
  c4->Divide(2,1);
  c4->cd(1)->SetBottomMargin(0.12);
  c4->cd(1)->SetLeftMargin(0.13);
  c4->cd(2)->SetBottomMargin(0.12);
  c4->cd(2)->SetLeftMargin(0.13);

  TCanvas *c4monly  = new TCanvas("c4monly","c4monly",730,550);
  c4monly->cd()->SetBottomMargin(0.12);
  c4monly->cd()->SetLeftMargin(0.13);
  c4monly->SetFillColor(10); 
  
  TCanvas *c5 = new TCanvas("c5","c5",1460,550);
  c5->SetFillColor(10);  
  c5->Divide(2,1);
  c5->cd(1)->SetBottomMargin(0.12);
  c5->cd(1)->SetLeftMargin(0.13);
  c5->cd(2)->SetBottomMargin(0.12);
  c5->cd(2)->SetLeftMargin(0.13);

  TCanvas *c5monly  = new TCanvas("c5monly","c5monly",730,550);
  c5monly->cd()->SetBottomMargin(0.12);
  c5monly->cd()->SetLeftMargin(0.13);
  c5monly->SetFillColor(10); 

  TCanvas *c4Eta = new TCanvas("c4Eta","c4Eta",1460,550);
  c4Eta->SetFillColor(10);  
  c4Eta->Divide(2,1);
  c4Eta->cd(1)->SetBottomMargin(0.12);
  c4Eta->cd(1)->SetLeftMargin(0.13);
  c4Eta->cd(2)->SetBottomMargin(0.12);
  c4Eta->cd(2)->SetLeftMargin(0.13);
  
  TCanvas *c5Eta = new TCanvas("c5Eta","c5Eta",1460,550);
  c5Eta->SetFillColor(10); 
  c5Eta->Divide(2,1);
  c5Eta->cd(1)->SetBottomMargin(0.12);
  c5Eta->cd(1)->SetLeftMargin(0.13);
  c5Eta->cd(2)->SetBottomMargin(0.12);
  c5Eta->cd(2)->SetLeftMargin(0.13);
               
  TLegend *lego = new TLegend(0.15,0.75,0.45,0.89);
  lego->SetFillColor(10);
  
  Int_t colors[8]={1,2,4,6,3,8,9,5};
  //Int_t colors[8]={6,2,4,1,3,8,9,5};
  Int_t styles[8]={20,4,21,22,25,28,4,29}; 

  TH1F *histosdxy[nOfFiles][nplots];
  TH1F *histosdz[nOfFiles][nplots];

  TH1F *histosEtadxy[nOfFiles][nplots];
  TH1F *histosEtadz[nOfFiles][nplots];

  //############################
  // Canvas of bin residuals
  //############################

  float phipitch = (2*TMath::Pi())/nplots;
  float etapitch = 5./nplots;
  //cout<<"etapitch= "<<etapitch;

  int nbins_ = 100;

  for(Int_t i=0; i<nplots; i++){
    for(Int_t j=0; j < nOfFiles; j++){

      //for 12 bins
      //float phiF=(-TMath::Pi()+i*0.5235)*(180/TMath::Pi());
      //float phiL=(-TMath::Pi()+(i+1)*0.5235)*(180/TMath::Pi());
      //for 6 bins
      //float phiF=(-TMath::Pi()+i*1.047)*(180/TMath::Pi());
      //float phiL=(-TMath::Pi()+(i+1)*1.047)*(180/TMath::Pi());

      //for general number of bins
  
      float phiF=(-TMath::Pi()+i*phipitch)*(180/TMath::Pi());
      float phiL=(-TMath::Pi()+(i+1)*phipitch)*(180/TMath::Pi());

      float etaF=-2.5+i*etapitch;
      float etaL=-2.5+(i+1)*etapitch;
      
      float dxymax_phi,dzmax_phi, dxymax_eta, dzmax_eta;
      
      if(useBS){
	dxymax_phi=1500;
	dzmax_phi=300000;
	dxymax_eta=1500;
	dzmax_eta=300000;
      } else {
       	dxymax_phi=1500;
	dzmax_phi=1500;
	dxymax_eta=1500;
	dzmax_eta=1500;
      }

      //###################### vs phi #####################

      char histotitledxy[129];
      sprintf(histotitledxy,"%.2f #circ<#phi<%.2f #circ;d_{xy} [#mum];tracks",phiF,phiL);
      char treeleg1[129];
      sprintf(treeleg1,"histodxyfile%i_plot%i",j,i);
      histosdxy[j][i] = new TH1F(treeleg1,histotitledxy,nbins_,-dxymax_phi,dxymax_phi);
      histosdxy[j][i]->SetLineColor(colors[j]);
      MakeNiceHistoStyle(histosdxy[j][i]);
      
      char histotitledz[129];
      sprintf(histotitledz,"%.2f #circ<#phi<%.2f #circ;d_{z} [#mum];tracks",phiF,phiL);
      char treeleg2[129];
      sprintf(treeleg2,"histodzfile%i_plot%i",j,i);
      histosdz[j][i] = new TH1F(treeleg2,histotitledz,nbins_,-dzmax_phi,dzmax_phi);
      histosdz[j][i]->SetLineColor(colors[j]);
      MakeNiceHistoStyle(histosdz[j][i]);

       //###################### vs eta ##################### 

      char histotitledxyeta[129];
      sprintf(histotitledxyeta,"%.2f <#eta<%.2f ;d_{xy} [#mum];tracks",etaF,etaL);
      char treeleg1eta[129];
      sprintf(treeleg1eta,"histodxyEtafile%i_plot%i",j,i);
      histosEtadxy[j][i] = new TH1F(treeleg1eta,histotitledxyeta,nbins_,-dxymax_eta,dxymax_eta);
      histosEtadxy[j][i]->SetLineColor(colors[j]);
      MakeNiceHistoStyle(histosEtadxy[j][i]);
      
      char histotitledzeta[129];
      sprintf(histotitledzeta,"%.2f <#eta<%.2f ;d_{z} [#mum];tracks",etaF,etaL);
      char treeleg2eta[129];
      sprintf(treeleg2eta,"histodzEtafile%i_plot%i",j,i);
      histosEtadz[j][i] = new TH1F(treeleg2eta,histotitledzeta,nbins_,-dzmax_eta,dzmax_eta);
      histosEtadz[j][i]->SetLineColor(colors[j]);
      MakeNiceHistoStyle(histosEtadz[j][i]);
    }
  }
  
  cout<<"Checkpoint 1: before anything"<<endl;

  //################ Phi #####################

  TObject *statObjsdxy[nOfFiles][nplots];
  TPaveStats *statsdxy[nOfFiles][nplots];
  TObject *statObjsdz[nOfFiles][nplots];
  TPaveStats *statsdz[nOfFiles][nplots];
  
  TH1F *histoMeansdxy[nOfFiles];  
  TH1F *histoMeansdz[nOfFiles];     
  TH1F *histoSigmasdxy[nOfFiles];   
  TH1F *histoSigmasdz[nOfFiles];    
  
  TF1 *tmpsdxy[nOfFiles][nplots];
  TF1 *tmpsdz[nOfFiles][nplots];
  
  TF1 *FitDzUp[nOfFiles];
  TF1 *FitDzDown[nOfFiles];
  
  TF1 *fleft[nOfFiles]; 
  TF1 *fright[nOfFiles];
  TF1 *fall[nOfFiles];  

  Double_t meansdxy[nOfFiles][nplots];
  Double_t meansdxyError[nOfFiles][nplots];
  //Double_t mediansdxy[nOfFiles][nplots];
  Double_t sigmasdxy[nOfFiles][nplots];
  Double_t sigmasdxyError[nOfFiles][nplots];

  Double_t meansdz[nOfFiles][nplots];
  Double_t meansdzError[nOfFiles][nplots];
  //Double_t mediansdz[nOfFiles][nplots];
  Double_t sigmasdz[nOfFiles][nplots];
  Double_t sigmasdzError[nOfFiles][nplots];

  //################ Eta #####################

  TObject *statObjsdxyEta[nOfFiles][nplots];
  TPaveStats *statsdxyEta[nOfFiles][nplots];
  TObject *statObjsdzEta[nOfFiles][nplots];
  TPaveStats *statsdzEta[nOfFiles][nplots]; 

  TH1F *histoMeansdxyEta[nOfFiles];  
  TH1F *histoMeansdzEta[nOfFiles];     
  TH1F *histoSigmasdxyEta[nOfFiles];   
  TH1F *histoSigmasdzEta[nOfFiles];    
  
  TF1 *tmpsdxyEta[nOfFiles][nplots];
  TF1 *tmpsdzEta[nOfFiles][nplots];
  
  Double_t meansdxyEta[nOfFiles][nplots];
  Double_t meansdxyEtaError[nOfFiles][nplots];
  //Double_t mediansdxy[nOfFiles][nplots];
  Double_t sigmasdxyEta[nOfFiles][nplots];
  Double_t sigmasdxyEtaError[nOfFiles][nplots];

  Double_t meansdzEta[nOfFiles][nplots];
  Double_t meansdzEtaError[nOfFiles][nplots];
  //Double_t mediansdz[nOfFiles][nplots];
  Double_t sigmasdzEta[nOfFiles][nplots];
  Double_t sigmasdzEtaError[nOfFiles][nplots];

  Double_t highedge=nplots-0.5;
  Double_t lowedge=-0.5;

  //#################### initialization limits ################################
  
  Double_t dxymean_phi, dxymean_eta, dxysigma_phi, dxysigma_eta, dzmean_phi, dzmean_eta, dzsigma_phi, dzsigma_eta ;

  if(useBS){
    dxymean_phi=100; 
    dxymean_eta=40; 
    dxysigma_phi=310; 
    dxysigma_eta=510; 
    dzmean_phi=2000; 
    dzmean_eta=2000; 
    dzsigma_phi=100000; 
    dzsigma_eta=100000;
  } else {
    dxymean_phi=200; 
    dxymean_eta=40; 
    dxysigma_phi=310; 
    dxysigma_eta=310; 
    dzmean_phi=150; 
    dzmean_eta=160; 
    dzsigma_phi=510; 
    dzsigma_eta=1000;
  }
  
  for(Int_t j=0; j < nOfFiles; j++) { 

    //################ phi #########################

    char meanslegdxy[129];
    sprintf(meanslegdxy,"histomeansdxy_file%i",j);
    histoMeansdxy[j]  = new TH1F(meanslegdxy,"#LT d_{xy} #GT vs #phi sector",nplots,lowedge,highedge);
    histoMeansdxy[j]->SetLineColor(colors[j]);
    histoMeansdxy[j]->GetYaxis()->SetRangeUser(-dxymean_phi,dxymean_phi);
    histoMeansdxy[j]->SetMarkerStyle(styles[j]);
    histoMeansdxy[j]->SetMarkerColor(colors[j]);
    histoMeansdxy[j]->GetXaxis()->SetTitle("#phi (sector) [degrees]");
    histoMeansdxy[j]->GetYaxis()->SetTitle ("#LT d_{xy} #GT [#mum]");
    MakeNiceTrendPlotStyle(histoMeansdxy[j]);
   
    char meanslegdz[129];
    sprintf(meanslegdz,"histomeansdz_file%i",j);
    histoMeansdz[j]   = new TH1F(meanslegdz,"#LT d_{z} #GT vs #phi sector",nplots,lowedge,highedge); 
    histoMeansdz[j]->SetLineColor(colors[j]);
    histoMeansdz[j]->GetYaxis()->SetRangeUser(-dzmean_phi,dzmean_phi);
    histoMeansdz[j]->SetMarkerStyle(styles[j]);
    histoMeansdz[j]->SetMarkerColor(colors[j]);
    histoMeansdz[j]->GetXaxis()->SetTitle("#phi (sector) [degrees]");
    histoMeansdz[j]->GetYaxis()->SetTitle ("#LT d_{z} #GT [#mum]");
    MakeNiceTrendPlotStyle(histoMeansdz[j]);

    char sigmaslegdxy[129];
    sprintf(sigmaslegdxy,"histosigmasdxy_file%i",j);
    histoSigmasdxy[j] = new TH1F(sigmaslegdxy,"#sigma_{d_{xy}} vs #phi sector",nplots,lowedge,highedge);
    histoSigmasdxy[j]->SetLineColor(colors[j]);
    histoSigmasdxy[j]->GetYaxis()->SetRangeUser(0,dxysigma_phi);
    histoSigmasdxy[j]->SetMarkerStyle(styles[j]);
    histoSigmasdxy[j]->SetMarkerColor(colors[j]);
    histoSigmasdxy[j]->GetXaxis()->SetTitle("#phi (sector) [degrees]");
    histoSigmasdxy[j]->GetYaxis()->SetTitle ("#sigma_{d_{xy}} [#mum]");
    MakeNiceTrendPlotStyle(histoSigmasdxy[j]);

    char sigmaslegdz[129];
    sprintf(sigmaslegdz,"histosigmasdz_file%i",j);
    histoSigmasdz[j]  = new TH1F(sigmaslegdz,"#sigma_{d_{z}} vs #phi sector",nplots,lowedge,highedge);
    histoSigmasdz[j]->SetLineColor(colors[j]);
    histoSigmasdz[j]->GetYaxis()->SetRangeUser(0,dzsigma_phi);
    histoSigmasdz[j]->SetMarkerStyle(styles[j]);
    histoSigmasdz[j]->SetMarkerColor(colors[j]);
    histoSigmasdz[j]->GetXaxis()->SetTitle("#phi (sector) [degrees]");
    histoSigmasdz[j]->GetYaxis()->SetTitle ("#sigma_{d_{z}} [#mum]");
    MakeNiceTrendPlotStyle(histoSigmasdz[j]);

    //################ eta #########################
  
    char meanslegdxyEta[129];
    sprintf(meanslegdxyEta,"histomeansEtadxy_file%i",j);
    histoMeansdxyEta[j]  = new TH1F(meanslegdxyEta,"#LT d_{xy} #GT vs #eta sector",nplots,lowedge,highedge);
    histoMeansdxyEta[j]->SetLineColor(colors[j]);
    histoMeansdxyEta[j]->GetYaxis()->SetRangeUser(-dxymean_eta,dxymean_eta);
    histoMeansdxyEta[j]->SetMarkerStyle(styles[j]);
    histoMeansdxyEta[j]->SetMarkerColor(colors[j]);
    histoMeansdxyEta[j]->GetXaxis()->SetTitle("#eta (sector)");
    histoMeansdxyEta[j]->GetYaxis()->SetTitle ("#LT d_{xy} #GT [#mum]");
    MakeNiceTrendPlotStyle(histoMeansdxyEta[j]);
   
    char meanslegdzEta[129];
    sprintf(meanslegdzEta,"histomeansEtadz_file%i",j);
    histoMeansdzEta[j]   = new TH1F(meanslegdzEta,"#LT d_{z} #GT vs #eta sector",nplots,lowedge,highedge); 
    histoMeansdzEta[j]->SetLineColor(colors[j]);
    histoMeansdzEta[j]->GetYaxis()->SetRangeUser(-dzmean_eta,dzmean_eta);
    histoMeansdzEta[j]->SetMarkerStyle(styles[j]);
    histoMeansdzEta[j]->SetMarkerColor(colors[j]);
    histoMeansdzEta[j]->GetXaxis()->SetTitle("#eta (sector)");
    histoMeansdzEta[j]->GetYaxis()->SetTitle ("#LT d_{z} #GT [#mum]");
    MakeNiceTrendPlotStyle(histoMeansdzEta[j]);
   
    char sigmaslegdxyEta[129];
    sprintf(sigmaslegdxyEta,"histosigmasEtadxy_file%i",j);
    histoSigmasdxyEta[j] = new TH1F(sigmaslegdxyEta,"#sigma_{d_{xy}} vs #eta sector",nplots,lowedge,highedge);
    histoSigmasdxyEta[j]->SetLineColor(colors[j]);
    histoSigmasdxyEta[j]->GetYaxis()->SetRangeUser(0,dxysigma_eta);
    histoSigmasdxyEta[j]->SetMarkerStyle(styles[j]);
    histoSigmasdxyEta[j]->SetMarkerColor(colors[j]);
    histoSigmasdxyEta[j]->GetXaxis()->SetTitle("#eta (sector)");
    histoSigmasdxyEta[j]->GetYaxis()->SetTitle ("#sigma_{d_{xy}} [#mum]");
    MakeNiceTrendPlotStyle(histoSigmasdxyEta[j]);

    char sigmaslegdzEta[129];
    sprintf(sigmaslegdzEta,"histosigmasEtadz_file%i",j);
    histoSigmasdzEta[j]  = new TH1F(sigmaslegdzEta,"#sigma_{d_{z}} vs #eta sector",nplots,lowedge,highedge);
    histoSigmasdzEta[j]->SetLineColor(colors[j]);
    histoSigmasdzEta[j]->GetYaxis()->SetRangeUser(0,dzsigma_eta);
    histoSigmasdzEta[j]->SetMarkerStyle(styles[j]);
    histoSigmasdzEta[j]->SetMarkerColor(colors[j]);
    histoSigmasdzEta[j]->GetXaxis()->SetTitle("#eta (sector)");
    histoSigmasdzEta[j]->GetYaxis()->SetTitle ("#sigma_{d_{z}} [#mum]");
    MakeNiceTrendPlotStyle(histoSigmasdzEta[j]);
  }

   cout<<"Checkpoint 2: starts fitting xy residuals"<<endl;

   //Histos maxima
   
   double hmaxDxyphi=0.;
   double hmaxDxyeta=0.;
   double hmaxDzphi=0.;
   double hmaxDzeta=0.;

   //######################## dxy Residuals #################################
   
   for(Int_t j=0; j < nOfFiles; j++) { 
     for(Int_t i=0; i<nplots; i++){
       
       //################### phi ########################
       
       char thePhiCutstring[128];
       sprintf(thePhiCutstring,"(phi>-TMath::Pi()+%.1i*(%.3f))&&(phi<-TMath::Pi()+%.1i*(%.3f))",i,phipitch,i+1,phipitch);
       TCut thePhiCut = thePhiCutstring;
       //cout<<thePhiCutstring<<endl;
       
       char treeleg3[129];

       if(useBS){
	 sprintf(treeleg3,"dxyBs*10000>>histodxyfile%i_plot%i",j,i);
       } else {
	 sprintf(treeleg3,"dxyFromMyVertex*10000>>histodxyfile%i_plot%i",j,i);
       }

       char phipositionString[129];
       float_t phiInterval = (360.)/nplots;
       float phiposition = (-180+i*phiInterval)+(phiInterval/2);
       sprintf(phipositionString,"%.f",phiposition);
       //cout<<"phiposition: "<<phiposition<<" phipositionString: "<<phipositionString<<endl;


       c2->cd(i+1);
       if(j==0){
	 trees[j]->Draw(treeleg3,theCut&&theDefaultCut&&thePhiCut);  
	 hmaxDxyphi=histosdxy[0][i]->GetMaximum();
       }
       else{ 
	 trees[j]->Draw(treeleg3,theCut&&theDefaultCut&&thePhiCut,"sames"); 
	 
	 if(histosdxy[j][i]->GetMaximum() >  hmaxDxyphi){
	   hmaxDxyphi = histosdxy[j][i]->GetMaximum();
	 }
       }
       
       histosdxy[0][i]->GetYaxis()->SetRangeUser(0.01,hmaxDxyphi*1.10);
       histosdxy[0][i]->Draw("sames");
       
       if(j!=0){
	 histosdxy[j][i]->Draw("sames");
       }

       fitResiduals(histosdxy[j][i]);   
       tmpsdxy[j][i] = (TF1*)histosdxy[j][i]->GetListOfFunctions()->FindObject("tmp");   
       if(tmpsdxy[j][i])  {
	 tmpsdxy[j][i]->SetLineColor(colors[j]);
	 tmpsdxy[j][i]->SetLineWidth(2); 
	 tmpsdxy[j][i]->Draw("same"); 
       

       meansdxy[j][i]=(tmpsdxy[j][i])->GetParameter(1);
       sigmasdxy[j][i]=(tmpsdxy[j][i])->GetParameter(2);
       
       meansdxyError[j][i]=(tmpsdxy[j][i])->GetParError(1);
       sigmasdxyError[j][i]=(tmpsdxy[j][i])->GetParError(2); 
      
       histoMeansdxy[j]->SetBinContent(i+1,meansdxy[j][i]); 
       histoMeansdxy[j]->SetBinError(i+1,meansdxyError[j][i]); 
       histoMeansdxy[j]->GetXaxis()->SetBinLabel(i+1,phipositionString);
       histoSigmasdxy[j]->SetBinContent(i+1,sigmasdxy[j][i]);
       histoSigmasdxy[j]->SetBinError(i+1,sigmasdxyError[j][i]); 
       histoSigmasdxy[j]->GetXaxis()->SetBinLabel(i+1,phipositionString);
       }

       else{
	 histoMeansdxy[j]->SetBinContent(i+1,0.); 
	 histoMeansdxy[j]->SetBinError(i+1,0.);
	 histoMeansdxy[j]->GetXaxis()->SetBinLabel(i+1,phipositionString);
	 histoSigmasdxy[j]->SetBinContent(i+1,0.);
	 histoSigmasdxy[j]->SetBinError(i+1,0.); 
	 histoSigmasdxy[j]->GetXaxis()->SetBinLabel(i+1,phipositionString); 
       }
       

       c2->Draw("sames");   
       c2->cd(i+1);
       
       statObjsdxy[j][i] = histosdxy[j][i]->GetListOfFunctions()->FindObject("stats");
       statsdxy[j][i] = static_cast<TPaveStats*>(statObjsdxy[j][i]);
       statsdxy[j][i]->SetLineColor(colors[j]);
       statsdxy[j][i]->SetTextColor(colors[j]);
       statsdxy[j][i]->SetShadowColor(10);
       statsdxy[j][i]->SetFillColor(10);
       statsdxy[j][i]->SetX1NDC(0.66);
       statsdxy[j][i]->SetY1NDC(0.80-j*0.19);
       statsdxy[j][i]->SetX2NDC(0.99);
       statsdxy[j][i]->SetY2NDC(0.99-j*0.19);
       statsdxy[j][i]->Draw("same");   
  
       //################### eta ########################
       
       char theEtaCutstring[128];
       sprintf(theEtaCutstring,"(eta>-2.5+%.1i*(%.3f))&&(eta<-2.5+%.1i*(%.3f))",i,etapitch,i+1,etapitch);
       TCut theEtaCut = theEtaCutstring;
       //cout<<theEtaCutstring<<endl;
       
       char treeleg3Eta[129];
       if(useBS){
	 sprintf(treeleg3Eta,"dxyBs*10000>>histodxyEtafile%i_plot%i",j,i);
       } else {
	 sprintf(treeleg3Eta,"dxyFromMyVertex*10000>>histodxyEtafile%i_plot%i",j,i);
       }
       
       char etapositionString[129];
       float etaposition = (-2.5+i*0.417)+(0.417/2);
       sprintf(etapositionString,"%.1f",etaposition);
       //cout<<"etaposition: "<<etaposition<<" etapositionString: "<<etapositionString<<endl;
       
       c2Eta->cd(i+1);
       if(j==0){
	 trees[j]->Draw(treeleg3Eta,theCut&&theDefaultCut&&theEtaCut);
	 hmaxDxyeta=histosEtadxy[0][i]->GetMaximum();  
       }
       else{ 
	 trees[j]->Draw(treeleg3Eta,theCut&&theDefaultCut&&theEtaCut,"sames"); 
	 if(histosEtadxy[j][i]->GetMaximum() >  hmaxDxyeta){
	   hmaxDxyeta = histosEtadxy[j][i]->GetMaximum();
	 }
       }

       histosEtadxy[0][i]->GetYaxis()->SetRangeUser(0.01,hmaxDxyeta*1.10);
       histosEtadxy[0][i]->Draw("sames");
       
       if(j!=0){
	 histosEtadxy[j][i]->Draw("sames");
       }

       
       fitResiduals(histosEtadxy[j][i]);   
       tmpsdxyEta[j][i] = (TF1*)histosEtadxy[j][i]->GetListOfFunctions()->FindObject("tmp");   
       if(tmpsdxyEta[j][i])  {
	 tmpsdxyEta[j][i]->SetLineColor(colors[j]);
	 tmpsdxyEta[j][i]->SetLineWidth(2); 
       
	 
       meansdxyEta[j][i]=(tmpsdxyEta[j][i])->GetParameter(1);
       sigmasdxyEta[j][i]=(tmpsdxyEta[j][i])->GetParameter(2);
       
       meansdxyEtaError[j][i]=(tmpsdxyEta[j][i])->GetParError(1);
       sigmasdxyEtaError[j][i]=(tmpsdxyEta[j][i])->GetParError(2);
       
       histoMeansdxyEta[j]->SetBinContent(i+1,meansdxyEta[j][i]); 
       histoMeansdxyEta[j]->SetBinError(i+1,meansdxyEtaError[j][i]); 
       histoMeansdxyEta[j]->GetXaxis()->SetBinLabel(i+1,etapositionString);
       histoSigmasdxyEta[j]->SetBinContent(i+1,sigmasdxyEta[j][i]);
       histoSigmasdxyEta[j]->SetBinError(i+1,sigmasdxyEtaError[j][i]); 
       histoSigmasdxyEta[j]->GetXaxis()->SetBinLabel(i+1,etapositionString);

       }
       else{

	 histoMeansdxyEta[j]->SetBinContent(i+1,0.); 
	 histoMeansdxyEta[j]->SetBinError(i+1,0.); 	 
	 histoMeansdxyEta[j]->GetXaxis()->SetBinLabel(i+1,etapositionString);
	 histoSigmasdxyEta[j]->SetBinContent(i+1,0.);
	 histoSigmasdxyEta[j]->SetBinError(i+1,0.);
	 histoSigmasdxyEta[j]->GetXaxis()->SetBinLabel(i+1,etapositionString);
       }

       c2Eta->Draw("sames");   
       c2Eta->cd(i+1);
       
       statObjsdxyEta[j][i] = histosEtadxy[j][i]->GetListOfFunctions()->FindObject("stats");
       statsdxyEta[j][i] = static_cast<TPaveStats*>(statObjsdxyEta[j][i]);
       statsdxyEta[j][i]->SetLineColor(colors[j]);
       statsdxyEta[j][i]->SetTextColor(colors[j]);
       statsdxyEta[j][i]->SetShadowColor(10);
       statsdxyEta[j][i]->SetFillColor(10);
       statsdxyEta[j][i]->SetX1NDC(0.66);
       statsdxyEta[j][i]->SetY1NDC(0.80-j*0.18);
       statsdxyEta[j][i]->SetX2NDC(0.99);
       statsdxyEta[j][i]->SetY2NDC(0.99-j*0.19);
       statsdxyEta[j][i]->Draw("same"); 
     } 
   }
   
   TString Canvas2Title ="VertexDxyResidualsPhiBin_"+label;
   TString Canvas2formatpng=Canvas2Title+".png"; 
   TString Canvas2formateps=Canvas2Title+".eps"; 
   c2->SaveAs(Canvas2formatpng);
   c2->Write();
   //c2->SaveAs(Canvas2formateps);

   TString Canvas2TitleEta ="VertexDxyResidualsEtaBin_"+label;
   TString Canvas2formatEtapng=Canvas2TitleEta+".png"; 
   TString Canvas2formatEtaeps=Canvas2TitleEta+".eps"; 
   c2Eta->SaveAs(Canvas2formatEtapng);
   c2Eta->Write();
   //c2Eta->SaveAs(Canvas2formatEtaeps);
   
   cout<<"Checkpoint 3: starts fitting z residuals"<<endl;
 
   //######################## dz Residuals #################################

   for(Int_t j=0; j < nOfFiles; j++) { 
     for(Int_t i=0; i<nplots; i++){
    
      char thePhiCutstring[128];
      sprintf(thePhiCutstring,"(phi>-TMath::Pi()+%.1i*(%.3f))&&(phi<-TMath::Pi()+%.1i*(%.3f))",i,phipitch,i+1,phipitch);
      TCut thePhiCut = thePhiCutstring;
      //cout<<thePhiCutstring<<endl;
      
      char treeleg4[129];
      if(useBS){
	sprintf(treeleg4,"dzBs*10000>>histodzfile%i_plot%i",j,i);
      } else {
	sprintf(treeleg4,"dzFromMyVertex*10000>>histodzfile%i_plot%i",j,i);
      }

      char phipositionString[129];
      float_t phiInterval = (360.)/nplots;
      float phiposition = (-180+i*phiInterval)+(phiInterval/2);
      sprintf(phipositionString,"%.f",phiposition);
      //cout<<"phiposition: "<<phiposition<<" phipositionString: "<<phipositionString<<endl;
	
      c3->cd(i+1);
      if(j==0){
	trees[j]->Draw(treeleg4,theCut&&theDefaultCut&&thePhiCut);
	hmaxDzphi=histosdz[0][i]->GetMaximum();  
      }
      
      else{ 
	trees[j]->Draw(treeleg4,theCut&&theDefaultCut&&thePhiCut,"sames"); 
	if(histosdz[j][i]->GetMaximum() >  hmaxDzphi){
	  hmaxDzphi = histosdz[j][i]->GetMaximum();
	}
      }

      histosdz[0][i]->GetYaxis()->SetRangeUser(0.01,hmaxDzphi*1.10);
      histosdz[0][i]->Draw("sames");
       
      if(j!=0){
	histosdz[j][i]->Draw("sames");
      }
      
      fitResiduals(histosdz[j][i]);   
      tmpsdz[j][i] = (TF1*)histosdz[j][i]->GetListOfFunctions()->FindObject("tmp");   
      if(tmpsdz[j][i])  {
	tmpsdz[j][i]->SetLineColor(colors[j]);
	tmpsdz[j][i]->SetLineWidth(2);
	tmpsdz[j][i]->Draw("same");  
      
        
      meansdz[j][i]=(tmpsdz[j][i])->GetParameter(1);
      sigmasdz[j][i]=(tmpsdz[j][i])->GetParameter(2);

      meansdzError[j][i]=(tmpsdz[j][i])->GetParError(1);
      sigmasdzError[j][i]=(tmpsdz[j][i])->GetParError(2);

      histoMeansdz[j]->SetBinContent(i+1,meansdz[j][i]); 
      histoMeansdz[j]->SetBinError(i+1,meansdzError[j][i]); 
      histoMeansdz[j]->GetXaxis()->SetBinLabel(i+1,phipositionString);
      histoSigmasdz[j]->SetBinContent(i+1,sigmasdz[j][i]);
      histoSigmasdz[j]->SetBinError(i+1,sigmasdzError[j][i]);
      histoSigmasdz[j]->GetXaxis()->SetBinLabel(i+1,phipositionString);
      }

      else{
	histoMeansdz[j]->SetBinContent(i+1,0.); 
	histoMeansdz[j]->SetBinError(i+1,0.); 
	histoMeansdz[j]->GetXaxis()->SetBinLabel(i+1,phipositionString);
	histoSigmasdz[j]->SetBinContent(i+1,0.);
	histoSigmasdz[j]->SetBinError(i+1,0.);
	histoSigmasdz[j]->GetXaxis()->SetBinLabel(i+1,phipositionString);
      }

      c3->Draw("sames");
      c3->cd(i+1);
	       
      statObjsdz[j][i] = histosdz[j][i]->GetListOfFunctions()->FindObject("stats");
      statsdz[j][i] = static_cast<TPaveStats*>(statObjsdz[j][i]);     
      statsdz[j][i]->SetLineColor(colors[j]);
      statsdz[j][i]->SetTextColor(colors[j]);
      statsdz[j][i]->SetFillColor(10);
      statsdz[j][i]->SetShadowColor(10);
      statsdz[j][i]->SetX1NDC(0.66);
      statsdz[j][i]->SetY1NDC(0.80-j*0.18);
      statsdz[j][i]->SetX2NDC(0.99);
      statsdz[j][i]->SetY2NDC(0.99-j*0.19);
      statsdz[j][i]->Draw("same");

      //################### eta ########################
      
      char theEtaCutstring[128];
      sprintf(theEtaCutstring,"(eta>-2.5+%.1i*(%.3f))&&(eta<-2.5+%.1i*(%.3f))",i,etapitch,i+1,etapitch);
      TCut theEtaCut = theEtaCutstring;
      //cout<<theEtaCutstring<<endl;
      
      char treeleg3Eta[129];
      if(useBS){
	sprintf(treeleg3Eta,"dzBs*10000>>histodzEtafile%i_plot%i",j,i);	
      } else {
	sprintf(treeleg3Eta,"dzFromMyVertex*10000>>histodzEtafile%i_plot%i",j,i);
      }
      
      char etapositionString[129];
      float etaposition = (-2.5+i*0.417)+(0.417/2);
      sprintf(etapositionString,"%.1f",etaposition);
      //cout<<"etaposition: "<<etaposition<<" etapositionString: "<<etapositionString<<endl;

      c3Eta->cd(i+1);
      if(j==0){
	trees[j]->Draw(treeleg3Eta,theCut&&theDefaultCut&&theEtaCut);
	hmaxDzeta=histosEtadz[0][i]->GetMaximum();    
      }
      else{ trees[j]->Draw(treeleg3Eta,theCut&&theDefaultCut&&theEtaCut,"sames");
	if(histosEtadz[j][i]->GetMaximum() >  hmaxDzeta){
	  hmaxDzeta = histosdz[j][i]->GetMaximum();
	}
      }

      histosEtadz[0][i]->GetYaxis()->SetRangeUser(0.01,hmaxDzeta*1.10);
      histosEtadz[0][i]->Draw("sames");
       
      if(j!=0){
	histosEtadz[j][i]->Draw("sames");
      }

      fitResiduals(histosEtadz[j][i]);   
      tmpsdzEta[j][i] = (TF1*)histosEtadz[j][i]->GetListOfFunctions()->FindObject("tmp");   
      if(tmpsdzEta[j][i])  {
	tmpsdzEta[j][i]->SetLineColor(colors[j]);
	tmpsdzEta[j][i]->SetLineWidth(2); 
      
      
      meansdzEta[j][i]=(tmpsdzEta[j][i])->GetParameter(1);
      sigmasdzEta[j][i]=(tmpsdzEta[j][i])->GetParameter(2);
      
      meansdzEtaError[j][i]=(tmpsdzEta[j][i])->GetParError(1);
      sigmasdzEtaError[j][i]=(tmpsdzEta[j][i])->GetParError(2);

      histoMeansdzEta[j]->SetBinContent(i+1,meansdzEta[j][i]); 
      histoMeansdzEta[j]->SetBinError(i+1,meansdzEtaError[j][i]); 
      histoMeansdzEta[j]->GetXaxis()->SetBinLabel(i+1,etapositionString);
      histoSigmasdzEta[j]->SetBinContent(i+1,sigmasdzEta[j][i]);
      histoSigmasdzEta[j]->SetBinError(i+1,sigmasdzEtaError[j][i]);
      histoSigmasdzEta[j]->GetXaxis()->SetBinLabel(i+1,etapositionString);

      }

      else {
	histoMeansdzEta[j]->SetBinContent(i+1,0.); 
	histoMeansdzEta[j]->SetBinError(i+1,0.); 
	histoMeansdzEta[j]->GetXaxis()->SetBinLabel(i+1,etapositionString);
	histoSigmasdzEta[j]->SetBinContent(i+1,0.);
	histoSigmasdzEta[j]->SetBinError(i+1,0.);
	histoSigmasdzEta[j]->GetXaxis()->SetBinLabel(i+1,etapositionString);
      }


      c3Eta->Draw("sames");   
      c3Eta->cd(i+1);

      statObjsdzEta[j][i] = histosEtadz[j][i]->GetListOfFunctions()->FindObject("stats");
      statsdzEta[j][i] = static_cast<TPaveStats*>(statObjsdzEta[j][i]);
      statsdzEta[j][i]->SetLineColor(colors[j]);
      statsdzEta[j][i]->SetTextColor(colors[j]);
      statsdzEta[j][i]->SetShadowColor(10);
      statsdzEta[j][i]->SetFillColor(10);
      statsdzEta[j][i]->SetX1NDC(0.66);
      statsdzEta[j][i]->SetY1NDC(0.80-j*0.18);
      statsdzEta[j][i]->SetX2NDC(0.99);
      statsdzEta[j][i]->SetY2NDC(0.99-j*0.19);
      statsdzEta[j][i]->Draw("same"); 
      
    }
  }
      
  TString Canvas3Title ="VertexDzResidualsPhiBin_"+label;
  TString Canvas3formatpng=Canvas3Title+".png";
  TString Canvas3formateps=Canvas3Title+".eps";  
  c3->SaveAs(Canvas3formatpng);
  c3->Write();
  //c3->SaveAs(Canvas3formateps);
  

  TString Canvas3TitleEta ="VertexDzResidualsEtaBin_"+label;
  TString Canvas3formatEtapng=Canvas3TitleEta+".png";
  TString Canvas3formatEtaeps=Canvas3TitleEta+".eps";
  c3Eta->SaveAs(Canvas3formatEtapng);
  c3Eta->Write();
  //c3Eta->SaveAs(Canvas3formatEtaeps);

  TPaveText *pt = new TPaveText(0.58,0.70,0.75,0.81,"NDC");
  pt->SetFillColor(10);
  pt->SetTextColor(1);
  //pt->SetTextSize(0.05);
  pt->SetTextFont(42);
  pt->SetTextAlign(11);
  TText *text1 = pt->AddText("CMS preliminary 2011");
  text1->SetTextSize(0.05);
  TText *text2 = pt->AddText("#sqrt{s}=7 TeV");
  text2->SetTextSize(0.05); 
  
  //######################################################
  //  Histograms Means and Widths
  //######################################################
  
  TLegend *legoHisto = new TLegend(0.14,0.71,0.58,0.84); 
  legoHisto->SetLineColor(10);
  
  TLegend *legoHistoSeparation = new TLegend(0.14,0.71,0.58,0.84); 
  legoHistoSeparation->SetLineColor(10);

  for(Int_t j=0; j < nOfFiles; j++) { 

    //###########################
    // Means
    //###########################

    //##################### phi ################## 

    c4->cd(1); 
    if(j==0){
      histoMeansdxy[j]->Draw("e1");
    }
    else  histoMeansdxy[j]->Draw("e1sames");
    
    legoHisto->AddEntry(histoMeansdxy[j],LegLabels[j]);
    legoHisto->SetFillColor(10);
    legoHisto->SetShadowColor(10);
    legoHisto->Draw();
    pt->Draw("sames");

    c5->cd(1);
    if(j==0){
      histoMeansdz[j]->Draw("e1");
    }
    else  histoMeansdz[j]->Draw("e1sames");
    legoHisto->Draw();
    pt->Draw("sames");

    //#########################################
    // 
    // Plots with means only
    //
    //#########################################

    c4monly->cd();

    if(j==0){
      histoMeansdxy[j]->SetMarkerSize(1.1);
      histoMeansdxy[j]->Draw("e1");
    }
    else  histoMeansdxy[j]->Draw("e1sames");
    
    legoHisto->Draw();
    pt->Draw("sames");

    c5monly->cd();
    
    Double_t deltaZ(0);
    Double_t sigmadeltaZ(-1);
    
    if(j==0){
      histoMeansdz[j]->SetMarkerSize(1.1);
      histoMeansdz[j]->Draw("e1");
    }
    else  histoMeansdz[j]->Draw("e1sames");
   
    TCanvas *theNewCanvas2 = new TCanvas("NewCanvas2","Fitting Canvas 2",800,600);
    theNewCanvas2->Divide(2,1);

    TH1F *hnewUp = (TH1F*)histoMeansdz[j]->Clone("hnewUp");
    TH1F *hnewDown = (TH1F*)histoMeansdz[j]->Clone("hnewDown");
    
    fleft[j] = new TF1(Form("fleft_%i",j),fULine,-0.5,2.5,1);
    fright[j] = new TF1(Form("fright_%i",j),fULine,8.5,11.5,1);
    fall[j] = new TF1(Form("fall_%i",j),fDLine,2.5,8.5,1);
    
    FitULine(hnewUp);  
    FitDzUp[j]   = (TF1*)hnewUp->GetListOfFunctions()->FindObject("lineUp"); 
    if(FitDzUp[j]){
      fleft[j]->SetParameters(FitDzUp[j]->GetParameters());
      fleft[j]->SetParErrors(FitDzUp[j]->GetParErrors());
      hnewUp->GetListOfFunctions()->Add(fleft[j]);
      fright[j]->SetParameters(FitDzUp[j]->GetParameters());
      fright[j]->SetParErrors(FitDzUp[j]->GetParErrors());
      hnewUp->GetListOfFunctions()->Add(fright[j]);
      FitDzUp[j]->Delete();

      theNewCanvas2->cd(1);
      makeNiceTF1Style(fright[j],colors[j]);
      makeNiceTF1Style(fleft[j],colors[j]);
      fright[j]->Draw("same");
      fleft[j]->Draw("same");
    }
    
    FitDLine(hnewDown);  
    FitDzDown[j] = (TF1*)hnewDown->GetListOfFunctions()->FindObject("lineDown");    
    
    if(FitDzDown[j]){
      fall[j]->SetParameters(FitDzDown[j]->GetParameters());
      fall[j]->SetParErrors(FitDzDown[j]->GetParErrors());
      hnewDown->GetListOfFunctions()->Add(fall[j]);
      FitDzDown[j]->Delete();
      
      theNewCanvas2->cd(2);
      makeNiceTF1Style(fall[j],colors[j]);
      fall[j]->Draw("same");
      c5monly->cd();
      fright[j]->Draw("sames");
      fleft[j]->Draw("same");
      fall[j]->Draw("same");
    }
    
    if(j==nOfFiles-1){
      theNewCanvas2->Close();
    }

    deltaZ=(fright[j]->GetParameter(0) - fall[j]->GetParameter(0));
    sigmadeltaZ=TMath::Sqrt(fright[j]->GetParError(0)*fright[j]->GetParError(0) + fall[j]->GetParError(0)*fall[j]->GetParError(0));
    TString COUT = Form(" #Delta z = %.f #pm %.f #mum",deltaZ,sigmadeltaZ);
    
    if(j==nOfFiles-1){ 
      outfile <<deltaZ<<"|"<<sigmadeltaZ<<endl;
    }

    legoHistoSeparation->AddEntry(histoMeansdxy[j],LegLabels[j]+COUT);
    legoHistoSeparation->SetFillColor(10);
    legoHistoSeparation->SetShadowColor(10);
    legoHistoSeparation->Draw();

    pt->Draw("sames");

    //##################### Eta ################## 

    c4Eta->cd(1); 
    if(j==0){
      histoMeansdxyEta[j]->Draw("e1");
    }
    else  histoMeansdxyEta[j]->Draw("e1sames");
    
    legoHisto->Draw();
    pt->Draw("sames");

    c5Eta->cd(1);
    if(j==0){
      histoMeansdzEta[j]->Draw("e1");
    }
    else  histoMeansdzEta[j]->Draw("e1sames");
    legoHisto->Draw();
    pt->Draw("sames");

    //###########################
    // Sigmas
    //###########################
    
    //##################### Phi ################## 

    c4->cd(2);
    if(j==0){
      histoSigmasdxy[j]->Draw("e1");
    }
    else histoSigmasdxy[j]->Draw("e1sames");
    legoHisto->Draw();
    pt->Draw("sames");

    c5->cd(2);
    if(j==0){
      histoSigmasdz[j]->Draw("e1");
    }
    else histoSigmasdz[j]->Draw("e1sames");
    legoHisto->Draw();
    pt->Draw("sames");

    //##################### Eta ##################  

    c4Eta->cd(2);
    if(j==0){
      histoSigmasdxyEta[j]->Draw("e1");
    }
    else histoSigmasdxyEta[j]->Draw("e1sames");
    legoHisto->Draw();
    pt->Draw("sames");

    c5Eta->cd(2);
    if(j==0){
      histoSigmasdzEta[j]->Draw("e1");
    }
    else histoSigmasdzEta[j]->Draw("e1sames");
    legoHisto->Draw();
    pt->Draw("sames");
      
  }
  
  TString Canvas4Title ="HistoMeansSigmasDxy_"+label;
  TString Canvas4formatpng=Canvas4Title+".png"; 
  TString Canvas4formateps=Canvas4Title+".eps"; 
  c4->SaveAs(Canvas4formatpng);
  c4->Write();
  // c4->SaveAs(Canvas4formateps);
  
  TString Canvas4Titlemonly ="HistoMeansDxy_"+label;
  TString Canvas4formatpngmonly=Canvas4Titlemonly+".png"; 
  TString Canvas4formatepsmonly=Canvas4Titlemonly+".eps"; 
  TString Canvas4formatrootmonly=Canvas4Titlemonly+".root"; 
  c4monly->SaveAs(Canvas4formatpngmonly);
  c4monly->Write();
  //c4monly->SaveAs(Canvas4formatepsmonly);
  //c4monly->SaveAs(Canvas4formatrootmonly);

  TString Canvas5Title ="HistoMeansSigmasDz_"+label;
  TString Canvas5formatpng=Canvas5Title+".png"; 
  TString Canvas5formateps=Canvas5Title+".eps";
  c5->SaveAs(Canvas5formatpng);
  c5->Write();
  //  c5->SaveAs(Canvas5formateps);

  TString Canvas5Titlemonly ="HistoMeansDzFit_"+label;
  TString Canvas5formatpngmonly=Canvas5Titlemonly+".png"; 
  TString Canvas5formatepsmonly=Canvas5Titlemonly+".eps"; 
  TString Canvas5formatrootmonly=Canvas5Titlemonly+".root"; 
  c5monly->SaveAs(Canvas5formatpngmonly);
  c5monly->Write();
  //c5monly->SaveAs(Canvas5formatepsmonly);
  //c5monly->SaveAs(Canvas5formatrootmonly);

  TString Canvas4TitleEta ="HistoMeansSigmasDxyEta_"+label;
  TString Canvas4formatEtapng=Canvas4TitleEta+".png"; 
  TString Canvas4formatEtaeps=Canvas4TitleEta+".eps"; 
  c4Eta->SaveAs(Canvas4formatEtapng);
  c4Eta->Write();
  //c4Eta->SaveAs(Canvas4formatEtaeps);
  
  TString Canvas5TitleEta ="HistoMeansSigmasDzEta_"+label;
  TString Canvas5formatEtapng=Canvas5TitleEta+".png";
  TString Canvas5formatEtaeps=Canvas5TitleEta+".eps";
  c5Eta->SaveAs(Canvas5formatEtapng);
  c5Eta->Write();
  //c5Eta->SaveAs(Canvas5formatEtaeps);

  outfile.close();

}

//##########################################
// Fitting Function
//##########################################

void fitResiduals(TH1 *hist)
{
  //float fitResult(9999);
  if (!hist || hist->GetEntries() < 20) return;
  
  float mean  = hist->GetMean();
  float sigma = hist->GetRMS();
  
  TF1 func("tmp", "gaus", mean - 1.5*sigma, mean + 1.5*sigma); 
  if (0 == hist->Fit(&func,"QNR")) { // N: do not blow up file by storing fit!
    mean  = func.GetParameter(1);
    sigma = func.GetParameter(2);
    // second fit: three sigma of first fit around mean of first fit
    func.SetRange(mean - 2*sigma, mean + 2*sigma);
      // I: integral gives more correct results if binning is too wide
      // L: Likelihood can treat empty bins correctly (if hist not weighted...)
    if (0 == hist->Fit(&func, "Q0LR")) {
      if (hist->GetFunction(func.GetName())) { // Take care that it is later on drawn:
	hist->GetFunction(func.GetName())->ResetBit(TF1::kNotDraw);
      }
    }
  }
  return;
}


void DoubleGausFitResiduals(TH1 *hist)
{
  //float fitResult(9999);
  // if (!hist || hist->GetEntries() < 20) return;
  
  float mean  = hist->GetMean();
  float sigma = hist->GetRMS();

  Double_t par[6];
  TF1 func1("gaus1", "gaus", mean - 1.5*sigma, mean + 1.5*sigma); 
  TF1 func2("gaus2", "gaus", mean - 2.5*sigma, mean + 2.5*sigma); 
  TF1 func("tmp", "gaus(0)+gaus(3)", mean - 2.5*sigma, mean + 2.5*sigma); 
  
  hist->Fit(&func1,"QNR");
  hist->Fit(&func2,"QNR+");
  func1.GetParameters(&par[0]);
  func2.GetParameters(&par[3]);
  cout<<"partials fit done!"<<endl;
  
  if(hist->GetEntries()>20) {
    
    cout<<"histo entries: "<<hist->GetEntries()<<endl; 

    func.SetParameters(par);
    if (0 == hist->Fit(&func,"QNR")) { // N: do not blow up file by storing fit!
      
      // cout<<"before extracting parameters"<<endl;
      
      mean  = func.GetParameter(1);
      sigma = func.GetParameter(5);
      
      cout<<"first total fit done!"<<endl;
      
      // second fit: three sigma of first fit around mean of first fit
      func.SetRange(mean - 2.5*sigma, mean + 2.5*sigma);
      // I: integral gives more correct results if binning is too wide
      // L: Likelihood can treat empty bins correctly (if hist not weighted...)
      if (0 == hist->Fit(&func, "Q0LR")) {
	if (hist->GetFunction(func.GetName())) { // Take care that it is later on drawn:
	  hist->GetFunction(func.GetName())->ResetBit(TF1::kNotDraw);
	}
 
	//cout<<"fit done!"<<endl;
	
      }
    }
    
  } else {
    cout<<"Unable to perform double gaussian fit "<<endl;
    // func.SetParameters(0);
    hist->Fit(&func, "Q0LR");
  }
  
  return;
}

Double_t fDLine(Double_t *x, Double_t *par)
{
  if (x[0] < 2.5 && x[0] > 8.5) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0];
}

Double_t fULine(Double_t *x, Double_t *par)
{
  if (x[0] >= 2.5 && x[0] <= 8.5) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0];
}

void FitULine(TH1 *hist)
{ 
  // define fitting function
  TF1 func1("lineUp",fULine,-0.5,11.5,1);
  //TF1 func1("lineUp","pol0",-0.5,11.5);
  
  if (0 == hist->Fit(&func1,"QR")) {
    if (hist->GetFunction(func1.GetName())) { // Take care that it is later on drawn:
      hist->GetFunction(func1.GetName())->ResetBit(TF1::kNotDraw);
    }
    cout<<"fit Up done!"<<endl;
  }
  
}

void FitDLine(TH1 *hist)
{
  // define fitting function
  // TF1 func1("lineDown",fDLine,-0.5,11.5,1);
  
  TF1 func2("lineDown","pol0",2.5,8.5);
  func2.SetRange(2.5,8.5);
  
  if (0 == hist->Fit(&func2,"QR")) {
    if (hist->GetFunction(func2.GetName())) { // Take care that it is later on drawn:
      hist->GetFunction(func2.GetName())->ResetBit(TF1::kNotDraw);
    }
    cout<<"fit Down done!"<<endl;
  } 
}

void makeNiceTF1Style(TF1 *f1,int color){
  f1->SetLineColor(color);
  f1->SetLineWidth(3);
  f1->SetLineStyle(2);
}


void MakeNiceHistoStyle(TH1 *hist){
  
  //hist->SetTitleSize(0.09); 
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetXaxis()->SetTitleOffset(0.9);
  hist->GetYaxis()->SetTitleSize(0.08);
  hist->GetXaxis()->SetTitleSize(0.07);
  hist->GetXaxis()->SetTitleFont(42); 
  hist->GetYaxis()->SetTitleFont(42);  
  hist->GetXaxis()->SetLabelSize(0.05); 
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->SetMarkerSize(1.2);
  
}

void MakeNiceTrendPlotStyle(TH1 *hist){
  
  hist->SetStats(kFALSE);  
  hist->SetLineWidth(2);
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  hist->GetXaxis()->SetTitleFont(42); 
  hist->GetYaxis()->SetTitleFont(42);  
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetTitleOffset(1.0);
  hist->GetYaxis()->SetTitleOffset(1.0);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(.06);
  hist->GetXaxis()->SetLabelSize(.06);
  hist->SetMarkerSize(1);
}
