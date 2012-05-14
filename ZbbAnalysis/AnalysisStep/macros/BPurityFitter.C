#include <TFile.h>
#include <TTree.h>
#include <iostream>
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
#include <TVectorF.h>
#include <TKey.h>
#include <TMath.h>
#include <TFractionFitter.h>
#include <Math/VectorUtil.h>
#include <fstream>

ofstream textfile_;
Float_t ch[3],ech[3],f_b[3],f_c[3],f_l[3],ef_b[3],ef_c[3],ef_l[3];
Int_t theStackColors[3]={628,409,856};
TString binlabels[3]={"ll","#mu#mu","ee"};

// per jet purity calculator
std::pair<Float_t,Float_t> BFractionFitter(TString tagger_,TString channel_,TString variable_, TFile* fileout);

// per event purity calculator
std::pair<Float_t,Float_t> EventPurityCalculator(TH1F *nbjetshisto,Float_t perJetPurity,Float_t perJetPurityerr);
Float_t purFunction(Float_t P,Int_t njets);

// graphics
void CMSPrel(Double_t Lumi);
void setStyle();
void MakeNiceHistoStyle(TH1F *hist, Int_t color=-999,Bool_t fillInside=true);
void MakeNiceGraphStyle(TGraph *hist, Int_t color);
TLegend* MakeSimpleTLegend(TObjArray *array,TString *LegLabels);
Double_t getMaximum(TObjArray *array);
TH1F* takeTheRatio(TH1F* h1,TH1F* h2);
void makeWaterMark(Double_t xNDC,Double_t yNDC, TString mytext);

// main method
void BPurityFitter(TString tagger_="ssvhpt",TString variable_="SVmass"){

  TString fileoutname = "fits_"+tagger_+variable_+".root";
  TString textoutname = "fits_"+tagger_+variable_+".txt";
  
  TFile* fileout = new TFile(fileoutname,"RECREATE");

  textfile_.open(textoutname);
  textfile_.precision(2);

  std::pair<Float_t,Float_t> llJetPur = BFractionFitter(tagger_,"",variable_,fileout);
  std::pair<Float_t,Float_t> mmJetPur = BFractionFitter(tagger_,"_mm",variable_,fileout);
  std::pair<Float_t,Float_t> eeJetPur = BFractionFitter(tagger_,"_ee",variable_,fileout);

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

  gr_b->GetYaxis()->SetRangeUser(0.,1.8);
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
  
  std::pair<Float_t,Float_t> llEventPur = EventPurityCalculator(nbjets_all,llJetPur.first,llJetPur.second);
  std::pair<Float_t,Float_t> mmEventPur = EventPurityCalculator(nbjets_mm ,mmJetPur.first,mmJetPur.second);
  std::pair<Float_t,Float_t> eeEventPur = EventPurityCalculator(nbjets_ee ,eeJetPur.first,eeJetPur.second);

  std::cout<<"=============== final report =============="<<std::endl;
  std::cout<<"Jet-purity ll = "<<llJetPur.first<<"+/-"<<llJetPur.second<<" Event-purity = "<< llEventPur.first <<"+/-"<<llEventPur.second<<std::endl;
  std::cout<<"Jet-purity mm = "<<mmJetPur.first<<"+/-"<<mmJetPur.second<<" Event-purity = "<< mmEventPur.first <<"+/-"<<mmEventPur.second<<std::endl;
  std::cout<<"Jet-purity ee = "<<eeJetPur.first<<"+/-"<<eeJetPur.second<<" Event-purity = "<< eeEventPur.first <<"+/-"<<eeEventPur.second<<std::endl;
  std::cout<<"==========================================="<<std::endl;

  textfile_<<"=============== final report =============="<<std::endl;
  textfile_<<"Jet-purity ll = "<<llJetPur.first<<"+/-"<<llJetPur.second<<" Event-purity = "<< llEventPur.first <<"+/-"<<llEventPur.second<<std::endl;
  textfile_<<"Jet-purity mm = "<<mmJetPur.first<<"+/-"<<mmJetPur.second<<" Event-purity = "<< mmEventPur.first <<"+/-"<<mmEventPur.second<<std::endl;
  textfile_<<"Jet-purity ee = "<<eeJetPur.first<<"+/-"<<eeJetPur.second<<" Event-purity = "<< eeEventPur.first <<"+/-"<<eeEventPur.second<<std::endl;
  textfile_<<"==========================================="<<std::endl;

  textfile_.close();

}

// calculation method
std::pair<Float_t,Float_t> BFractionFitter(TString tagger_,TString channel_,TString variable_,TFile* fileout){

  TString n_ch="b+Z#rightarrow ll";
  if (channel_=="_mm"){
    n_ch = "b+Z#rightarrow #mu#mu";
  } else if (channel_=="_ee"){
    n_ch = "b+Z#rightarrow ee";
  }
    
  setStyle();
  
  TCanvas *c1 =new TCanvas("c1"+channel_,"c1",1000,500);
  c1->Divide(2,1);
  
  TCanvas *c2 =new TCanvas("c2"+channel_,"c2",1200,600);
  TPad *c2pad1a = new TPad("c2pad1a","The pad",0.  ,0.24,0.33,1.);
  TPad *c2pad2a = new TPad("c2pad2a","The pad",0.33,0.24,0.66,1.);
  TPad *c2pad3a = new TPad("c2pad3a","The pad",0.66,0.24,0.99,1.);
  TPad *c2pad1b = new TPad("c2pad1b","The pad",0.  ,0.01 ,0.33,0.3);
  TPad *c2pad2b = new TPad("c2pad2b","The pad",0.33,0.01 ,0.66,0.3);
  TPad *c2pad3b = new TPad("c2pad3b","The pad",0.66,0.01 ,0.99,0.3);

  c2pad1a->Draw();
  c2pad2a->Draw();
  c2pad3a->Draw();

  c2pad1b->SetBottomMargin(0.35);
  c2pad2b->SetBottomMargin(0.35);
  c2pad3b->SetBottomMargin(0.35);

  c2pad1b->Draw();
  c2pad2b->Draw();
  c2pad3b->Draw();

  // c2->Divide(3,1);
 
  TCanvas *c3 =new TCanvas("c3"+channel_,"c3",800,600);
  TPad *c3pada = new TPad("c3pada","The pad",0.,0.24,1.,1.);
  TPad *c3padb = new TPad("c3padb","The pad",0.,0.01,1.,0.3);  
  c3pada->Draw();
  c3padb->Draw();
  c3padb->SetBottomMargin(0.35);

  TString filename =" templates_"+tagger_+".root";
  TFile* file = TFile::Open(filename);

  TH1F *data = (TH1F*)file->Get(variable_+"data"+channel_);  // data histogram
  TH1F *mcB  = (TH1F*)file->Get(variable_+"_mc_b"+channel_); // first MC histogram
  TH1F *mcC  = (TH1F*)file->Get(variable_+"_mc_c"+channel_); // second MC histogram
  TH1F *mcL  = (TH1F*)file->Get(variable_+"_mc_l"+channel_); // third MC histogram

  TObjArray *mc = new TObjArray(3);                // MC histograms are put in this array
  mc->Add(mcB);
  mc->Add(mcC);
  mc->Add(mcL);

  // first part of the canvas (MC without fitting)
  c1->cd(1);
  MakeNiceHistoStyle(data);
  MakeNiceHistoStyle(mcB,theStackColors[0],false);
  MakeNiceHistoStyle(mcC,theStackColors[1],false);
  MakeNiceHistoStyle(mcL,theStackColors[2],false);

  TObjArray *all = new TObjArray(4);                // all histograms are put in this array
  all->Add(data);
  all->Add(mcB);
  all->Add(mcC);
  all->Add(mcL);
  
  THStack *hs=new THStack("hs","test stacked histograms");
  hs->Add(mcB);
  hs->Add(mcC);
  hs->Add(mcL);

  TObjArray *dmc = new TObjArray(2);
  dmc->Add(data);
  dmc->Add(hs->GetStack()->Last());
  Double_t themax1 = getMaximum(dmc);

  data->SetMaximum(themax1*1.10);
  data->Draw("Ep");
  hs->Draw("HISTsame");
  data->Draw("Epsame");
 
  TString LegLabels[4]={"CMS DATA 2011","b jets","c jets","light jets"};

  TLegend* leg = MakeSimpleTLegend(all,LegLabels);
  
  leg->Draw("same");
  CMSPrel(2.2);
  
  makeWaterMark(0.2,0.93,"out-of-the-box MC");
  makeWaterMark(0.2,0.89,n_ch);

  // fit computation
  Double_t fmc1=0., fmc1err=0.;
  Double_t fmc2=0., fmc2err=0.;
  Double_t fmc3=0., fmc3err=0.;

  TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
  //fit->Constrain(1,0.0,1.0);             // constrain fraction 1 to be between 0 and 1
  Int_t status = fit->Fit();               // perform the fit
  cout << "fit status: " << status << endl;
  
  if (status == 0) {                       // check on fit status
    
    fit->GetResult(0,fmc1,fmc1err);
    fit->GetResult(1,fmc2,fmc2err);
    fit->GetResult(2,fmc3,fmc3err);
    
    TString output = Form("#chi^{2}=%5.3f ndf=%3d [P=%5.3f]",fit->GetChisquare(),fit->GetNDF(),fit->GetProb());
    
    // canvas of results
    c1->cd(2);
    data->Draw("Ep");
    TH1F* result = (TH1F*) fit->GetPlot();
    result->SetLineWidth(2.);
    
    TH1F* resultB=(TH1F*) fit->GetMCPrediction(0);
    TH1F *resultB_new=(TH1F*)resultB->Clone();
    resultB_new->Scale(fmc1*data->Integral()/resultB->GetSumOfWeights());
    MakeNiceHistoStyle(resultB,theStackColors[0]);
    MakeNiceHistoStyle(resultB_new,theStackColors[0]);
    
    TH1F* resultC=(TH1F*) fit->GetMCPrediction(1);
    TH1F *resultC_new=(TH1F*)resultC->Clone();
    resultC_new->Scale(fmc2*data->Integral()/resultC->GetSumOfWeights());
    MakeNiceHistoStyle(resultC,theStackColors[1]);
    MakeNiceHistoStyle(resultC_new,theStackColors[1]);

    TH1F* resultL=(TH1F*) fit->GetMCPrediction(2);
    TH1F *resultL_new=(TH1F*)resultL->Clone();
    resultL_new->Scale(fmc3*data->Integral()/resultL->GetSumOfWeights());
    MakeNiceHistoStyle(resultL,theStackColors[2]);
    MakeNiceHistoStyle(resultL_new,theStackColors[2]);

    TObjArray *fractions = new TObjArray(4);   
    fractions->Add(data);
    fractions->Add(resultB_new);
    fractions->Add(resultC_new);
    fractions->Add(resultL_new);

    THStack *hs2=new THStack("hs2","test stacked histograms");
    hs2->Add(resultB_new);
    hs2->Add(resultC_new);
    hs2->Add(resultL_new);
    
    hs2->Draw("HISTsame");
    data->Draw("Epsame");
    CMSPrel(2.2);

    TString LegLabels2[4]={"CMS Data",Form("f_{b}: %5.3f #pm %5.3f",fmc1,fmc1err),Form("f_{c}: %5.3f #pm %5.3f",fmc2,fmc2err),Form("f_{l}: %5.3f #pm %5.3f",fmc3,fmc3err)};
    TLegend* leg2 = MakeSimpleTLegend(fractions,LegLabels2);
    leg2->Draw("same");

    makeWaterMark(0.2,0.93,"after fit");
    makeWaterMark(0.2,0.89,n_ch);

    // canvas of comparison and result 
    c3->cd();
    c3pada->cd();
    data->Draw("Ep");
    result->Draw("same");
    CMSPrel(2.2);
    makeWaterMark(0.17,0.91,n_ch);
    
    TLatex *ll = new TLatex();
    ll->SetTextSize(0.04);
    ll->SetTextFont(42);
    ll->SetLineWidth(2);
    ll->SetNDC();
    ll->DrawLatex(0.6,0.68,Form("f_{b} = %5.3f #pm %5.3f",fmc1,fmc1err));
    ll->DrawLatex(0.6,0.62,Form("f_{c} = %5.3f #pm %5.3f",fmc2,fmc2err));
    ll->DrawLatex(0.6,0.56,Form("f_{l} = %5.3f #pm %5.3f",fmc3,fmc3err));
    ll->DrawLatex(0.6,0.50,output);
    
    textfile_<<"=================== "<<n_ch<<" ================="<<endl;
    textfile_<<Form("f_{b} = %5.3f +/- %5.3f",fmc1,fmc1err)<<endl;
    textfile_<<Form("f_{c} = %5.3f +/- %5.3f",fmc2,fmc2err)<<endl;
    textfile_<<Form("f_{l} = %5.3f +/- %5.3f",fmc3,fmc3err)<<endl;
    textfile_<<endl;
    textfile_<<output<<endl;

    if(channel_=="") {
      f_b[0]=fmc1;
      f_c[0]=fmc2; 
      f_l[0]=fmc3;
      ef_b[0]=fmc1err;
      ef_c[0]=fmc2err; 
      ef_l[0]=fmc3err;
    } else if (channel_=="_mm"){
      f_b[1]=fmc1;
      f_c[1]=fmc2; 
      f_l[1]=fmc3;
      ef_b[1]=fmc1err;
      ef_c[1]=fmc2err; 
      ef_l[1]=fmc3err;
    } else {
      f_b[2]=fmc1;
      f_c[2]=fmc2; 
      f_l[2]=fmc3;
      ef_b[2]=fmc1err;
      ef_c[2]=fmc2err; 
      ef_l[2]=fmc3err;
    }
    
    c3padb->cd();
    TH1F* hratioDATAMC = takeTheRatio(data,result);
    hratioDATAMC->Draw("PE1");

    // canvas of template comparisons
    //------------------------------------
    c2pad1a->cd();
    //c2->cd(1);
    TH1F *mcB_new=(TH1F*)mcB->Clone();
    mcB_new->SetMarkerStyle(20);
    mcB_new->SetMarkerSize(1.);
    mcB_new->SetLineColor(1);

    TObjArray *compB = new TObjArray(2); 
    compB->Add(mcB_new);
    compB->Add(resultB);

    Double_t themaxB = getMaximum(compB);
    
    resultB->SetMaximum(themaxB*1.20);
    resultB->Draw("HIST");
    mcB_new->Draw("Psame");

    TString LegLabels3[2]={"b- template","b- after fit"};
    TLegend* leg3 = MakeSimpleTLegend(compB,LegLabels3);
    leg3->Draw();
    CMSPrel(2.2);
    makeWaterMark(0.195,0.91,n_ch);

    c2pad1b->cd();
    TH1F* hratioB = takeTheRatio(mcB_new,resultB);
    hratioB->Draw("PE1");

    //------------------------------------
    c2pad2a->cd();
    // c2->cd(2);
    TH1F *mcC_new=(TH1F*)mcC->Clone();
    mcC_new->SetMarkerStyle(20);
    mcC_new->SetMarkerSize(1.);
    mcC_new->SetLineColor(1);

    TObjArray *compC = new TObjArray(2); 
    compC->Add(mcC_new);
    compC->Add(resultC);

    Double_t themaxC = getMaximum(compC);
    
    resultC->SetMaximum(themaxC*1.20);
    resultC->Draw("HIST");
    mcC_new->Draw("Psame");

    TString LegLabels4[2]={"c- template","c- after fit"};
    TLegend* leg4 = MakeSimpleTLegend(compC,LegLabels4);
    leg4->Draw();
    CMSPrel(2.2);
    makeWaterMark(0.195,0.91,n_ch);

    c2pad2b->cd();
    TH1F* hratioC = takeTheRatio(mcC_new,resultC);
    hratioC->Draw("PE1");

    //------------------------------------
    c2pad3a->cd();
    //c2->cd(3);
    TH1F *mcL_new=(TH1F*)mcL->Clone();
    mcL_new->SetMarkerStyle(20);
    mcL_new->SetMarkerSize(1.);
    mcL_new->SetLineColor(1);

    TObjArray *compL = new TObjArray(2); 
    compL->Add(mcL_new);
    compL->Add(resultL);
    
    Double_t themaxL = getMaximum(compL);
    
    resultL->SetMaximum(themaxL*1.20);
    resultL->Draw("HIST");
    mcL_new->Draw("Psame");

    TString LegLabels5[2]={"l- template","l- after fit"};
    TLegend* leg5 = MakeSimpleTLegend(compL,LegLabels5);
    leg5->Draw();
    CMSPrel(2.2);
    makeWaterMark(0.195,0.91,n_ch);
    
    c2pad3b->cd();
    TH1F* hratioL = takeTheRatio(mcL_new,resultL);
    hratioL->Draw("PE1");

  }

  c1->SaveAs("purity"+tagger_+variable_+channel_+".png");
  c2->SaveAs("templatecontrol"+tagger_+variable_+channel_+".png");
  c3->SaveAs("finalfit"+tagger_+variable_+channel_+".png");
  
  fileout->cd();
  c1->Write();
  c2->Write();
  c3->Write();

  std::pair<Float_t,Float_t> result;
  result = make_pair(fmc1,fmc1err);
  
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
    }
  }  else {
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.0);
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

  latexLabel->DrawLatex(0.55,0.9,"CMS VHF Preliminary");
  latexLabel->DrawLatex(0.55,0.84,"anti-k_{T} (R = 0.5) PF Jets ");
  latexLabel->DrawLatex(0.55,0.78,(TString)Form("#sqrt{s} = 7 TeV  #int Ldt = %.1f fb^{-1}",Lumi));  

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
std::pair<Float_t,Float_t> EventPurityCalculator(TH1F *nbjetshisto,Float_t perJetPurity,Float_t perJetPurityerr){

  cout.precision(3);

  Float_t P_b = perJetPurity;
  std::vector<std::pair<Int_t,Float_t> > mult_and_frac;
  
  for(Int_t n=0; n<nbjetshisto->GetNbinsX(); n++){
    if(nbjetshisto->GetBinContent(n)>0.){ 
      std::pair<Int_t,Float_t> bin_and_frac;
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
    std::cout<<"bin= "<<fixed<< mult_and_frac[i].first<<" f("<<mult_and_frac[i].first<<"b)="<<fixed<<mult_and_frac[i].second<<
      " P("<<mult_and_frac[i].first<<"b)="<<fixed<<purFunction(P_b,mult_and_frac[i].first)<<" Pur="<<fixed<<perEventPurity<<endl;
  }

  std::cout<<"........."<<std::endl;

  perEventPurityError=TMath::Sqrt(TMath::Power(diff,2)*TMath::Power(perJetPurityerr,2));

  std::pair<Float_t,Float_t> result;
  result = make_pair(perEventPurity,perEventPurityError);

  return result;

}

////////////////////////////////////////////////////////////////////
Float_t purFunction(Float_t P,Int_t njets){
  
  Float_t P_e = 1.-TMath::Power((1.-P),njets);
  return P_e;
}
