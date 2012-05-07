#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TLatex.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TArrow.h>
#include <iomanip>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <TLine.h>
#include <TPaveText.h>
#include <TFractionFitter.h>

using namespace std;
using namespace TMath;

////////////////////////////////////////
void setStyle();

////////////////////////////////////////
void makeNiceHistoStyle(TH1F *h, Int_t color=0);

////////////////////////////////////////
void CMSPrel(Float_t Lumi, TString _decaychannel);

////////////////////////////////////////
TH1F* takeTheRatio(TH1F* h1,TH1F* h2);

////////////////////////////////////////
Float_t getWeightedBinomialError(TFile *file, Float_t myCounts, Float_t mySquaredCounts);

////////////////////////////////////////
Float_t GetTheIntegral(TH1F *histo,Int_t _firstbin, Int_t _lastbin,Bool_t getTheSquare=false);

////////////////////////////////////////
Double_t getMaximum(TObjArray *array);

////////////////////////////////////////
Float_t CalculateError(pair<Float_t,Float_t> R_tt_pair, 
		       pair<Float_t,Float_t> R_Z_pair, 
		       pair<Float_t,Float_t> N_obs_in_pair, 
		       pair<Float_t,Float_t> N_obs_out_pair);

////////////////////////////////////////
Float_t propErrRatio(Float_t num,Float_t err_num,Float_t den,Float_t err_den){
  
  return (num/den)*Sqrt( Power((err_num/num),2) + Power((err_den/den),2));
};


////////////////////////////////////////
// main method
////////////////////////////////////////
void estimatefTTBar(TString foldername = "JetBTag",TString theALGOWP = "ssvhem",
		    Int_t rebin = 1,Bool_t LogScale = false,Bool_t isMuonChannel = true,Float_t maxRegion=200.,Bool_t useMyAlgo=true){

  setStyle();
  
  ofstream outfile_;
  TString  outfilename_="report_fttbar.txt";

  outfile_.open(outfilename_,ios::out | ios::app);

  TString histoname;
  TString _decaychannel;

  if(isMuonChannel){ 
    histoname ="MmumuNOMASSCUT";
    _decaychannel="Z #rightarrow #mu#mu";
  } else{ 
    histoname ="MeeNOMASSCUT";
    _decaychannel="Z-#rightarrow ee";
  }

  //***************
  //Get files
  //***************
  Int_t nOfFiles = 5;
  TFile *files[nOfFiles];
  TString filenames[nOfFiles];
  filenames[0] = "analyzePAT_DATA2011.root";
  filenames[1] = "analyzePAT_MC_TTJets_All.root";
  filenames[2] = "analyzePAT_MC_ZlJets_All.root";
  filenames[3] = "analyzePAT_MC_ZcToLL_All.root";
  filenames[4] = "analyzePAT_MC_Zb5fToLL_All.root";

  for (Int_t i = 0;i < nOfFiles;i++) files[i] = TFile::Open(filenames[i]);
  
  //****************
  //Compute weights
  //****************
  Float_t lumi = 2130.;
  Float_t weights[nOfFiles],xsections[nOfFiles],nevts[nOfFiles];
  nevts[0] = 1;
  for (Int_t i = 1; i < nOfFiles; i++) nevts[i] = ((TH1F*)files[i]->Get("analyzePat/Selevents/SelectedEvts"))->GetBinContent(1);
  
  xsections[0] = 1.;
  xsections[1] = 165.;
  xsections[4] = 3048.;
  xsections[3] = 3048.;
  xsections[2] = 3048.;

  weights[0] = 1;
  for (Int_t i = 1; i < nOfFiles; i++) weights[i] = lumi*xsections[i]/nevts[i];
  
  //*****************
  //Get histograms
  //*****************
  //Get histograms from files
  TH1F *histos[nOfFiles];   
  for (Int_t i = 0; i < nOfFiles; i++){
    histos[i] = (TH1F*)files[i]->Get("finaldistros_"+theALGOWP+"/"+foldername+"/"+foldername+"_"+histoname);
    //  histos[i] = (TH1F*)files[i]->Get("steps_before_btag/"+foldername+"/"+foldername+"_"+histoname);
    histos[i]->Scale(weights[i]);
  }

  //Eventually rebin histograms
  for (Int_t i = 0; i < nOfFiles; i++) histos[i]->Rebin(rebin);

  TObjArray *arrayhistos = new TObjArray(5);    
  for (Int_t i = 0; i < nOfFiles; i++){
    arrayhistos->Add(histos[i]);
  }//loop over histograms
  
  Float_t maxYvalue = getMaximum(arrayhistos);
 
  cout<<"Histogram maxBinContent = "<<maxYvalue<<endl;

  //****************
  //Histogram style
  //****************
  histos[0]->SetMarkerStyle(20);
  histos[0]->SetMarkerSize(1.0);
  // histos[1]->SetFillColor(393);  
  // histos[2]->SetFillColor(856);
  // histos[3]->SetFillColor(409);
  // histos[4]->SetFillColor(628);

  Int_t theStackColors[5] ={0,393,856,409,628};

  for (Int_t i = 0; i < nOfFiles; i++){ 
    makeNiceHistoStyle(histos[i],theStackColors[i]);
    histos[i]->GetYaxis()->SetRangeUser(0.5,maxYvalue + 1.5*maxYvalue);
    histos[i]->GetYaxis()->SetTitle("N_{evts}"); 
  }

  //***************
  //TLegend
  //***************
  TLegend* leg = new TLegend(0.68,0.54,0.85,0.75);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(10);
  leg->SetTextFont(42);
  leg->AddEntry(histos[0],"CMS DATA","P");
  leg->AddEntry(histos[1],"MC t#bar{t}","f");
  leg->AddEntry(histos[2],"MC Z+l","f");
  leg->AddEntry(histos[3],"MC Z+c","f");
  leg->AddEntry(histos[4],"MC Z+b","f");

  //***************
  //THStack
  //***************

  THStack *hs = new THStack("hs",foldername+"_"+theALGOWP+"_"+histoname);
  for (Int_t i = 1; i < nOfFiles; i++) hs->Add(histos[i]);
  
  TCanvas *canv1 = new TCanvas("canv1","canv1",800,600);
  canv1->SetFillColor(0);

  canv1->cd();
  histos[0]->Draw("Pe");
  hs->Draw("HISTSAME");
  histos[0]->Draw("esame");
  leg->Draw("same");
  
  CMSPrel(2.2,_decaychannel);

  if (LogScale)canv1->cd()->SetLogy();
  canv1->Draw();
  canv1->Update();

  // ------> Second Canvas <------

  TH1F *h_tt  = (TH1F*)histos[1]->Clone(); 
  TH1F *h_Z   = (TH1F*)histos[2]->Clone(); 
 
  h_Z->Add(histos[3]);
  h_Z->Add(histos[4]);

  THStack *hs2 = new THStack("hs2",foldername+"_"+theALGOWP+"_"+histoname);
  hs2->Add(h_tt);
  hs2->Add(h_Z);

  TCanvas *canv2 = new TCanvas("canv2","canv2",800,600);
  canv2->SetFillColor(0);

  TLegend* leg2 = new TLegend(0.48,0.72,0.65,0.93);
  leg2->SetFillColor(0);
  leg2->SetLineColor(0);
  leg2->SetShadowColor(10);
  leg2->SetTextFont(42);
  leg2->AddEntry(histos[0],"CMS DATA","P");
  leg2->AddEntry(h_tt,"MC t#bar{t}","f");
  leg2->AddEntry(h_Z,"MC Z+jets","f");

  canv2->cd();
  histos[0]->Draw("Pe");
  hs2->Draw("HISTSAME");
  histos[0]->Draw("esame");

  CMSPrel(2.2,_decaychannel);
  leg2->Draw("same");

  // TArrow *l1=new TArrow(c->GetUxmin(),3.0,c->GetUxmax(),3.0);
  TArrow *l1=new TArrow(60.0,canv1->GetUymin(),60.0,maxYvalue + 1.5*maxYvalue,0.1,"<|");
  TArrow *l2=new TArrow(120.0,canv1->GetUymin(),120.0,maxYvalue + 1.5*maxYvalue,0.1,"<|");
  TArrow *l3=new TArrow(maxRegion,canv1->GetUymin(),maxRegion,maxYvalue + 1.5*maxYvalue,0.1,"<|");
  l1->SetLineColor(kBlue);
  l2->SetLineColor(kBlue);
  l3->SetLineColor(kRed); 
  l1->SetLineStyle(9);
  l2->SetLineStyle(9);
  l3->SetLineStyle(9);  
  l1->SetLineWidth(2);
  l2->SetLineWidth(2);
  l3->SetLineWidth(2);
  l1->Draw("same");
  l2->Draw("same");
  l3->Draw("same");

  if (LogScale)canv2->cd()->SetLogy();
  canv2->Draw();
  canv2->Update();

  // ------> Third Canvas <------

  TCanvas *canv3 =new TCanvas("canv3","canv3",1000,500);
  canv3->Divide(2,1);

  TCanvas *canv4 =new TCanvas("canv4","canv4",800,600);
  TPad *c4pada = new TPad("c3pada","The pad",0.,0.24,1.,1.);
  TPad *c4padb = new TPad("c3padb","The pad",0.,0.01,1.,0.3);  
  c4pada->Draw();
  c4padb->Draw();
  c4padb->SetBottomMargin(0.35);

  TObjArray *mc = new TObjArray(3);        // MC histograms are put in this array
  mc->Add(h_tt);
  mc->Add(h_Z);
 
  // fit computation
  Double_t fmc1=0., fmc1err=0.;
  Double_t fmc2=0., fmc2err=0.;
  
  TFractionFitter* fit = new TFractionFitter(histos[0],mc,"Q");
  Int_t status = fit->Fit();  

  if (status == 0) {                       // check on fit status
     
    fit->GetResult(0,fmc1,fmc1err);
    fit->GetResult(1,fmc2,fmc2err);
    
    TString output = Form("#chi^{2}=%5.3f ndf=%3d [P=%5.3f]",fit->GetChisquare(),fit->GetNDF(),fit->GetProb());

    canv3->cd(1)->SetLogy();
    histos[0]->Draw("Ep");
    h_tt->Draw("HISTSAME");
    h_Z->Draw("HISTSAME");
    histos[0]->Draw("Epsame");
    CMSPrel(2.2,_decaychannel);

    // canvas of results
    canv3->cd(2)->SetLogy();

    histos[0]->Draw("Ep");
    TH1F* result = (TH1F*) fit->GetPlot();
    result->SetLineWidth(2.);
    
    TH1F* resultTop=(TH1F*) fit->GetMCPrediction(0);
    TH1F *resultTop_new=(TH1F*)resultTop->Clone();
    resultTop_new->Scale(fmc1*histos[0]->Integral()/resultTop->GetSumOfWeights());
    makeNiceHistoStyle(resultTop);
    makeNiceHistoStyle(resultTop_new);
    
    TH1F* resultZ=(TH1F*) fit->GetMCPrediction(1);
    TH1F *resultZ_new=(TH1F*)resultZ->Clone();
    resultZ_new->Scale(fmc2*histos[0]->Integral()/resultZ->GetSumOfWeights());
    makeNiceHistoStyle(resultZ);
    makeNiceHistoStyle(resultZ_new);

    TObjArray *fractions = new TObjArray(4);   
    fractions->Add(histos[0]);
    fractions->Add(resultTop_new);
    fractions->Add(resultZ_new);

    THStack *hs3=new THStack("hs3","test stacked histograms");
    hs3->Add(resultTop_new);
    hs3->Add(resultZ_new);

    hs3->Draw("HISTsame");
    histos[0]->Draw("Epsame");
    CMSPrel(2.2,_decaychannel);

    //TString LegLabels2[4]={"CMS Data",Form("f_{b}: %5.3f #pm %5.3f",fmc1,fmc1err),Form("f_{c}: %5.3f #pm %5.3f",fmc2,fmc2err),Form("f_{l}: %5.3f #pm %5.3f",fmc3,fmc3err)};
    
    // canvas of comparison and result 
    canv4->cd();
    c4pada->cd()->SetLogy();
    histos[0]->Draw("Ep");
    result->Draw("same");
    CMSPrel(2.2,_decaychannel);
    
    TLatex *ll = new TLatex();
    ll->SetTextSize(0.04);
    ll->SetTextFont(42);
    ll->SetLineWidth(2);
    ll->SetNDC();
    ll->DrawLatex(0.6,0.68,Form("f_{t#bar{t}} = %5.3f #pm %5.3f",fmc1,fmc1err));
    ll->DrawLatex(0.6,0.62,Form("f_{Z} = %5.3f #pm %5.3f",fmc2,fmc2err));
    ll->DrawLatex(0.6,0.50,output);

    c4padb->cd();
    TH1F* hratioDATAMC = takeTheRatio(histos[0],result);
    hratioDATAMC->Draw("PE1");

  }

  /* 
     calculation of fttbar from https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=176106
    
     the formula is given by:
     N_tt(in) =  ( (R_tt)/(R_Z -R_tt) ) * (R_Z * N_obs(out) - N_obs(in) ) 

     R_tt = N_tt(in)/N_tt(out)      R_Z = N_Z(in)/N_Z(out)
   
  */

  Float_t minSignalRegion=60.;
  Float_t maxSignalRegion=120.;

  Int_t _firstbin      = histos[0]->FindBin(minSignalRegion);
  Int_t _lastbin       = histos[0]->FindBin(maxSignalRegion);
  Int_t _verylastbin   = histos[0]->FindBin(maxRegion);

  // counts and squared counts
  Float_t N_tt_in(0.), N_tt_out(0.), N_Z_in(0.), N_Z_out(0.), N_obs_in(0.), N_obs_out(0.);
  Float_t sq_N_tt_in(0.), sq_N_tt_out(0.), sq_N_Z_in(0.), sq_N_Z_out(0.);

  // errors on the counts and the squared counts
  Float_t err_N_tt_in(0.), err_N_tt_out(0.), err_N_Z_in(0.), err_N_Z_out(0.), err_N_obs_in(0.), err_N_obs_out(0.);
  Float_t err_R_tt(0.), err_R_Z(0.);

  if(useMyAlgo){
     N_tt_in     = GetTheIntegral(h_tt,_firstbin,_lastbin); 
     N_tt_out    = GetTheIntegral(h_tt,_lastbin,_verylastbin);
     
     sq_N_tt_in  = GetTheIntegral(h_tt,_firstbin,_lastbin,true); 
     sq_N_tt_out = GetTheIntegral(h_tt,_lastbin,_verylastbin,true);
  } else {
     N_tt_in  =  h_tt->Integral(_firstbin,_lastbin);
     N_tt_out =  h_tt->Integral(_lastbin,_verylastbin);
  }

  //err_N_tt_in  = getWeightedBinomialError(files[1],N_tt_in,sq_N_tt_in);
  //err_N_tt_out = getWeightedBinomialError(files[1],N_tt_out,sq_N_tt_out);

  err_N_tt_in  = TMath::Sqrt(N_tt_in);
  err_N_tt_out = TMath::Sqrt(N_tt_out);

  Float_t R_tt = N_tt_in/N_tt_out;
  err_R_tt = propErrRatio(N_tt_in,err_N_tt_in,N_tt_out,err_N_tt_out);

  pair<Float_t,Float_t> R_tt_pair = make_pair(R_tt,err_R_tt);
  
  if(useMyAlgo){
    N_Z_in  = GetTheIntegral(h_Z,_firstbin,_lastbin);   
    N_Z_out = GetTheIntegral(h_Z,_lastbin,_verylastbin);
    
    sq_N_Z_in  = GetTheIntegral(h_Z,_firstbin,_lastbin,true); 
    sq_N_Z_out = GetTheIntegral(h_Z,_lastbin,_verylastbin,true);
  } else {
    N_Z_in  = h_Z->Integral(_firstbin,_lastbin);
    N_Z_out = h_Z->Integral(_lastbin,_verylastbin);
  }


  // err_N_Z_in  = getWeightedBinomialError(files[2],N_Z_in,sq_N_Z_in);
  // err_N_Z_out = getWeightedBinomialError(files[2],N_Z_out,sq_N_Z_out);

  err_N_Z_in  = TMath::Sqrt(N_Z_in); 
  err_N_Z_out =	TMath::Sqrt(N_Z_out);

  Float_t R_Z = N_Z_in/N_Z_out;
  err_R_Z = propErrRatio(N_Z_in,err_N_Z_in,N_Z_out,err_N_Z_out);

  pair<Float_t,Float_t> R_Z_pair = make_pair(R_Z,err_R_Z);

  if(useMyAlgo){
    N_obs_in  = GetTheIntegral(histos[0],_firstbin,_lastbin);    
    N_obs_out = GetTheIntegral(histos[0],_lastbin,_verylastbin);
  } else {
    N_obs_in  =  histos[0]->Integral(_firstbin,_lastbin);
    N_obs_out =  histos[0]->Integral(_lastbin,_verylastbin);
  }

  err_N_obs_in  = Sqrt(N_obs_in);
  err_N_obs_out = Sqrt(N_obs_out);

  pair<Float_t,Float_t> N_obs_in_pair  = make_pair(N_obs_in,err_N_obs_in);
  pair<Float_t,Float_t> N_obs_out_pair = make_pair(N_obs_out,err_N_obs_out);
    
  // final calculation of fttbar

  Float_t N_tt_in_Meas = ((R_tt)/(R_Z - R_tt))*(R_Z * N_obs_out - N_obs_in);
  Float_t err_N_tt_in_Meas = CalculateError(R_tt_pair,R_Z_pair,N_obs_in_pair,N_obs_out_pair);

  Float_t err_f_tt = propErrRatio(N_tt_in_Meas,err_N_tt_in_Meas,N_obs_in,err_N_obs_in);

  //*************************
  // Printing of the factors
  //*************************

  canv2->cd();
  TLatex *report = new TLatex();
  report->SetTextSize(0.035);
  report->SetTextFont(42);
  report->SetLineWidth(2);
  report->SetNDC();
  report->DrawLatex(0.67,0.71,Form("N^{TOT}_{DATA}=%.f",N_obs_in+N_obs_out));
  report->DrawLatex(0.67,0.65,Form("N^{in}_{DATA}=%.f N^{out}_{DATA}=%.f",N_obs_in,N_obs_out));
  report->DrawLatex(0.67,0.59,Form("N^{in}_{Z}=%.f N^{out}_{Z}=%.f",N_Z_in,N_Z_out));
  report->DrawLatex(0.67,0.53,Form("N^{in}_{t#bar{t}}=%.f N^{out}_{t#bar{t}}=%.f",N_tt_in,N_tt_out));
  report->DrawLatex(0.67,0.46,Form("f^{est}_{t#bar{t}}=%.3f f^{pred}_{t#bar{t}}=%.3f",(N_tt_in_Meas/N_obs_in),(N_tt_in/(N_Z_in+N_tt_in))));
  
  canv1->SaveAs(foldername+"_"+theALGOWP+"_"+histoname+".png");

  using namespace std;
  cout.precision(3);

  // for screen
  _decaychannel.ReplaceAll("#rightarrow","->");
  _decaychannel.ReplaceAll("#","");

  cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"| Printing report: "<<foldername<<" "<<theALGOWP<<" "<<_decaychannel<<endl;
  cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"| N_tt_in         = "<< setw(10) << fixed <<N_tt_in   <<" +/- "<< err_N_tt_in   <<endl;
  cout<<"| N_tt_out        = "<< setw(10) << fixed <<N_tt_out  <<" +/- "<< err_N_tt_out  <<endl;
  cout<<"| R_tt            = "<< setw(10) << fixed <<R_tt      <<" +/- "<< err_R_tt      <<endl;
  cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl; 
  cout<<"| N_Z_in          = "<< setw(10) << fixed <<N_Z_in    <<" +/- "<< err_N_Z_in    <<endl;
  cout<<"| N_Z_out         = "<< setw(10) << fixed <<N_Z_out   <<" +/- "<< err_N_Z_out   <<endl;
  cout<<"| R_Z             = "<< setw(10) << fixed <<R_Z       <<" +/- "<< err_R_Z       <<endl;
  cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"| N_obs_in        = "<< setw(10) << fixed <<N_obs_in  <<" +/- "<< err_N_obs_in  <<endl;
  cout<<"| N_obs_out       = "<< setw(10) << fixed <<N_obs_out <<" +/- "<< err_N_obs_out <<endl;
  cout<<"-----------------------------------------------------"<<endl;
  cout<<"| N_tt_in_Meas    = "<< setw(10) << fixed <<N_tt_in_Meas <<" +/- "<<err_N_tt_in_Meas<<endl;
  cout<<"-----------------------------------------------------"<<endl;
  cout<<"| f_ttbar         = "<< setw(10) << fixed <<N_tt_in_Meas/N_obs_in <<" +/- "<<err_f_tt<<endl;
  cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

  // for outfile
  outfile_<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  outfile_<<"| Printing report: "<<foldername<<" "<<theALGOWP<<" "<<_decaychannel<<endl;
  outfile_<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  outfile_<<"| N_tt_in         = "<< setw(10) << fixed <<N_tt_in   <<" +/- "<< err_N_tt_in   <<endl;
  outfile_<<"| N_tt_out        = "<< setw(10) << fixed <<N_tt_out  <<" +/- "<< err_N_tt_out  <<endl;
  outfile_<<"| R_tt            = "<< setw(10) << fixed <<R_tt      <<" +/- "<< err_R_tt      <<endl;
  outfile_<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl; 
  outfile_<<"| N_Z_in          = "<< setw(10) << fixed <<N_Z_in    <<" +/- "<< err_N_Z_in    <<endl;
  outfile_<<"| N_Z_out         = "<< setw(10) << fixed <<N_Z_out   <<" +/- "<< err_N_Z_out   <<endl;
  outfile_<<"| R_Z             = "<< setw(10) << fixed <<R_Z       <<" +/- "<< err_R_Z       <<endl;
  outfile_<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  outfile_<<"| N_obs_in        = "<< setw(10) << fixed <<N_obs_in  <<" +/- "<< err_N_obs_in  <<endl;
  outfile_<<"| N_obs_out       = "<< setw(10) << fixed <<N_obs_out <<" +/- "<< err_N_obs_out <<endl;
  outfile_<<"-----------------------------------------------------"<<endl;
  outfile_<<"| N_tt_in_Meas    = "<< setw(10) << fixed <<N_tt_in_Meas <<" +/- "<<err_N_tt_in_Meas<<endl;
  outfile_<<"-----------------------------------------------------"<<endl;
  outfile_<<"| f_ttbar         = "<< setw(10) << fixed <<N_tt_in_Meas/N_obs_in <<" +/- "<<err_f_tt<<endl;
  outfile_<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
 
  outfile_.close();

}

////////////////////////////////////////////////////////////////////
void CMSPrel(Float_t Lumi,TString _decaychannel){
  
  TLatex *latexLabel = new TLatex();
  latexLabel->SetTextSize(0.035);
  latexLabel->SetTextFont(42);
  latexLabel->SetLineWidth(2);
  latexLabel->SetNDC();

  latexLabel->DrawLatex(0.096,0.965,"CMS VHF Preliminary");
  latexLabel->DrawLatex(0.85,0.965,"#sqrt{s} = 7 TeV"); 
  latexLabel->DrawLatex(0.68,0.89,"anti-k_{T} (R = 0.5) PF Jets ");
  latexLabel->DrawLatex(0.68,0.83,(TString)Form("#int Ldt = %.1f fb^{-1}",Lumi));  
  latexLabel->DrawLatex(0.68,0.77,_decaychannel);
}



////////////////////////////////////////////////////////////////////
Float_t GetTheIntegral(TH1F *histo,Int_t _firstbin, Int_t _lastbin,Bool_t getTheSquare){

  Float_t sum(0.);

  if(!getTheSquare){
    for(Int_t nbin=_firstbin; nbin<_lastbin; nbin++) {
      sum+=histo->GetBinContent(nbin);
      //sum+=Sqrt(histo->GetSumw2()->At(nbin));
    }
  } else {
    for(Int_t nbin=_firstbin; nbin<_lastbin; nbin++) {
      sum+=Power(histo->GetBinContent(nbin),2);
      //sum+=histo->GetSumw2()->At(nbin);
    }
  }
  
  return sum;

}

////////////////////////////////////////////////////////////////////
Float_t CalculateError(pair<Float_t,Float_t> R_tt_pair, 
			pair<Float_t,Float_t> R_Z_pair, 
			pair<Float_t,Float_t> N_obs_in_pair, 
			pair<Float_t,Float_t> N_obs_out_pair){
  
  Float_t _error(0.);
  Float_t _dFdR_tt(0.), _dFdR_Z(0.), _dFdN_obs_in(0.), _dFdN_obs_out(0.);

  /* 
     calculation of fttbar from https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=176106
     
     the formula is given by:
     N_tt(in) =  ( (R_tt)/(R_Z -R_tt) ) * (R_Z * N_obs(out) - N_obs(in) ) 

     R_tt = N_tt(in)/N_tt(out)      R_Z = N_Z(in)/N_Z(out)
     
  */

  // calculation of the partial derivatives
  
  _dFdR_tt      = (R_Z_pair.first*(R_Z_pair.first*N_obs_out_pair.first - N_obs_in_pair.first) ) / Power(R_tt_pair.first - R_Z_pair.first,2);
  _dFdR_Z       = (R_tt_pair.first*(- R_tt_pair.first*N_obs_out_pair.first + N_obs_in_pair.first) ) / Power(R_tt_pair.first - R_Z_pair.first,2);
  _dFdN_obs_in  = (R_tt_pair.first) / (R_tt_pair.first - R_Z_pair.first); 
  _dFdN_obs_out = (R_tt_pair.first*R_Z_pair.first)/(R_tt_pair.first - R_Z_pair.first);
  
  _error = Sqrt( Power(_dFdR_tt*R_tt_pair.second,2) + 
			Power(_dFdR_Z*R_Z_pair.second,2) + 
			Power(_dFdN_obs_in*N_obs_in_pair.second,2) + 
			Power(_dFdN_obs_out*N_obs_out_pair.second,2) );

  return _error;
  
}

////////////////////////////////////////////////////////////////////
Float_t getWeightedBinomialError(TFile *file, Float_t myCounts, Float_t mySquaredCounts){
  
  Int_t nAOD   = ((TH1F*)file->Get("analyzePat/Selevents/SelectedEvts"))->GetBinContent(1);

  TString href = "finaldistros_ssvhem/EventYields/NPileupReweight";

  Float_t avg_sumw  = ((TH1F*)file->Get(href))->GetSumOfWeights()/((TH1F*)file->Get(href))->GetEntries(); 
  Float_t avg_sumw2 = ((TH1F*)file->Get(href))->GetSumw2()->GetSum()/((TH1F*)file->Get(href))->GetEntries();

  Float_t sumPass  = myCounts;
  Float_t sum2Pass = mySquaredCounts;
  Float_t sumFail  = nAOD*avg_sumw - sumPass;
  Float_t sum2Fail = (nAOD*avg_sumw2) - sum2Pass;
  
  Float_t err_eff  = Sqrt(sum2Pass*sumFail*sumFail + sum2Fail*sumPass*sumPass) / (nAOD*avg_sumw*nAOD*avg_sumw);

  cout<<"<------------------>"<<endl;
  cout<<"nAOD:     "<<nAOD<<endl;
  cout<<"avg_sumw: "<<avg_sumw<<endl;
  cout<<"avg_sumw2:"<<avg_sumw2<<endl;
  cout<<"sumPass:  "<<sumPass<<endl;
  cout<<"sum2Pass: "<<sum2Pass<<endl;

  return  nAOD*err_eff;

}

/////////////////////////////////////////////////////////////////////
void setStyle(){
  
  //***************
  //Graphics
  //***************
  
  // For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(600); //Height of canvas
  gStyle->SetCanvasDefW(600); //Width of canvas
  gStyle->SetCanvasDefX(0);   //Position on screen
  gStyle->SetCanvasDefY(0);

  // For the Pad:
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  
  // For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);
  
  // For the title:
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);

  // For the margins:
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.03);
  
}

void makeNiceHistoStyle(TH1F *h, Int_t color){
  if(color!=0) h->SetLineColor(color);
  h->SetFillColor(color);
  h->SetTitleSize(0.09);
  h->SetTitleFont(42); 
  h->SetLineWidth(2);
  h->GetYaxis()->SetTitleOffset(1.);
  h->GetXaxis()->SetTitleOffset(1.);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleFont(42); 
  h->GetYaxis()->SetTitleFont(42);  
  h->GetXaxis()->SetLabelSize(0.05); 
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelFont(42); 
  h->SetStats(false);
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
