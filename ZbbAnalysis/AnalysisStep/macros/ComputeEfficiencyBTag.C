//-----------------------------------------//
//howTo use it
//root [0] .L ComputeEfficiencyBTag.C++
//root [1] ComputeEfficiencyBTag("fileC.root","fileB.root","theTagger (SSVHEM or SSVHPT or TCHEL or TCHEM or TCHPT)",isWeighted (true or false))
//------------------------------------------//

#include <string>
#include <iostream>
#include <Rtypes.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"

//**********************************************************************
TLegend* MakeTLegend(TObjArray *array,TString *LegLabels){
//**********************************************************************
  TLegend* leg = new TLegend(0.61,0.69,0.88,0.94);
  leg->SetFillColor(10);
  leg->SetShadowColor(10);
  leg->SetTextFont(42);
  for(Int_t j=0; j < array->GetSize(); j++){
    if(j==2){
      leg->AddEntry(static_cast<TH1*>(array->At(j)),LegLabels[j],"L");
    } else {
      leg->AddEntry(static_cast<TH1*>(array->At(j)),LegLabels[j],"P");
    }
  }
  return leg;
}

//**********************************************************************
void MakeNiceLatexLabel(TLatex *latexLabel){  
//**********************************************************************
  latexLabel->SetTextSize(0.05);
  latexLabel->SetTextFont(62);
  latexLabel->SetLineWidth(2);
  latexLabel->SetNDC();
  latexLabel->DrawLatex(0.20, 0.89, "CMS Z+b WG");
  latexLabel->DrawLatex(0.20, 0.84, "Simulation");
}
//**********************************************************************
void MakeNiceHistoStyle(TH1F *hist){
//**********************************************************************
  hist->SetTitleSize(0.09);
  hist->SetTitleFont(42); 
  hist->SetLineWidth(2);
  hist->GetXaxis()->SetNdivisions(305);
  hist->GetYaxis()->SetTitleOffset(1.11);
  hist->GetXaxis()->SetTitleOffset(1);
  hist->GetXaxis()->CenterTitle(kTRUE);
  hist->GetYaxis()->CenterTitle(kTRUE);
  hist->GetYaxis()->SetTitleSize(0.065);
  hist->GetXaxis()->SetTitleSize(0.065);
  hist->GetXaxis()->SetTitleFont(42); 
  hist->GetYaxis()->SetTitleFont(42);  
  hist->GetXaxis()->SetLabelSize(0.055); 
  hist->GetYaxis()->SetLabelSize(0.055);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
}

/*************************************************************************************/
TH1F* GetTH1Efficiency(const TH1F* h_passed, const TH1F* h_total, std::string newname)
/*************************************************************************************/
{
  TH1F *h_eff = (TH1F*)h_total->Clone(newname.c_str()); 
  h_eff->Reset();
  Int_t nbin=h_total->GetNbinsX();
  double passed, total, eff, erreff;
  double sum_passed(0), sum_total(0);
  for (Int_t ibin=1; ibin<=nbin; ibin++){
    total=h_total->GetBinContent(ibin);
    passed=h_passed->GetBinContent(ibin);
    //std::cout << ibin << " , " << passed <<  " / " << total <<std::endl; 
    if ( ibin>1) {
      sum_total+=total;
      sum_passed+=passed;
    }
    if ( total>0 ) {
      eff=passed/total;
      erreff=TMath::Sqrt(eff*(1-eff)/total);   
      h_eff->SetBinContent(ibin,eff);
      h_eff->SetBinError(ibin,erreff);
    }       
    //    std::cout << ibin << " , " << eff <<  " / " << erreff <<std::endl; 
  }

  // now grand-total
  eff=0; erreff=0;
  if ( sum_total>0 ) {
      eff=sum_passed/sum_total;
      erreff=TMath::Sqrt(eff*(1-eff)/sum_total)*(sum_passed/h_passed->GetEntries()); //implementing the correcting factor
  }

  // printing efficiencies
  TString namecase=newname;
  TObjArray *infos =namecase.Tokenize("_");
  TString theEff = (((TObjString*)(infos->At(2)))->String()).ReplaceAll("Over","/") +" "+((TObjString*)(infos->At(4)))->String();  
  TString efftype = namecase.Contains("b1")  ? " b-eff: " : " c-eff: ";
  std::cout << theEff << efftype << eff <<"+/-"<< erreff <<std::endl;
  return (h_eff);
}

void ComputeEfficiencyBTag(std::string InputFileC, std::string InputFileB, std::string tagger,Bool_t isWeighted){
  
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
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
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.02);
  
  TString foutname="TH1F_"+tagger+".root";
  TFile *fout = new TFile(foutname,"RECREATE");

  TString title_b_eff = "p_{T} leading b in the event (recoJet) all; p_{T} leading b-jet (GeV);"+tagger+" b-tag efficiency #epsilon_{b}";
  TString title_c_eff = "p_{T} leading c in the event (recoJet) all; p_{T} leading c-jet (GeV);"+tagger+" b-tag efficiency #epsilon_{b}";

  TFile *finC = new TFile(InputFileC.c_str(),"READ");
  TFile *finB = new TFile(InputFileB.c_str(),"READ");
  
  TString tagpart, wgtpart;

  // weight part
  if(isWeighted){
    wgtpart="_wgt";
  } else {
    wgtpart="_unw";
  }

  // tagger part
  if(tagger=="SSVHEM"){
    tagpart="_ssvhem";
  } else if(tagger=="SSVHPT"){
    tagpart="_ssvhpt";
  } else if(tagger=="TCHEL"){
    tagpart="_tchel";
  } else if(tagger=="TCHEM"){
    tagpart="_tchem";
  } else if(tagger=="TCHPT"){
    tagpart="_tchpt";
  } else {
    std::cout<<"Using unkown tagger"<<endl;
    return;
  }

  TString myFinalDistroString = "finaldistros"+tagpart;

  TH1F *h_ptb1_brl_zll  =(TH1F*)finB->Get(myFinalDistroString+"/ZLL/ZLL_RECO_pt_b1_etaBrl"); 
  h_ptb1_brl_zll ->SetTitle(title_b_eff);
  TH1F *h_ptb1_fwd_zll  =(TH1F*)finB->Get(myFinalDistroString+"/ZLL/ZLL_RECO_pt_b1_etaFwd"); 
  h_ptb1_fwd_zll ->SetTitle(title_b_eff);
  TH1F *h_ptb1_brl_goodjet=(TH1F*)finB->Get(myFinalDistroString+"/GoodJet/GoodJet_RECO_pt_b1_etaBrl"); 
  h_ptb1_brl_goodjet ->SetTitle(title_b_eff);
  TH1F *h_ptb1_fwd_goodjet=(TH1F*)finB->Get(myFinalDistroString+"/GoodJet/GoodJet_RECO_pt_b1_etaFwd"); 
  h_ptb1_fwd_goodjet ->SetTitle(title_b_eff);
  TH1F *h_ptb1_brl_btag  =(TH1F*)finB->Get(myFinalDistroString+"/JetBTag/JetBTag_RECO_pt_b1_etaBrl"); 
  h_ptb1_brl_btag ->SetTitle(title_b_eff);
  TH1F *h_ptb1_fwd_btag  =(TH1F*)finB->Get(myFinalDistroString+"/JetBTag/JetBTag_RECO_pt_b1_etaFwd"); 
  h_ptb1_fwd_btag ->SetTitle(title_b_eff);

  TH1F *h_ptc1_brl_zll  =(TH1F*)finC->Get(myFinalDistroString+"/ZLL/ZLL_RECO_pt_c1_etaBrl"); 
  h_ptc1_brl_zll->SetTitle(title_c_eff);
  TH1F *h_ptc1_fwd_zll  =(TH1F*)finC->Get(myFinalDistroString+"/ZLL/ZLL_RECO_pt_c1_etaFwd"); 
  h_ptc1_fwd_zll->SetTitle(title_c_eff);
  TH1F *h_ptc1_brl_goodjet=(TH1F*)finC->Get(myFinalDistroString+"/GoodJet/GoodJet_RECO_pt_c1_etaBrl"); 
  h_ptc1_brl_goodjet->SetTitle(title_c_eff);
  TH1F *h_ptc1_fwd_goodjet=(TH1F*)finC->Get(myFinalDistroString+"/GoodJet/GoodJet_RECO_pt_c1_etaFwd"); 
  h_ptc1_fwd_goodjet->SetTitle(title_c_eff);
  TH1F *h_ptc1_brl_btag =(TH1F*)finC->Get(myFinalDistroString+"/JetBTag/JetBTag_RECO_pt_c1_etaBrl"); 
  h_ptc1_brl_btag->SetTitle(title_c_eff);
  TH1F *h_ptc1_fwd_btag  =(TH1F*)finC->Get(myFinalDistroString+"/JetBTag/JetBTag_RECO_pt_c1_etaFwd"); 
  h_ptc1_fwd_btag->SetTitle(title_c_eff);

  // full eta range
  TH1F *h_ptb1_all_zll    =(TH1F*)h_ptb1_brl_zll->Clone("h_ptb1_all_zll");    
  h_ptb1_all_zll->Add(h_ptb1_fwd_zll);
  TH1F *h_ptb1_all_goodjet=(TH1F*)h_ptb1_brl_goodjet->Clone("h_ptb1_all_goodjet");
  h_ptb1_all_goodjet->Add(h_ptb1_fwd_goodjet);
  TH1F *h_ptb1_all_btag   =(TH1F*)h_ptb1_brl_btag->Clone("h_ptb1_all_btag");
  h_ptb1_all_btag->Add(h_ptb1_fwd_btag);
  
  TH1F *h_ptc1_all_zll  =(TH1F*)h_ptc1_brl_zll->Clone("h_ptc1_all_zll");
  h_ptc1_all_zll->Add(h_ptc1_fwd_zll);
  TH1F *h_ptc1_all_goodjet=(TH1F*)h_ptc1_brl_goodjet->Clone("h_ptc1_all_goodjet");
  h_ptc1_all_goodjet->Add(h_ptc1_fwd_goodjet);
  TH1F *h_ptc1_all_btag   =(TH1F*)h_ptc1_brl_btag->Clone("h_ptc1_all_btag");
  h_ptc1_all_btag->Add(h_ptc1_fwd_btag);

  // values for eff_b from BTV-11-001
  TH1F *h_eff_ptb1_pas = (TH1F*)h_ptb1_brl_zll->Clone("h_eff_ptb1_pas");
  h_eff_ptb1_pas->Reset();
  const int nbin(11);
  // SSVHEM
  double eff_pas_ssvhem[nbin]={0.,0.137526,0.330398,0.472956,0.571908,0.625577,0.657442,0.680922, 0.6826,0.6826,0.652411};
  // SSVHPT
  double eff_pas_ssvhpt[nbin]={0., 0.0181275, 0.11992, 0.235657, 0.333267, 0.389044, 0.432271, 0.462948, 0.48247, 0.478287, 0.457371};
  /*
  // TCHEL
  double eff_pas_tchel[nbin]={0., 0.0181275, 0.11992, 0.235657, 0.333267, 0.389044, 0.432271, 0.462948, 0.48247, 0.478287, 0.457371};
  // TCHPT
  double eff_pas_tchem[nbin]={0., 0.0181275, 0.11992, 0.235657, 0.333267, 0.389044, 0.432271, 0.462948, 0.48247, 0.478287, 0.457371};
  // TCHPT
  double eff_pas_tchpt[nbin]={0., 0.0181275, 0.11992, 0.235657, 0.333267, 0.389044, 0.432271, 0.462948, 0.48247, 0.478287, 0.457371};
  */

  for (Int_t ibin=1; ibin<=nbin; ibin++){
      if(tagger=="SSVHEM"){
	h_eff_ptb1_pas->SetBinContent(ibin,eff_pas_ssvhem[ibin-1]);  
      }  else if(tagger=="SSVHPT"){
	h_eff_ptb1_pas->SetBinContent(ibin,eff_pas_ssvhpt[ibin-1]);  
      }
      /*
	else if(tagger=="TCHEL"){
	h_eff_ptb1_pas->SetBinContent(ibin,eff_pas_tchel[ibin-1]);  
	} else if(tagger=="TCHEM"){
	h_eff_ptb1_pas->SetBinContent(ibin,eff_pas_tchem[ibin-1]);  
	} else if(tagger=="TCHPT"){
	h_eff_ptb1_pas->SetBinContent(ibin,eff_pas_tchpt[ibin-1]);  
	}
      */
      else {
	std::cout<<"Using unkown tagger"<<endl;
	return;
      }   
  }
  
  h_eff_ptb1_pas->SetLineColor(kMagenta);
  MakeNiceHistoStyle(h_eff_ptb1_pas);

  // now compute efficiencies 
  // 1. wrt ZLL step
  TH1F* h_eff_bTagOverZLL_ptb1_brl = GetTH1Efficiency(h_ptb1_brl_btag,h_ptb1_brl_zll,"h_eff_bTagOverZLL_ptb1_brl");
  TH1F* h_eff_bTagOverZLL_ptb1_fwd = GetTH1Efficiency(h_ptb1_fwd_btag,h_ptb1_fwd_zll,"h_eff_bTagOverZLL_ptb1_fwd");
  TH1F* h_eff_bTagOverZLL_ptb1_all = GetTH1Efficiency(h_ptb1_all_btag,h_ptb1_all_zll,"h_eff_bTagOverZLL_ptb1_all");

  TH1F* h_eff_bTagOverZLL_ptc1_brl = GetTH1Efficiency(h_ptc1_brl_btag,h_ptc1_brl_zll,"h_eff_bTagOverZLL_ptc1_brl");
  TH1F* h_eff_bTagOverZLL_ptc1_fwd = GetTH1Efficiency(h_ptc1_fwd_btag,h_ptc1_fwd_zll,"h_eff_bTagOverZLL_ptc1_fwd");
  TH1F* h_eff_bTagOverZLL_ptc1_all = GetTH1Efficiency(h_ptc1_all_btag,h_ptc1_all_zll,"h_eff_bTagOverZLL_ptc1_all");

  h_eff_bTagOverZLL_ptb1_all->SetMarkerStyle(kOpenCircle); 
  h_eff_bTagOverZLL_ptb1_all->SetMarkerColor(kRed);     
  h_eff_bTagOverZLL_ptb1_all->SetLineColor(kRed);     
  h_eff_bTagOverZLL_ptb1_all->GetYaxis()->SetRangeUser(0.,1.55);     
  MakeNiceHistoStyle(h_eff_bTagOverZLL_ptb1_all);

  h_eff_bTagOverZLL_ptb1_brl->SetMarkerStyle(kOpenCircle); 
  h_eff_bTagOverZLL_ptb1_brl->SetMarkerColor(kRed);     
  h_eff_bTagOverZLL_ptb1_brl->SetLineColor(kRed);     
  h_eff_bTagOverZLL_ptb1_brl->GetYaxis()->SetRangeUser(0.,1.55);     
  MakeNiceHistoStyle(h_eff_bTagOverZLL_ptb1_brl);

  h_eff_bTagOverZLL_ptb1_fwd->SetMarkerStyle(kOpenCircle); 
  h_eff_bTagOverZLL_ptb1_fwd->SetMarkerColor(kRed);       
  h_eff_bTagOverZLL_ptb1_fwd->SetLineColor(kRed);       
  h_eff_bTagOverZLL_ptb1_fwd->GetYaxis()->SetRangeUser(0.,1.55);     
  MakeNiceHistoStyle(h_eff_bTagOverZLL_ptb1_fwd);

  h_eff_bTagOverZLL_ptc1_all->SetMarkerStyle(kOpenCircle); 
  h_eff_bTagOverZLL_ptc1_all->SetMarkerColor(kRed);     
  h_eff_bTagOverZLL_ptc1_all->SetLineColor(kRed);     
  h_eff_bTagOverZLL_ptc1_all->GetYaxis()->SetRangeUser(0.,0.45);       
  MakeNiceHistoStyle(h_eff_bTagOverZLL_ptc1_all);

  h_eff_bTagOverZLL_ptc1_brl->SetMarkerStyle(kOpenCircle); 
  h_eff_bTagOverZLL_ptc1_brl->SetMarkerColor(kRed);     
  h_eff_bTagOverZLL_ptc1_brl->SetLineColor(kRed);     
  h_eff_bTagOverZLL_ptc1_brl->GetYaxis()->SetRangeUser(0.,0.45);       
  MakeNiceHistoStyle(h_eff_bTagOverZLL_ptc1_brl);

  h_eff_bTagOverZLL_ptc1_fwd->SetMarkerStyle(kOpenCircle); 
  h_eff_bTagOverZLL_ptc1_fwd->SetMarkerColor(kRed);     
  h_eff_bTagOverZLL_ptc1_fwd->SetLineColor(kRed);     
  h_eff_bTagOverZLL_ptc1_fwd->GetYaxis()->SetRangeUser(0.,0.45);       
  MakeNiceHistoStyle(h_eff_bTagOverZLL_ptc1_fwd);

  // 2. wrt GoodJet step
  TH1F* h_eff_bTagOverGoodJet_ptb1_brl = GetTH1Efficiency(h_ptb1_brl_btag,h_ptb1_brl_goodjet,"h_eff_bTagOverGoodJet_ptb1_brl");
  TH1F* h_eff_bTagOverGoodJet_ptb1_fwd = GetTH1Efficiency(h_ptb1_fwd_btag,h_ptb1_fwd_goodjet,"h_eff_bTagOverGoodJet_ptb1_fwd");
  TH1F* h_eff_bTagOverGoodJet_ptb1_all = GetTH1Efficiency(h_ptb1_all_btag,h_ptb1_all_goodjet,"h_eff_bTagOverGoodJet_ptb1_all");

  TH1F* h_eff_bTagOverGoodJet_ptc1_brl = GetTH1Efficiency(h_ptc1_brl_btag,h_ptc1_brl_goodjet,"h_eff_bTagOverGoodJet_ptc1_brl");
  TH1F* h_eff_bTagOverGoodJet_ptc1_fwd = GetTH1Efficiency(h_ptc1_fwd_btag,h_ptc1_fwd_goodjet,"h_eff_bTagOverGoodJet_ptc1_fwd");
  TH1F* h_eff_bTagOverGoodJet_ptc1_all = GetTH1Efficiency(h_ptc1_all_btag,h_ptc1_all_goodjet,"h_eff_bTagOverGoodJet_ptc1_all");

  h_eff_bTagOverGoodJet_ptb1_all->SetMarkerStyle(kOpenSquare); 
  h_eff_bTagOverGoodJet_ptb1_all->SetMarkerColor(kBlue);     
  h_eff_bTagOverGoodJet_ptb1_all->SetLineColor(kBlue);     
  h_eff_bTagOverGoodJet_ptb1_all->GetYaxis()->SetRangeUser(0.,1.55);     
  MakeNiceHistoStyle(h_eff_bTagOverGoodJet_ptb1_all);
 
  h_eff_bTagOverGoodJet_ptb1_brl->SetMarkerStyle(kOpenSquare); 
  h_eff_bTagOverGoodJet_ptb1_brl->SetMarkerColor(kBlue);     
  h_eff_bTagOverGoodJet_ptb1_brl->SetLineColor(kBlue);     
  h_eff_bTagOverGoodJet_ptb1_brl->GetYaxis()->SetRangeUser(0.,1.55);     
  MakeNiceHistoStyle(h_eff_bTagOverGoodJet_ptb1_brl);

  h_eff_bTagOverGoodJet_ptb1_fwd->SetMarkerStyle(kOpenSquare); 
  h_eff_bTagOverGoodJet_ptb1_fwd->SetMarkerColor(kBlue);     
  h_eff_bTagOverGoodJet_ptb1_fwd->SetLineColor(kBlue);       
  h_eff_bTagOverGoodJet_ptb1_fwd->GetYaxis()->SetRangeUser(0.,1.55);     
  MakeNiceHistoStyle(h_eff_bTagOverGoodJet_ptb1_fwd);

  h_eff_bTagOverGoodJet_ptc1_all->SetMarkerStyle(kOpenSquare); 
  h_eff_bTagOverGoodJet_ptc1_all->SetMarkerColor(kBlue);     
  h_eff_bTagOverGoodJet_ptc1_all->SetLineColor(kBlue);     
  h_eff_bTagOverGoodJet_ptc1_all->GetYaxis()->SetRangeUser(0.,0.45);      
  MakeNiceHistoStyle( h_eff_bTagOverGoodJet_ptc1_all); 

  h_eff_bTagOverGoodJet_ptc1_brl->SetMarkerStyle(kOpenSquare); 
  h_eff_bTagOverGoodJet_ptc1_brl->SetMarkerColor(kBlue);     
  h_eff_bTagOverGoodJet_ptc1_brl->SetLineColor(kBlue);     
  h_eff_bTagOverGoodJet_ptc1_brl->GetYaxis()->SetRangeUser(0.,0.45);       
  MakeNiceHistoStyle( h_eff_bTagOverGoodJet_ptc1_brl); 

  h_eff_bTagOverGoodJet_ptc1_fwd->SetMarkerStyle(kOpenSquare);
  h_eff_bTagOverGoodJet_ptc1_fwd->SetMarkerColor(kBlue);     
  h_eff_bTagOverGoodJet_ptc1_fwd->SetLineColor(kBlue);     
  h_eff_bTagOverGoodJet_ptc1_fwd->GetYaxis()->SetRangeUser(0.,0.45);       
  MakeNiceHistoStyle(h_eff_bTagOverGoodJet_ptc1_fwd); 

  TCanvas *c1 = new TCanvas("c1","MC eff b/c",0,0,1200,700);
  c1->SetFillColor(10);
  c1->Divide(3,2);
  
  // latex labels
  TLatex *latexLabel1 = new TLatex();
  TLatex *latexLabel2 = new TLatex();
  TLatex *latexLabel3 = new TLatex();
  
  // b-efficiency
  c1->cd(1); h_eff_bTagOverGoodJet_ptb1_all->Draw("PE"); h_eff_bTagOverZLL_ptb1_all->Draw("PESAME"); h_eff_ptb1_pas->Draw("SAME");
  
  TObjArray *arrayHistos = new TObjArray();
  arrayHistos->Expand(3);
  arrayHistos->Add(h_eff_bTagOverGoodJet_ptb1_all);  arrayHistos->Add(h_eff_bTagOverZLL_ptb1_all);  arrayHistos->Add(h_eff_ptb1_pas);  
  
  // legends
  TString LegLabels[3] = {"N^{Z+b}_{ll+tag}/N^{Z+b}_{ll}","N^{Z+b}_{ll+tag}/N^{Z+b}_{ll+j}","#epsilon_{b}(PAS)"};
  TLegend *legend = MakeTLegend(arrayHistos,LegLabels);
  legend->Draw("same");
	
  MakeNiceLatexLabel(latexLabel1);  
  latexLabel1->DrawLatex(0.20,0.79,"|#eta_{jet}|<2.1");
  c1->cd(2); h_eff_bTagOverGoodJet_ptb1_brl->Draw("PE"); h_eff_bTagOverZLL_ptb1_brl->Draw("PESAME"); h_eff_ptb1_pas->Draw("SAME");
  MakeNiceLatexLabel(latexLabel2);  
  latexLabel2->DrawLatex(0.20,0.79,"|#eta_{jet}|<1.2");
  legend->Draw("same");
  c1->cd(3); h_eff_bTagOverGoodJet_ptb1_fwd->Draw("PE"); h_eff_bTagOverZLL_ptb1_fwd->Draw("PESAME"); h_eff_ptb1_pas->Draw("SAME");
  MakeNiceLatexLabel(latexLabel3);  
  latexLabel3->DrawLatex(0.20,0.79,"1.2<|#eta_{jet}|<2.1");
  legend->Draw("same");
 
  // c-efficiency
  c1->cd(4); h_eff_bTagOverGoodJet_ptc1_all->Draw("PE"); h_eff_bTagOverZLL_ptc1_all->Draw("PESAME");

  TObjArray *arrayHistos2 = new TObjArray();
  arrayHistos2->Expand(2);
  arrayHistos2->Add(h_eff_bTagOverGoodJet_ptc1_all);  arrayHistos2->Add(h_eff_bTagOverZLL_ptc1_all); 
  
  // legends
  TString LegLabels2[2] = {"N^{Z+c}_{ll+tag}/N^{Z+c}_{ll}","N^{Z+c}_{ll+tag}/N^{Z+b}_{ll+j}"};
  TLegend *legend2 = MakeTLegend(arrayHistos2,LegLabels2);
  legend2->Draw("same");

  MakeNiceLatexLabel(latexLabel1);  
  latexLabel1->DrawLatex(0.20,0.79,"|#eta_{jet}|<2.1");
  c1->cd(5); h_eff_bTagOverGoodJet_ptc1_brl->Draw("PE"); h_eff_bTagOverZLL_ptc1_brl->Draw("PESAME");
  MakeNiceLatexLabel(latexLabel2);  
  latexLabel1->DrawLatex(0.20,0.79,"|#eta_{jet}|<1.2");
  legend2->Draw("same");
  c1->cd(6); h_eff_bTagOverGoodJet_ptc1_fwd->Draw("PE"); h_eff_bTagOverZLL_ptc1_fwd->Draw("PESAME");
  MakeNiceLatexLabel(latexLabel3);  
  latexLabel1->DrawLatex(0.20,0.79,"1.2<|#eta_{jet}|<2.1");
  legend2->Draw("same");

  c1->SaveAs("c1"+tagpart+wgtpart+".root");
  c1->SaveAs("c1"+tagpart+wgtpart+".png");
  
  // write ouput file
  fout->cd();
  h_eff_bTagOverGoodJet_ptb1_all->Write(); 
  h_eff_bTagOverGoodJet_ptb1_brl->Write(); 
  h_eff_bTagOverGoodJet_ptb1_fwd->Write(); 
  h_eff_bTagOverGoodJet_ptc1_all->Write();
  h_eff_bTagOverGoodJet_ptc1_brl->Write();
  h_eff_bTagOverGoodJet_ptc1_fwd->Write();
  h_eff_ptb1_pas->Write();
  fout->Close();

  return;
}


