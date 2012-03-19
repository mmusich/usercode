#include "JetByJetComparisonHistos.h"
#include "JetInfo.h"
#include <TColor.h>
#include <TCanvas.h>
#include <TClass.h>
#include <TObjArray.h>
#include <iostream>
#include <TStyle.h>
#include <TMath.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TPaveStats.h>
#include <tdrStyle.C>

void cmsPrel(const double& intLumi,TString obj1name_,TString obj2name_,Bool_t isMC) {
  
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.045);
  latex->SetTextFont(42); //22

  latex->SetTextAlign(13);
  if(isMC){
    latex->DrawLatex(0.2,0.99,"CMS - MC Simulation 2011, #sqrt{s} = 7 TeV");
  } else {
    latex->DrawLatex(0.12, 0.99, Form("CMS Preliminary 2011,     #sqrt{s} = 7 TeV,  L = %.2g pb^{ -1}",intLumi));
  }

  TLatex *latex2 = new TLatex();
  latex2->SetNDC();
  latex2->SetTextSize(0.045);
  latex2->SetTextFont(72); //22
  latex2->SetTextColor(); //22
  latex2->DrawLatex(0.20,0.85,obj1name_+" vs "+obj2name_);
}


void DrawAllHistos(TString filename,Bool_t isMC,TString format){
  
  //  setTDRStyle("logredblue");

  std::vector<TH1F*> h1vec;
  std::vector<TH2F*> h2vec;
  
  h1vec.clear();
  h2vec.clear();

  TFile *file = TFile::Open(filename);

  // to get the names of the conditions
  TObjArray *conditions_from_name = filename.Tokenize("_"); 
  TString sa = conditions_from_name->At(1)->GetName();
  TString sb = conditions_from_name->At(3)->GetName();
  sb.ReplaceAll(".root","");

  file->cd();
  
  TDirectory *stardir = gDirectory;
  TObject *thesourcedir;
  TIter nextdir(stardir->GetListOfKeys());

  while((thesourcedir=nextdir())){
    
    TString dirName = thesourcedir->GetName();
    
    stardir->cd(dirName);
    TDirectory *current_sourcedir = gDirectory;
    TH1::AddDirectory(kFALSE);

    std::cout << "*************************" <<std::endl;
    std::cout << "Reading Directory: " << dirName <<std::endl;
    
    TObject *obj;
    TIter next(current_sourcedir->GetListOfKeys());

    while ((obj=next())) {
      
      TString objName =obj->GetName();

      if(objName.Contains("pfx")) continue;

      if (objName.Contains("h2") && !objName.Contains("pfx")) {

	//std::cout << "Reading: " << obj->GetName() <<std::endl;

	TH2F* h2 = (TH2F*)file->Get(dirName+"/"+objName);
	h2->SetName(objName+dirName);
	h2vec.push_back(h2);

	TCanvas *c = new TCanvas(h2->GetName(),h2->GetName(),800,600);
	c->cd();
	gPad->SetTopMargin(0.08);
	gPad->SetRightMargin(0.15);
	h2->SetStats(kFALSE);
	h2->Draw("colz");
	c->Draw();
	
	c->cd();
	TProfile *hpfx_tmp = (TProfile*) h2->ProfileX("_pfx",1,-1,"o");
	hpfx_tmp->SetStats(kFALSE);
	//hpfx_tmp->SetMarkerColor(kBlack);
	hpfx_tmp->SetMarkerColor(kRed);
	hpfx_tmp->SetMarkerSize(1.2); 
	hpfx_tmp->SetMarkerStyle(20); 
	hpfx_tmp->Draw("psame");
	
	c->Draw();
	cmsPrel(60.,sa,sb, isMC);
	
	TString canvName = h2->GetName();
	c->cd()->SetLogz();
	c->SaveAs(canvName+"."+format);
	  
      } else { 
	
	TH1F* h1 = (TH1F*)file->Get(dirName+"/"+objName);
	h1->SetName(objName+dirName);
	h1vec.push_back(h1);
	
	TCanvas *c = new TCanvas(h1->GetName(),h1->GetName(),600,600);
	c->cd()->SetLogy();
	
	h1->SetMarkerColor(kBlack);
	h1->SetMarkerStyle(20);
	h1->SetLineWidth(1.5); 
	h1->SetFillColor(393);
	//h1->SetFillStyle(3005);
	h1->Draw("hist");
	h1->Draw("e1same");
	c->Draw();
	
	TObject    *statObj;
	TPaveStats *stats;
  
	statObj = h1->GetListOfFunctions()->FindObject("stats");
	stats= static_cast<TPaveStats*>(statObj);
	stats->SetFillColor(10);
	stats->SetLineWidth(1);
	stats->SetShadowColor(0);
	stats->SetTextFont(42);
	stats->SetTextSize(0.025);
	//stats->SetLineColor(LineColors);
	//stats->SetTextColor(LineColors);
	stats->SetX1NDC(0.75);
	stats->SetY1NDC(0.72);
	stats->SetX2NDC(0.97);
	stats->SetY2NDC(0.92);
	stats->Draw("same"); 
		
	cmsPrel(60.,sa,sb,isMC);
	TString canvName = h1->GetName();
	c->SaveAs(canvName+"."+format);
       	
      }
    }
  }
}

