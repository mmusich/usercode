#include "TH1F.h"
#include "TFile.h"
#include "TObject.h"
#include "TKey.h"
#include "TDirectory.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>

void ExtractBTagEfficiencies(std::string filename1_="TH1F_SSVHEM.root",std::string filename2_="TH1F_SSVHPT.root"){

  ofstream outfile_;
  outfile_.open ("effbctables_cff.py");
  outfile_.precision(3);
  
  TFile *f_SSVHEM = new TFile(filename1_.c_str()); 
  TFile *f_SSVHPT = new TFile(filename2_.c_str());  

  outfile_<<"import FWCore.ParameterSet.Config as cms"<<std::endl;
  
  std::vector<TH1*> HistVect_SSVHEM;
  std::vector<TH1*> HistVect_SSVHPT;
  TString effstrings[4]={"effb_barrel=cms.untracked.double(","effb_forward=cms.untracked.double(","effc_barrel=cms.untracked.double(","effc_forward=cms.untracked.double("};

  f_SSVHEM->cd();
  HistVect_SSVHEM.push_back((TH1F*)f_SSVHEM->Get("h_eff_bTagOverGoodJet_ptb1_brl"));
  HistVect_SSVHEM.push_back((TH1F*)f_SSVHEM->Get("h_eff_bTagOverGoodJet_ptb1_fwd"));
  HistVect_SSVHEM.push_back((TH1F*)f_SSVHEM->Get("h_eff_bTagOverGoodJet_ptc1_brl"));
  HistVect_SSVHEM.push_back((TH1F*)f_SSVHEM->Get("h_eff_bTagOverGoodJet_ptc1_brl"));     

  f_SSVHPT->cd();
  HistVect_SSVHPT.push_back((TH1F*)f_SSVHPT->Get("h_eff_bTagOverGoodJet_ptb1_brl"));	 
  HistVect_SSVHPT.push_back((TH1F*)f_SSVHPT->Get("h_eff_bTagOverGoodJet_ptb1_fwd"));	 
  HistVect_SSVHPT.push_back((TH1F*)f_SSVHPT->Get("h_eff_bTagOverGoodJet_ptc1_brl"));	 
  HistVect_SSVHPT.push_back((TH1F*)f_SSVHPT->Get("h_eff_bTagOverGoodJet_ptc1_brl"));    
   
  //-------------------- SSVHEM constants ----------------------//

  outfile_<< "## MC EFF b/c SSVHEM"<<std::endl;
  outfile_<<"mc_effbc_ssvhem_2011 = cms.VPSet()"<<std::endl;
  outfile_<<"mc_effbc_ssvhem_2011.extend(["<<std::endl;
 
  Int_t maxNbin = HistVect_SSVHEM[0]->GetNbinsX();
  for(Int_t nbin=3; nbin<=maxNbin; ++nbin){
    outfile_<<"cms.PSet(ptrange=cms.untracked.vdouble("<<HistVect_SSVHEM[0]->GetBinLowEdge(nbin) <<","<<HistVect_SSVHEM[0]->GetBinLowEdge(nbin)+HistVect_SSVHEM[0]->GetBinWidth(nbin)<<"),";     
    for(UInt_t i=0; i<HistVect_SSVHEM.size(); i++){
      outfile_<<effstrings[i]<<HistVect_SSVHEM[i]->GetBinContent(nbin);
      TString opt1 = (i==HistVect_SSVHEM.size()-1) ? ")" : ")," ;
      outfile_<<opt1;
    }//loop on histograms 
    TString opt2 = (nbin==maxNbin) ? ")" : ")," ;
    outfile_<<opt2<<std::endl;
  }// loop on pt bins
  
  outfile_<<"])"<<std::endl;
  outfile_<<std::endl;
  outfile_<<std::endl;

  //-------------------- SSVHPT constants ----------------------//
  
  outfile_<<"## MC EFF b/c SSVHPT"<<std::endl;
  outfile_<<"mc_effbc_ssvhpt_2011 = cms.VPSet()"<<std::endl;
  outfile_<<"mc_effbc_ssvhpt_2011.extend(["<<std::endl;
  
   for(Int_t nbin=3; nbin<=maxNbin; ++nbin){
     outfile_<<"cms.PSet(ptrange=cms.untracked.vdouble("<<HistVect_SSVHPT[0]->GetBinLowEdge(nbin) <<","<<HistVect_SSVHPT[0]->GetBinLowEdge(nbin)+HistVect_SSVHPT[0]->GetBinWidth(nbin)<<"),";     
     for(UInt_t i=0; i<HistVect_SSVHPT.size(); i++){
       outfile_<<effstrings[i]<<HistVect_SSVHPT[i]->GetBinContent(nbin);
      TString opt1 = (i==HistVect_SSVHPT.size()-1) ? ")" : ")," ;
      outfile_<<opt1;
     }//loop on histograms 
     TString opt2 = (nbin==maxNbin) ? ")" : ")," ;
     outfile_<<opt2<<std::endl;
  }// loop on pt bins
   
  outfile_<<"])"<<std::endl;

  outfile_.close();

}

