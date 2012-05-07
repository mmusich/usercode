#include "ZbbAnalysis/AnalysisStep/interface/EventCategory.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbHistos.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbUtils.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbTypeDefs.h"
#include "ZbbAnalysis/AnalysisStep/interface/Zbbstruct4JEC.h"
#include "ZbbAnalysis/AnalysisStep/interface/VtxAssociatorsUtils.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometryVector/interface/Phi.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include "TCut.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TProfile.h"
#include "TString.h"
#include "TLorentzVector.h"
#include <Math/VectorUtil.h>

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include <iostream>
#include <algorithm>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ZbbCandidate methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ZbbCandidate::book(){
  edm::Service<TFileService> fs;
  std::string DirName = &(*N_);
  TFileDirectory HistoVars = fs->mkdir(DirName);

  h_eventWeight_ = HistoVars.make<TH1F>(N_+"_eventWeight","event weight; w_{event}",100,0.5,1.5);
  
  TH1F::SetDefaultSumw2(kTRUE);

  // lepton specific
  h_zllmass_         = HistoVars.make<TH1F>(N_+"_zmass","l^{+}l^{-} mass; M(l^{+}l^{-}) (GeV)",30,60,120);
  h_zllpt_           = HistoVars.make<TH1F>(N_+"_zpt","l^{+}l^{-} p_{T}; p_{T}(l^{+}l^{-}) (GeV)",50,0,150);
  h_zlleta_          = HistoVars.make<TH1F>(N_+"_zeta","l^{+}l^{-} #eta; #eta(l^{+}l^{-}) (GeV)",60,-3,3);
  h_lldeltaPhi_      = HistoVars.make<TH1F>(N_+"_lldeltaPhi","#Delta #phi of di-lepton ; #Delta#phi(l^{+}l^{-}) (rad)",60,0.,TMath::Pi());  
  h_lldeltaR_        = HistoVars.make<TH1F>(N_+"_lldeltaR","#Delta R of di-lepton; #DeltaR(l^{+}l^{-})",100,0,6); 
  h_zllCosThetaCS_   = HistoVars.make<TH1F>(N_+"_llCosThetaCS","Cosine of Collins-Soper frame angle; cos(#theta_{CS,l^{+}})",20,-1.,1.); 
  h_zllCosThetaStar_ = HistoVars.make<TH1F>(N_+"_llCosThetaStar","Cosine of Z decay frame l^{+} angle; cos(#theta*_{l^{+}})",20,-1.,1.); 
  h_zllAcop_         = HistoVars.make<TH1F>(N_+"_llAcop","di-lepton acoplanarity; acoplanarity #Delta#phi (rad)",50,0.,TMath::Pi());
  h_zllMt_           = HistoVars.make<TH1F>(N_+"_llMt","Transverse mass; transverse mass M_{T}(l^{+}l^{-}) (GeV)",150,0.,300.);

  // muon specific
  h_zmumumass_   = HistoVars.make<TH1F>(N_+"_zmassMu","#mu^{+}#mu^{-} mass; M(#mu^{+}#mu^{-}) (GeV)",30,60,120);
  h_zmumupt_     = HistoVars.make<TH1F>(N_+"_zptMu","#mu^{+}#mu^{-} p_{T}; p_{T}(#mu^{+}#mu^{-}) (GeV)",50,0,150);
  h_zmumueta_    = HistoVars.make<TH1F>(N_+"_zetaMu","#mu^{+}#mu^{-} #eta; #eta(#mu^{+}#mu^{-}) (GeV)",60,-3,3);
  h_mu1pt_       = HistoVars.make<TH1F>(N_+"_mu1pt","leading muon p_{T}; p^{lead #mu}_{T} (GeV)",50,0,150);
  h_mu2pt_       = HistoVars.make<TH1F>(N_+"_mu2pt","subleading muon p_{T}; p^{2nd #mu}_{T} (GeV)",50,0,150);
  h_mu1eta_      = HistoVars.make<TH1F>(N_+"_mu1eta","leading muon #eta; #eta^{lead #mu}",25,-2.5,2.5);
  h_mu2eta_      = HistoVars.make<TH1F>(N_+"_mu2eta","subleading muon #eta; #eta^{2nd #mu}",25,-2.5,2.5);
  h_mumudeltaPhi_= HistoVars.make<TH1F>(N_+"_mumudeltaPhi","#Delta #phi of di-muon ; #Delta#phi(#mu^{+}#mu^{-}) (rad)",60,0.,TMath::Pi());  
  h_mumudeltaR_  = HistoVars.make<TH1F>(N_+"_mumudeltaR","#Delta R of di-muon; #DeltaR(#mu^{+}#mu^{-})",100,0,6); 
  h_zmumuCosThetaCS_   = HistoVars.make<TH1F>(N_+"_mumuCosThetaCS","Cosine of Collins-Soper frame angle; cos(#theta_{CS,#mu^{+}})",20,-1.,1.); 
  h_zmumuCosThetaStar_ = HistoVars.make<TH1F>(N_+"_mumuCosThetaStar","Cosine of Z decay frame #mu^{+} angle; cos(#theta*_{#mu^{+}})",20,-1.,1.); 
  h_zmumuAcop_         = HistoVars.make<TH1F>(N_+"_mumuAcop","di-muon acoplanarity; acoplanarity #Delta#phi (rad)",50,0.,TMath::Pi());
  h_zmumuMt_           = HistoVars.make<TH1F>(N_+"_mumuMt","Transverse mass; transverse mass M_{T}(#mu^{+}#mu^{-}) (GeV)",150,0.,300.);

  // electron specific
  h_zeemass_     = HistoVars.make<TH1F>(N_+"_zmassEle","e^{+}e^{-} mass; M(e^{+}e^{-}) (GeV)",30,60,120);
  h_zeept_       = HistoVars.make<TH1F>(N_+"_zptEle","e^{+}e^{-} p_{T}; p_{T}(e^{+}e^{-}) (GeV)",50,0,150);
  h_zeeeta_      = HistoVars.make<TH1F>(N_+"_zetaEle","e^{+}e^{-} #eta; #eta(e^{+}e^{-}) (GeV)",60,-3,3);
  h_ele1pt_      = HistoVars.make<TH1F>(N_+"_ele1pt","leading electron p_{T};  p^{lead e}_{T} (GeV)",50,0,150);
  h_ele2pt_      = HistoVars.make<TH1F>(N_+"_ele2pt","subleading electron p_{T};  p^{2nd e}_{T} (GeV)",50,0,150);
  h_ele1eta_     = HistoVars.make<TH1F>(N_+"_ele1eta","leading electron #eta; #eta^{lead e}",25,-2.5,2.5);
  h_ele2eta_     = HistoVars.make<TH1F>(N_+"_ele2eta","subleading electron #eta; #eta^{2nd e} ",25,-2.5,2.5);
  h_eedeltaPhi_  = HistoVars.make<TH1F>(N_+"_eedeltaPhi","#Delta #phi of di-electron; #Delta#phi(e^{+}e^{-}) (rad)",60,0.,TMath::Pi());  
  h_eedeltaR_    = HistoVars.make<TH1F>(N_+"_eedeltaR","#Delta R of di-electron; #DeltaR(e^{+}e^{-})",100,0,6); 
  h_zeeCosThetaCS_   = HistoVars.make<TH1F>(N_+"_eeCosThetaCS","Cosine of Collins-Soper frame angle; cos(#theta_{CS,e^{+}})",20,-1.,1.); 
  h_zeeCosThetaStar_ = HistoVars.make<TH1F>(N_+"_eeCosThetaStar","Cosine of Z decay frame e^{+} angle; cos(#theta*_{e^{+}})",20,-1.,1.); 
  h_zeeAcop_         = HistoVars.make<TH1F>(N_+"_eeAcop","di-electron acoplanarity; acoplanarity #Delta#phi (rad)",50,0.,TMath::Pi());
  h_zeeMt_           = HistoVars.make<TH1F>(N_+"_eeMt","Transverse mass; transverse mass M_{T}(e^{+}e^{-}) (GeV)",150,0.,300.);

  // jet specific
  h_nj_            = HistoVars.make<TH1F>(N_+"_nj","jet count; N_{jets}",15,-0.5,14.5);
  h_SSVHEdisc_     = HistoVars.make<TH1F>(N_+"_SSVHEdisc","SSVHEdisc; SSV HE discriminant",100,-2.,10.);
  h_SSVHPdisc_     = HistoVars.make<TH1F>(N_+"_SSVHPdisc","SSVHPdisc; SSV HP discriminant",100,-2.,10.);
  h_TCHEdisc_      = HistoVars.make<TH1F>(N_+"_TCHEdisc","TCHEdisc; TC HE discriminant",100,-10,10);
  h_TCHPdisc_      = HistoVars.make<TH1F>(N_+"_TCHPdisc","TCHPdisc; TC HP discriminant",100,-10,10);
  h_CSVdisc_       = HistoVars.make<TH1F>(N_+"_CSVdisc","CSVdisc; CSV discriminant",48,-0.1,1.1);
  h_JPdisc_        = HistoVars.make<TH1F>(N_+"_JPdisc","JPdisc; JP discriminant",48,0.1,1.1);
  h_met_           = HistoVars.make<TH1F>(N_+"_MET","MET; #slash{E}_{T} (GeV)",100,0,200);
  h_phimet_        = HistoVars.make<TH1F>(N_+"_METphi","MET #phi; MET #phi (rad)",60,-TMath::Pi(),TMath::Pi());
  h_jetpt_         = HistoVars.make<TH1F>(N_+"_jetpt","All jets p_{T}; Jet p_{T} (GeV)",50,15,215);
  h_jeteta_        = HistoVars.make<TH1F>(N_+"_jeteta","All jets #eta; Jet #eta",25,-2.5, 2.5);
  h_jetphi_        = HistoVars.make<TH1F>(N_+"_jetphi","All jet #phi; Jet #phi (rad)",60,-TMath::Pi(),TMath::Pi());
  h_jetoverlapLept_= HistoVars.make<TH1F>(N_+"_jetoverlapLept","jets overlaps with leptons from Z",2,-0.5,1.5);

  TString jetOverlapLeptBinLabels[2] ={"no overlap","has overlap"};

  for(UInt_t bin=1; bin<=2; bin++){
    h_jetoverlapLept_->GetXaxis()->SetBinLabel(bin,jetOverlapLeptBinLabels[bin-1]);    
  }

  h_jet1pt_        = HistoVars.make<TH1F>(N_+"_jet1pt","leading jet p_{T}; leading jet p_{T} (GeV)",50,15,215);
  h_jet1eta_       = HistoVars.make<TH1F>(N_+"_jet1eta","leading jet #eta; leading jet #eta",25,-2.5,2.5);
  h_SSVHEdisc1_    = HistoVars.make<TH1F>(N_+"_SSVHEdisc1","SSVHEdisc 1^{st} b-jet; SSV HE discriminant 1^{st} b-jet",100,0.,10);
  h_SSVHPdisc1_    = HistoVars.make<TH1F>(N_+"_SSVHPdisc1","SSVHPdisc 1^{st} b-jet; SSV HP discriminant 1^{st} b-jet",100,0.,10);
  h_TCHEdisc1_     = HistoVars.make<TH1F>(N_+"_TCHEdisc1","TCHEdisc 1^{st} b-jet; TC HE discriminant 1^{st} b-jet",100,-10,10);
  h_TCHPdisc1_     = HistoVars.make<TH1F>(N_+"_TCHPdisc1","TCHPdisc 1^{st} b-jet; TC HP discriminant 1^{st} b-jet",100,-10,10);
  h_CSVdisc1_      = HistoVars.make<TH1F>(N_+"_CSVdisc1","CSVdisc 1^{st} b-jet; CSV discriminant 1^{st} b-jet",48,-0.1,1.1);
  h_JPdisc1_       = HistoVars.make<TH1F>(N_+"_JPdisc1","JPdisc 1^{st} b-jet; JP discriminant 1^{st} b-jet",48,0.,1.1);

  //ideally replace DISCR w/ bTagAlgoWP as in the cfg file
  //float bTag_pTbins[12]={0, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 150};
  //p_bTagAlgoWP_SFb_pt_    = HistoVars.make<TProfile>(N_+"_DISCR_SF b_pt",  "SF_b vs b-jet pT DISCR; pT b-jet; SF_b DISCR",11,bTag_pTbins);
  //p_bTagAlgoWP_EffbMC_pt_ = HistoVars.make<TProfile>(N_+"_DISCR_EFFbMC_pt","Eff_b MC vs b-jet pT DISCR; pT b-jet; Eff_b MC DISCR",11,bTag_pTbins);

  h_bTagAlgoWP_MCflav_ = HistoVars.make<TH1F>(N_+"_DISCR_MCflav","DISCR MC flavour; MC flavour",22,-0.5,21.5); 

  h_jet2pt_        = HistoVars.make<TH1F>(N_+"_jet2pt","subleading jet p_{T}; subleading jet p_{T} (GeV)",50,15,215);
  h_jet2eta_       = HistoVars.make<TH1F>(N_+"_jet2eta","subleading jet #eta; subleading jet #eta",25,-2.5,2.5);
  h_SSVHEdisc2_    = HistoVars.make<TH1F>(N_+"_SSVHEdisc2","SSVHEdisc; SSV HE discriminant 2^{nd} b-jet",100,0.,10);
  h_SSVHPdisc2_    = HistoVars.make<TH1F>(N_+"_SSVHPdisc2","SSVHPdisc; SSV HP discriminant 2^{nd} b-jet",100,0.,10);
  h_TCHEdisc2_     = HistoVars.make<TH1F>(N_+"_TCHEdisc2","TCHEdisc; TC HE discriminant 2^{nd} b-jet",100,-10,10);
  h_TCHPdisc2_     = HistoVars.make<TH1F>(N_+"_TCHPdisc2","TCHPdisc; TC HP discriminant 2^{nd} b-jet",100,-10,10);
  h_CSVdisc2_      = HistoVars.make<TH1F>(N_+"_CSVdisc2","CSVdisc 2^{nd} b-jet; CSV discriminant 2^{nd} b-jet",48,-0.1,1.1);
  h_JPdisc2_       = HistoVars.make<TH1F>(N_+"_JPdisc2","JPdisc 2^{st} b-jet; JP discriminant 2^{st} b-jet",48,0.,1.1);

  h_nhf_           = HistoVars.make<TH1F>(N_+"_nhf","neutral hadron energy fraction; E_{n. had}/E_{tot}",50,0,1.01);
  h_nef_           = HistoVars.make<TH1F>(N_+"_nef","neutral EmEnergy fraction; E_{n. em}/E_{tot}",50,0,1.01);
  h_nconstituents_ = HistoVars.make<TH1F>(N_+"_npf","total multiplicity; jet total multiplicity",50,-0.5,49.5);
  h_chf_           = HistoVars.make<TH1F>(N_+"_chf","charged hadron energy fraction; E_{ch. had}/E_{tot}",50,0,1.01);
  h_nch_           = HistoVars.make<TH1F>(N_+"_nch","charged multiplicity; jet charged multiplicity",10,-0.5,9.5);
  h_cef_           = HistoVars.make<TH1F>(N_+"_cef","charged EmEnergy fraction; E_{ch. em}/E_{tot}",50,0,1.01);
  h_alljetDeltaR_  = HistoVars.make<TH1F>(N_+"_alljetDeltaR","#Delta R between all jet pairs; #DeltaR(j,j)",100,0,10); 
  h_alljetDeltaPhi_= HistoVars.make<TH1F>(N_+"_alljetDeltaPhi","#Delta #phi of all jet pairs; #Delta#phi(j,j) (rad)",60,0.,TMath::Pi());  
  
  h_jetid_         = HistoVars.make<TH1F>(N_+"_jetid","Jet Id level (none, loose, medium, tight)",4,-0.5,3.5);
  
  TString jetIdBinLabels[4] ={"none","loose","medium","tight"};
   
  for(UInt_t bin=1; bin<=4; bin++){
    h_jetid_->GetXaxis()->SetBinLabel(bin,jetIdBinLabels[bin-1]);    
  }

  //b jet stuff
  h_nb_          = HistoVars.make<TH1F>(N_+"_nb","b-tagged jet count; N_{b-jets}",5,-0.5,4.5);
  h_njb_         = HistoVars.make<TH2I>(N_+"_njVsb","number of bjets vs number of jets;N_{jets};N_{b-jets}",15,0,15,5,0,5);
  h_sumbHtOversumHt_ =HistoVars.make<TH1F>(N_+"_sumbHtOversumHt","Scalar sum b-jet p_{T}/ Scalar sum of all jet p_{T};sum_{bj}H_{T}/sum_{j}H_{T}",100,0.,1.);
  h_bjetpt_      = HistoVars.make<TH1F>(N_+"_bjetpt","all b-jetd p_{T}; all b-jets p_{T} (GeV)",50,15,215);
  h_bjeteta_     = HistoVars.make<TH1F>(N_+"_bjeteta"," all b-jets #eta; all b-jets #eta",50,-2.5,2.5);
  h_bjetmass_    = HistoVars.make<TH1F>(N_+"_bjetmass","all b-jet mass; all b-jets m_{bj} (GeV)",100,0.,100.);
  h_bjetOtherJetsDeltaR_  = HistoVars.make<TH1F>(N_+"_bjetOtherJetsDeltaR","#Delta R between b-jet and all other jets; #DeltaR(j_{b},j)",100,0,10); 
  h_bjetOtherJetsDeltaPhi_= HistoVars.make<TH1F>(N_+"_bjetOtherJetsDeltaPhi","#Delta #phi of b-jet with other jets; #Delta#phi(j_{b},j) (rad)",60,0.,TMath::Pi());  

  //general SV infos

  h_SVnVertices_ = HistoVars.make<TH1F>(N_+"_SVnVertices","number of secondary vertices; N_{SV}",5,-0.5,4.5);
  h_SVJetdeltaR_ = HistoVars.make<TH1F>(N_+"_SVdeltaR","#DeltaR between vertex direction and jet direction; #DeltaR(#vec{SV},#vec{j}_{b})",100,0.,0.5);
  h_SVmass_      = HistoVars.make<TH1F>(N_+"_SVmass","vertex mass; M_{SV} (GeV)",100,0.,10.);
  h_SVdist_      = HistoVars.make<TH1F>(N_+"_SVdist","transverse distance between PV and SV; d_{xy}(PV,SV) (cm)",100,0.,2.);
  h_SVdistErr_   = HistoVars.make<TH1F>(N_+"_SVdistErr","transverse distance error between PV and SV; #sigma_{d(PV,SV)} (cm)",100,0.,0.5);
  h_SVdistSig_   = HistoVars.make<TH1F>(N_+"_SVdistSig","transverse distance significance between PV and SV; d_{xy}(PV,SV)/#sigma_{d(PV,SV)}",100, 0., 50.);
  h_SVnTracks_   = HistoVars.make<TH1F>(N_+"_SVnTracks","number of tracks at secondary vertex; N_{tracks}^{SV}",20,-0.5, 19.5);
  h_SVchi2_      = HistoVars.make<TH1F>(N_+"_SVchi2","secondary vertex fit #chi^{2}; #chi^{2}(SV)",100,0., 50.);
  
  //jet ordering specific
  h_bjet1pt_     = HistoVars.make<TH1F>(N_+"_bjet1pt","leading bjet p_{T}; leading b-jet p_{T} (GeV)",50,15,215);
  h_bjet1eta_    = HistoVars.make<TH1F>(N_+"_bjet1eta","leading bjet #eta; leading b-jet #eta",25,-2.5,2.5);
  h_bjet1mass_   = HistoVars.make<TH1F>(N_+"_bjet1mass","leading b-jet mass; leading b-jet m_{bj} (GeV)",100,0.,100.);
  h_bjet2pt_     = HistoVars.make<TH1F>(N_+"_bjet2pt","subleading bjet p_{T}; subleading b-jet p_{T} (GeV)",50,15,215);
  h_bjet2eta_    = HistoVars.make<TH1F>(N_+"_bjet2eta","subleading bjet #eta; subleading b-jet #eta",25,-2.5,2.5);
  h_bjet2mass_   = HistoVars.make<TH1F>(N_+"_bjet2mass","subleading b-jet mass; leading  b-jet m_{bj} (GeV)",100,0.,100.);
  h_scaldptZbj1_ = HistoVars.make<TH1F>(N_+"_scaldptZbj1","Scalar difference of Z and leading b-jet p_{T}; p_{T}^{Z} - p_{T}^{lead b} (GeV)",100,0,500);
  h_drZbj1_      = HistoVars.make<TH1F>(N_+"_drZbj1","#DeltaR between Z and leading b-jet; #Delta R(Z,b)",50,0,5);
  h_vecdptZbj1_  = HistoVars.make<TH1F>(N_+"_vecdptZbj1","Vector difference of Z and leading b-jet p_{T}; | #vec{p}_{T}^{Z} - #vec{p}_{T}^{lead b}| (GeV)",50,0,500);
  h_dphiZbj1_    = HistoVars.make<TH1F>(N_+"_dphiZbj1","#Delta#phi between Z and leading b-jet; #Delta #phi(Z,b) (rad)",60,0,TMath::Pi());
  h_bbM_         = HistoVars.make<TH1F>(N_+"_dijetM","b#bar{b} invariant mass; M(b#bar{b}) (GeV)",50,0,1000);
  h_bbPt_        = HistoVars.make<TH1F>(N_+"_dijetPt","b#bar{b} p_{T}; p_{T}(b#bar{b}) (GeV)",50,0,500);
  h_bbDeltaR_    = HistoVars.make<TH1F>(N_+"_bbDeltaR","#Delta R between b-jets for b#bar{b} events; #DeltaR(j_{b},j_{#bar{b}})",100,0,10); 
  h_ZbM_         = HistoVars.make<TH1F>(N_+"_ZbM","Z+b invariant mass; M(Z+b) (GeV)",50,0,1000);
  h_ZbPt_        = HistoVars.make<TH1F>(N_+"_ZbPt","Z+b p_{T}; p_{T}(Z+b) (GeV)",50,0,500);
  h_ZbbM_        = HistoVars.make<TH1F>(N_+"_ZbbM","Zb#bar{b} invariant mass; M(Z+b#bar{b}) (GeV)",50,0,1000);
  h_ZbbPt_       = HistoVars.make<TH1F>(N_+"_ZbbPt","Zb#bar{b} p_{T}; p_{T}(Z+b#bar{b}) (GeV)",50,0,500);
  h_ZbbM2D_      = HistoVars.make<TH2F>(N_+"_ZbbM2D","Zb#bar{b} mass vs b#bar{b} mass; M(b#bar{b}) (GeV); M(Z+b#bar{b}) (GeV)",50,0,1000,50,0,1000);

}

void ZbbCandidate::fillMu(const EventCategory& ec_mu){

  //std::cout<<"ZbbCandidate::fillMu"<<std::endl;
  Bool_t thecut_= ZbbUtils::ProduceTheCut(ec_mu,N_,"");
  Double_t w = ZbbUtils::getTheWeight(ec_mu,N_);

  if(thecut_){
    h_zllmass_->Fill(ec_mu.bestZcandidate_.mass(),w);  
    h_zllpt_->Fill(ec_mu.bestZcandidate_.pt(),w);  
    h_zlleta_->Fill(ec_mu.bestZcandidate_.eta(),w);
    h_zmumueta_->Fill(ec_mu.bestZcandidate_.eta(),w);
    h_zmumumass_->Fill(ec_mu.bestZcandidate_.mass(),w);  
    h_zmumupt_->Fill(ec_mu.bestZcandidate_.pt(),w);  
    h_mu1pt_->Fill(ec_mu.bestZcandidate_.daughter(0)->pt(),w);    
    h_mu2pt_->Fill(ec_mu.bestZcandidate_.daughter(1)->pt(),w);    
    h_mu1eta_->Fill(ec_mu.bestZcandidate_.daughter(0)->eta(),w);   
    h_mu2eta_->Fill(ec_mu.bestZcandidate_.daughter(1)->eta(),w);
    double deltaRll   = ROOT::Math::VectorUtil::DeltaR(ec_mu.bestZcandidate_.daughter(0)->momentum(),ec_mu.bestZcandidate_.daughter(1)->momentum());
    double deltaPhill = ZbbUtils::myDeltaPhi(ec_mu.bestZcandidate_.daughter(0)->phi(),ec_mu.bestZcandidate_.daughter(1)->phi());
    h_lldeltaR_->Fill(deltaRll,w);
    h_lldeltaPhi_->Fill(deltaPhill,w); 
    h_mumudeltaR_->Fill(deltaRll,w);  
    h_mumudeltaPhi_->Fill(deltaPhill,w);   
    h_zllCosThetaCS_->Fill(ZbbUtils::CostCS(ec_mu.bestZcandidate_),w);   
    h_zllCosThetaStar_->Fill(ZbbUtils::CosThetaStar(ec_mu.bestZcandidate_),w); 
    h_zmumuCosThetaCS_->Fill(ZbbUtils::CostCS(ec_mu.bestZcandidate_),w);      
    h_zmumuCosThetaStar_->Fill(ZbbUtils::CosThetaStar(ec_mu.bestZcandidate_),w);

    Geom::Phi<double> deltaphi_(ec_mu.bestZcandidate_.daughter(0)->phi()-atan2(ec_mu.bestZcandidate_.daughter(1)->py(),ec_mu.bestZcandidate_.daughter(1)->px()));
    double acop = deltaphi_.value();
    if (acop<0) acop = - acop;
    acop = TMath::Pi() - acop;

    h_zmumuAcop_->Fill(acop);      
    h_zmumuMt_->Fill(ec_mu.bestZcandidate_.mt(),w);        
    h_zllAcop_->Fill(acop);      
    h_zllMt_->Fill(ec_mu.bestZcandidate_.mt(),w);    

  }
  return;  
}

void ZbbCandidate::fillEle(const EventCategory& ec_ele){

  //std::cout<<"ZbbCandidate::fillEle"<<std::endl;
  Bool_t thecut_ = ZbbUtils::ProduceTheCut(ec_ele,N_,"");
  Double_t w = ZbbUtils::getTheWeight(ec_ele,N_);

  if(thecut_){
    h_zllmass_->Fill(ec_ele.bestZcandidate_.mass(),w);  
    h_zllpt_->Fill(ec_ele.bestZcandidate_.pt(),w);  
    h_zlleta_->Fill(ec_ele.bestZcandidate_.eta(),w);
    h_zeeeta_->Fill(ec_ele.bestZcandidate_.eta(),w);
    h_zeemass_->Fill(ec_ele.bestZcandidate_.mass(),w);  
    h_zeept_->Fill(ec_ele.bestZcandidate_.pt(),w);  
    h_ele1pt_->Fill(ec_ele.bestZcandidate_.daughter(0)->pt(),w);    
    h_ele2pt_->Fill(ec_ele.bestZcandidate_.daughter(1)->pt(),w);    
    h_ele1eta_->Fill(ec_ele.bestZcandidate_.daughter(0)->eta(),w);   
    h_ele2eta_->Fill(ec_ele.bestZcandidate_.daughter(1)->eta(),w);  
    double deltaRll   = ROOT::Math::VectorUtil::DeltaR(ec_ele.bestZcandidate_.daughter(0)->momentum(),ec_ele.bestZcandidate_.daughter(1)->momentum());
    double deltaPhill = ZbbUtils::myDeltaPhi(ec_ele.bestZcandidate_.daughter(0)->phi(),ec_ele.bestZcandidate_.daughter(1)->phi());
    h_lldeltaR_->Fill(deltaRll,w);
    h_lldeltaPhi_->Fill(deltaPhill,w); 
    h_eedeltaR_->Fill(deltaRll,w);  
    h_eedeltaPhi_->Fill(deltaPhill,w);   
    h_zllCosThetaCS_->Fill(ZbbUtils::CostCS(ec_ele.bestZcandidate_),w);   
    h_zllCosThetaStar_->Fill(ZbbUtils::CosThetaStar(ec_ele.bestZcandidate_),w); 
    h_zeeCosThetaCS_->Fill(ZbbUtils::CostCS(ec_ele.bestZcandidate_),w);      
    h_zeeCosThetaStar_->Fill(ZbbUtils::CosThetaStar(ec_ele.bestZcandidate_),w); 

    Geom::Phi<double> deltaphi_(ec_ele.bestZcandidate_.daughter(0)->phi()-atan2(ec_ele.bestZcandidate_.daughter(1)->py(),ec_ele.bestZcandidate_.daughter(1)->px()));
    double acop = deltaphi_.value();
    if (acop<0) acop = - acop;
    acop = TMath::Pi() - acop;

    h_zeeAcop_->Fill(acop);      
    h_zeeMt_->Fill(ec_ele.bestZcandidate_.mt(),w);        
    h_zllAcop_->Fill(acop);      
    h_zllMt_->Fill(ec_ele.bestZcandidate_.mt(),w);  

  }
  return;
}

void ZbbCandidate::fillJets(const EventCategory& ec_mu,const EventCategory& ec_ele, edm::Handle<edm::View<pat::Jet> > jets, edm::Handle<edm::View<pat::MET> > mets, std::string bTagAlgoWP,std::map<EffMCpair,std::vector<double> > theAssociativeMapEffbMC, std::map<EffMCpair,std::vector<double> > theAssociativeMapEffcMC){

  //std::cout<<"ZbbCandidate::fillJets"<<std::endl;
  int nJets = 0;
  int nBjets = 0;
  
  Bool_t thecut_;   
  EventCategory ec; 
  
  if(ec_mu.isZLL()){
    thecut_= ZbbUtils::ProduceTheCut(ec_mu,N_,"");
    ec = ec_mu;
  } else if(ec_ele.isZLL()) {
    thecut_ = ZbbUtils::ProduceTheCut(ec_ele,N_,"");
    ec = ec_ele;
  } else {
    return;
  }
  
  Double_t w = ZbbUtils::getTheWeight(ec,N_);

  // vector of b-tagged jets
  std::vector<pat::Jet> theBtaggedJets;
  // vector of good Jets
  std::vector<pat::Jet> theGoodJets;

  if(thecut_){

    // fill the event-weight after selection
    h_eventWeight_->Fill(w);

    for(edm::View<pat::Jet>::const_iterator jetcand=jets->begin(); jetcand!=jets->end(); ++jetcand){
       
      if( ZbbUtils::hasOverlap(jetcand->eta(),jetcand->phi(),ec.bestZcandidate_.daughter(0)->eta(),ec.bestZcandidate_.daughter(0)->phi(),0.5) && 
	  ZbbUtils::hasOverlap(jetcand->eta(),jetcand->phi(),ec.bestZcandidate_.daughter(1)->eta(),ec.bestZcandidate_.daughter(1)->phi(),0.5) ){
	h_jetoverlapLept_->Fill(1.,w);
      }  
      else {
	h_jetoverlapLept_->Fill(0.,w);	
      }

      if ( ZbbUtils::isJetIdOk((*jetcand),"tight")){
	h_jetid_->Fill(3.,w);
      } else if (ZbbUtils::isJetIdOk((*jetcand),"medium")){
	h_jetid_->Fill(2.,w);
      } else if (ZbbUtils::isJetIdOk((*jetcand),"loose")){
	h_jetid_->Fill(1.,w);
      } else{ 
	h_jetid_->Fill(0.,w);
      }

      if(ZbbUtils::isJetIdOk((*jetcand),"loose") && ZbbUtils::isGoodJet((*jetcand),ec.bestZcandidate_,lCuts_)){
	nJets++;
	
	theGoodJets.push_back(*jetcand);

	h_jetpt_->Fill(jetcand->pt(),w);         
	h_jeteta_->Fill(jetcand->eta(),w);        
	h_jetphi_->Fill(jetcand->phi(),w);       
	h_nhf_->Fill((jetcand->neutralHadronEnergy() + jetcand->HFHadronEnergy() ) / jetcand->energy(),w);           
	h_nef_->Fill(jetcand->neutralEmEnergyFraction(),w);           
	h_nconstituents_->Fill(jetcand->numberOfDaughters(),w); 
	h_chf_->Fill(jetcand->chargedHadronEnergyFraction(),w);          
	h_nch_->Fill(jetcand->chargedHadronEnergyFraction(),w);           
	h_cef_->Fill(jetcand->chargedEmEnergyFraction(),w);           
 
	h_SSVHEdisc_->Fill(jetcand->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"),w);     
	h_SSVHPdisc_->Fill(jetcand->bDiscriminator("simpleSecondaryVertexHighPurBJetTags"),w);     
	h_TCHEdisc_->Fill(jetcand->bDiscriminator("trackCountingHighEffBJetTags"),w);    
        h_TCHPdisc_->Fill(jetcand->bDiscriminator("trackCountingHighPurBJetTags"),w);  
	h_CSVdisc_->Fill(jetcand->bDiscriminator("combinedSecondaryVertexBJetTags"),w);
	h_JPdisc_->Fill(jetcand->bDiscriminator("jetProbabilityBJetTags"),w);

	if(ZbbUtils::isBJet((*jetcand),bTagAlgoWP)){
	  theBtaggedJets.push_back(*jetcand);
	  nBjets++;
	  if ( ec.isMC() ) {
	    h_bTagAlgoWP_MCflav_->Fill(abs(jetcand->partonFlavour()));
	    // ZbbUtils::FLAVOUR theJetflavour = ZbbUtils::getPartonFlavour(*jetcand);
	    // EM 2011.07.15 //if (theJetflavour==ZbbUtils::B_JETS) p_bTagAlgoWP_SFb_pt_->Fill(std::min(jetcand->pt(),149.0),ZbbUtils::getbEffScaleFactor(theJetflavour,*jetcand,iSetup)); 
	  } 
	}     	
      }
    } // close the for on jets 

    h_met_->Fill((*mets)[0].et(),w);           
    h_phimet_->Fill((*mets)[0].phi(),w);  
    h_nj_->Fill(nJets,w); 
    h_nb_->Fill(nBjets,w); 
    h_njb_->Fill(nJets,nBjets,w);
    
    //if good jets

    if(nJets>0){
      h_jet1pt_->Fill(theGoodJets[0].pt(),w);       
      h_jet1eta_->Fill(theGoodJets[0].eta(),w);      
      if(nJets>1){
	h_jet2pt_->Fill(theGoodJets[1].pt(),w);       
	h_jet2eta_->Fill(theGoodJets[1].eta(),w);
      }
    }
	
    //if b-tagged jets

    Double_t allJetsSumPt_(0);
    Double_t allbJetsSumPt_(0);

    for(UInt_t i=0; i<theGoodJets.size();i++){
      allJetsSumPt_+=theGoodJets[i].pt();
    }

    if(nBjets>0){

      for(UInt_t i=0; i<theBtaggedJets.size();i++){
	h_bjetpt_->Fill(theBtaggedJets[i].pt(),w); 
	h_bjeteta_->Fill(theBtaggedJets[i].eta(),w);
	h_bjetmass_->Fill(theBtaggedJets[i].mass(),w);
	allbJetsSumPt_+=theBtaggedJets[i].pt();
      }

      h_sumbHtOversumHt_->Fill((allbJetsSumPt_/allJetsSumPt_),w);

      h_bjet1pt_->Fill(theBtaggedJets[0].pt(),w); 
      h_bjet1eta_->Fill(theBtaggedJets[0].eta(),w);
      h_bjet1mass_->Fill(theBtaggedJets[0].mass(),w);

      h_SSVHEdisc1_->Fill(theBtaggedJets[0].bDiscriminator("simpleSecondaryVertexHighEffBJetTags"),w);     
      h_SSVHPdisc1_->Fill(theBtaggedJets[0].bDiscriminator("simpleSecondaryVertexHighPurBJetTags"),w);     
      h_TCHEdisc1_->Fill(theBtaggedJets[0].bDiscriminator("trackCountingHighEffBJetTags"),w);    
      h_TCHPdisc1_->Fill(theBtaggedJets[0].bDiscriminator("trackCountingHighPurBJetTags"),w);  
      h_CSVdisc1_->Fill(theBtaggedJets[0].bDiscriminator("combinedSecondaryVertexBJetTags"),w);
      h_JPdisc1_->Fill(theBtaggedJets[0].bDiscriminator("jetProbabilityBJetTags"),w);

      TLorentzVector b1P4_(theBtaggedJets[0].px(),theBtaggedJets[0].py(),theBtaggedJets[0].pz(),theBtaggedJets[0].energy());
      TLorentzVector ZP4_(ec.bestZcandidate_.px(),ec.bestZcandidate_.py(),ec.bestZcandidate_.pz(),ec.bestZcandidate_.energy());
      TLorentzVector ZbP4_ = ZP4_+b1P4_;

      h_scaldptZbj1_->Fill(ZbP4_.Pt(),w);     
      h_vecdptZbj1_->Fill(ec.bestZcandidate_.pt() - theBtaggedJets[0].pt(),w);  
      h_ZbM_->Fill(ZbP4_.M(),w); 
      h_ZbPt_->Fill(ZbP4_.Pt(),w);
      h_dphiZbj1_->Fill(deltaPhi(ec.bestZcandidate_.phi(),theBtaggedJets[0].phi()),w);    
      h_drZbj1_->Fill(deltaR(ec.bestZcandidate_.eta(),ec.bestZcandidate_.phi(),theBtaggedJets[0].eta(),theBtaggedJets[0].phi()),w);  

      if(nBjets>1){
	h_bjet2pt_->Fill(theBtaggedJets[1].pt(),w);        
	h_bjet2eta_->Fill(theBtaggedJets[1].eta(),w);
	h_bjet2mass_->Fill(theBtaggedJets[1].mass(),w);
	
	h_SSVHEdisc2_->Fill(theBtaggedJets[1].bDiscriminator("simpleSecondaryVertexHighEffBJetTags"),w);     
	h_SSVHPdisc2_->Fill(theBtaggedJets[1].bDiscriminator("simpleSecondaryVertexHighPurBJetTags"),w);     
	h_TCHEdisc2_->Fill(theBtaggedJets[1].bDiscriminator("trackCountingHighEffBJetTags"),w);    
	h_TCHPdisc2_->Fill(theBtaggedJets[1].bDiscriminator("trackCountingHighPurBJetTags"),w); 
	h_CSVdisc2_->Fill(theBtaggedJets[1].bDiscriminator("combinedSecondaryVertexBJetTags"),w);
	h_JPdisc2_->Fill(theBtaggedJets[1].bDiscriminator("jetProbabilityBJetTags"),w);

	TLorentzVector b2P4_(theBtaggedJets[1].px(),theBtaggedJets[1].py(),theBtaggedJets[1].pz(),theBtaggedJets[1].energy());
	TLorentzVector bbP4_ = b1P4_ + b2P4_;
        h_bbM_->Fill(bbP4_.M(),w);
	h_bbPt_->Fill(bbP4_.Pt(),w);
	h_bbDeltaR_->Fill(ROOT::Math::VectorUtil::DeltaR(theBtaggedJets[0].momentum(),theBtaggedJets[1].momentum()),w);
       	TLorentzVector ZbbP4_ = ZbP4_ + b2P4_;
	h_ZbbM_->Fill(ZbbP4_.M(),w);
	h_ZbbPt_->Fill(ZbbP4_.Pt(),w);
	h_ZbbM2D_->Fill(ZbbP4_.M(),bbP4_.M(),w);
	
      } // close if more than 1 btag
    } // close if btag

    // loop for SV informations
    for(std::vector<pat::Jet>::const_iterator bjet=theBtaggedJets.begin(); bjet!=theBtaggedJets.end(); ++bjet){
      
      // retrieve the "secondary vertex tag infos"
      // this is the output of the b-tagging reconstruction code
      // and contains secondary vertices in the jets
      const reco::SecondaryVertexTagInfo &svTagInfo = *bjet->tagInfoSecondaryVertex("secondaryVertex");

      // count the number of secondary vertices
      h_SVnVertices_ ->Fill(svTagInfo.nVertices());
      
       // ignore jets without SV from now on
      if (svTagInfo.nVertices() < 1)
	continue;
      
      // pick the first secondary vertex (the "best" one)
      const reco::Vertex &sv = svTagInfo.secondaryVertex(0);
      
      // and plot number of tracks and chi^2
      h_SVnTracks_->Fill(sv.tracksSize());
      h_SVchi2_->Fill(sv.chi2());
      
      // the precomputed transverse distance to the primary vertex
      Measurement1D distance = svTagInfo.flightDistance(0, true);
      h_SVdist_->Fill(distance.value());
      h_SVdistErr_->Fill(distance.error());       
      h_SVdistSig_->Fill(distance.significance());
      
      // the precomputed direction with respect to the primary vertex
      GlobalVector dir = svTagInfo.flightDirection(0);
      
      // unfortunately CMSSW hsa all kinds of vectors, and sometimes we need to convert them *sigh*
      math::XYZVector dir2(dir.x(), dir.y(), dir.z());
      
      // compute a few variables that we are plotting
      double deltaR = ROOT::Math::VectorUtil::DeltaR(bjet->momentum(), dir2);
      h_SVJetdeltaR_->Fill(deltaR);
      
      // compute the invariant mass from a four-vector sum
      math::XYZTLorentzVector trackFourVectorSum;
      
      // loop over all tracks in the vertex
      for(reco::Vertex::trackRef_iterator track = sv.tracks_begin(); track != sv.tracks_end(); ++track) {
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> > vec;
	vec.SetPx((*track)->px());
	vec.SetPy((*track)->py());
	vec.SetPz((*track)->pz());
	vec.SetM(0.13957);	// pion mass
	trackFourVectorSum += vec;
      }
      
      // get the invariant mass: sqrt(E^2 - px^2 - py^2 - pz^2)
      double vertexMass = trackFourVectorSum.M();
      h_SVmass_->Fill(vertexMass);
    } //closes loop on bjets

    for(std::vector<pat::Jet>::const_iterator jet1=theGoodJets.begin(); jet1!= theGoodJets.end(); ++jet1){
      for(std::vector<pat::Jet>::const_iterator jet2=jet1+1; jet2!=theGoodJets.end(); ++jet2){
	if(jet2!=jet1){
	  double deltaRjj   = ROOT::Math::VectorUtil::DeltaR(jet1->momentum(),jet2->momentum());
	  double deltaPhijj = ZbbUtils::myDeltaPhi(jet1->phi(),jet2->phi());
	  h_alljetDeltaR_->Fill(deltaRjj,w);  
	  h_alljetDeltaPhi_->Fill(deltaPhijj,w);  
	} // check self
      } // loop 1 on jets
    } // loop 2 on jets
    
    if(nBjets>0){
      for(std::vector<pat::Jet>::const_iterator bjet=theBtaggedJets.begin(); bjet!=theBtaggedJets.end(); ++bjet){
	for(std::vector<pat::Jet>::const_iterator jet1=theGoodJets.begin(); jet1!= theGoodJets.end(); ++jet1){
	  if(jet1==bjet) 	   // check self
	    continue;
	  double deltaRbj   = ROOT::Math::VectorUtil::DeltaR(bjet->momentum(),jet1->momentum());
	  double deltaPhibj = ZbbUtils::myDeltaPhi(bjet->phi(),jet1->phi());
	  if(deltaRbj!=0 && deltaPhibj!=0){
	    h_bjetOtherJetsDeltaR_->Fill(deltaRbj,w);    
	    h_bjetOtherJetsDeltaPhi_->Fill(deltaPhibj,w);
	  }
	} // loop on b-jets 
      } // if bjet
    } //extra loop on jets

  } // if right category
  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ZbbCandidateNtuple methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ZbbCandidateNtuple::book(){
  edm::Service<TFileService> fs;
  std::string DirName = &(*N_);
  TFileDirectory HistoVars = fs->mkdir(DirName);

  //Reset the structures
  Event.reset(); Z.reset(); bjet.reset(); jet2.reset(); MET.reset();
  btagjets.reset();

  //Define the Tree and the branches
  ZbbNtuple = HistoVars.make<TTree>("ZbbNtuple","Data for JEC/SV studies");
  ZbbNtuple->Branch("Event",&Event,"run_Event/I:lumisec_Event:num_Event:njet_Event:nbjet_Event:nvert_Event:weight_Event/F:isMuon/O");
  ZbbNtuple->Branch("Z",&Z,"pT_Z/F:pX_Z:pY_Z:eta_Z:phi_Z:eta_d1:phi_d1:pt_d1:eta_d2:phi_d2:pt_d2");
  ZbbNtuple->Branch("bjet",&bjet,"pT_bjet/F:eta_bjet:phi_bjet:dzDilept_bjet:frac2b_bjet:frac1_bjet:frac2_bjet:frac3_bjet");
  ZbbNtuple->Branch("jet2",&jet2,"pT_jet2/F:eta_jet2:phi_jet2:dzDilept_jet2:frac2b_jet2:frac1_jet2:frac2_jet2:frac3_jet2");
  ZbbNtuple->Branch("MET",&MET,"pT_MET/F:px_MET:py_MET");
  ZbbNtuple->Branch("btagjets",&btagjets,"trueflavour_bjets[6]/I:pT_bjets[6]/F:eta_bjets[6]:phi_bjets[6]:ssvhediscr_bjets[6]:ssvhpdiscr_bjets[6]:csvdiscr_bjets[6]:jpdiscr_bjets[6]:jpdiscr_fromtags_bjets[6]:tchediscr_bjets[6]:tchpdiscr_bjets[6]:svmass_bjets[6]");
  
}


void ZbbCandidateNtuple::fill(const EventCategory& ec_mu, const EventCategory& ec_ele, edm::Handle<edm::View<pat::Jet> > jets, std::string bTagAlgoWP, edm::Handle<edm::View<pat::MET> > mets, const edm::Event& iEvent, Int_t nvvertex){

  //Reset the structures
  Event.reset(); Z.reset(); bjet.reset(); jet2.reset(); MET.reset(); 

  Bool_t thecut_;   
  EventCategory ec; 

  // vector of good Jets
  std::vector<pat::Jet> theJets;
  std::vector<pat::Jet> theBTaggedJets;

  if(ec_mu.isZLL()){
    thecut_= ZbbUtils::ProduceTheCut(ec_mu,N_,"");
    ec = ec_mu;
    Event.isMuon = true;
  } else if(ec_ele.isZLL()) {
    thecut_ = ZbbUtils::ProduceTheCut(ec_ele,N_,"");
    ec = ec_ele;
    Event.isMuon = false;
  } else {
    return;
  }

  // returns the Z vertex from EventCategory
  // const reco::Vertex& theZVertex = Zvertexes->at(0);
  const reco::Vertex& theZVertex = ec.theZvertex_;
  
  if ( !thecut_ ) return; // event not selected according to the right category

  // Get b tag information (JP)
  edm::Handle<reco::JetTagCollection> bTagHandleJP;
  iEvent.getByLabel("MyJetProbabilityBJetTags", bTagHandleJP);
  const reco::JetTagCollection & bTagsJP = *(bTagHandleJP.product());

  // create the list of jets for the JEC study
  Double_t sigcut = 4;
  Double_t etcut  = 0.15;

  for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
    if( ZbbUtils::isJetIdOk((*jet),"loose")                                        // 1. loose jet id
        //&& VtxAssociatorsUtils::jetVertex(theZVertex, *jet, 4, sigcut, etcut)    // 2. matched vith dilepton vtx; algo 2b (cut on Sum_pTvtx/Sum_pT_All associated by ref)
	) theJets.push_back(*jet);
  }

  //if ( theJets.size() < 1 ) return; // no jets associated to the dilepton vertex

  // sort by pT
  std::vector<pat::Jet> theSortedJets = ZbbUtils::sortJetsBypT(theJets);

  //if ( !( ZbbUtils::isGoodJet(theSortedJets[0],ec.bestZcandidate_,lCuts_) && ZbbUtils::isBJet(theSortedJets[0],bTagAlgoWP) ) ) return; // the leading jet is not b-tagged

  //fill Event stuff
  Event.run_Event=iEvent.eventAuxiliary().run();
  Event.lumisec_Event=iEvent.eventAuxiliary().luminosityBlock();
  Event.num_Event=iEvent.eventAuxiliary().id().event();
  Event.weight_Event=ZbbUtils::getTheWeight(ec,N_);
  Event.njet_Event=theSortedJets.size();
  Event.nvert_Event=nvvertex;

  //fill Z stuff
  Z.pT_Z=ec.bestZcandidate_.pt();
  Z.pX_Z=ec.bestZcandidate_.px();
  Z.pY_Z=ec.bestZcandidate_.py();
  Z.eta_Z=ec.bestZcandidate_.eta();
  Z.phi_Z=ec.bestZcandidate_.phi();
  Z.eta_d1=ec.bestZcandidate_.daughter(0)->eta();
  Z.phi_d1=ec.bestZcandidate_.daughter(0)->phi();
  Z.pt_d1=ec.bestZcandidate_.daughter(0)->pt();
  Z.eta_d2=ec.bestZcandidate_.daughter(1)->eta();
  Z.phi_d2=ec.bestZcandidate_.daughter(1)->phi();
  Z.pt_d2=ec.bestZcandidate_.daughter(1)->pt();

  Double_t bjet_ptsum = 0;
  Double_t bjet_ptsum_nosigcut = 0;
  Double_t bjet_ptsumx = 0;
  Double_t bjet_ptsumy = 0;
  Double_t bjet_ptsumall = 0;
  Double_t bjet_et = 0;
  
  UInt_t theBjetIntex_ = -100;

  // loop on the jets to extract the b-jets
  for (UInt_t j = 0; j < theSortedJets.size(); ++j){ 
    if ( ZbbUtils::isGoodJet(theSortedJets[j],ec.bestZcandidate_,lCuts_) && ZbbUtils::isBJet(theSortedJets[j],bTagAlgoWP) ) { //select the b-jet and be sure it's a good jet
      theBTaggedJets.push_back(theSortedJets[j]);
    }
  }
  
  // number of b-jets in the event
  Event.nbjet_Event = theBTaggedJets.size();

  Int_t nb = theBTaggedJets.size();
  // filling the branched vector
  for(Int_t k=0;k<nb;k++){

    btagjets.pT_bjets[k]  = theBTaggedJets[k].pt();
    btagjets.eta_bjets[k] = theBTaggedJets[k].eta();
    btagjets.phi_bjets[k] = theBTaggedJets[k].phi();
    
    // discriminants
    btagjets.ssvhediscr_bjets[k]  = theBTaggedJets[k].bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    btagjets.ssvhpdiscr_bjets[k]  = theBTaggedJets[k].bDiscriminator("simpleSecondaryVertexHighPurBJetTags");

    btagjets.tchediscr_bjets[k]   = theBTaggedJets[k].bDiscriminator("trackCountingHighEffBJetTags");
    btagjets.tchpdiscr_bjets[k]   = theBTaggedJets[k].bDiscriminator("trackCountingHighPurBJetTags");

    btagjets.csvdiscr_bjets[k]    = theBTaggedJets[k].bDiscriminator("combinedSecondaryVertexBJetTags");
    btagjets.jpdiscr_bjets[k]   = theBTaggedJets[k].bDiscriminator("jetProbabilityBJetTags");
   
    btagjets.jpdiscr_fromtags_bjets[k] = ZbbUtils::getDiscriminatorFromTags(bTagsJP,theBTaggedJets[k]);
   
    if ( ec.isMC() ) {
      btagjets.trueflavour_bjets[k] = theBTaggedJets[k].partonFlavour();
    }
    
    // if(ZbbUtils::isJetPartonMatched(*theBTagJet,genParticlesCollection,genJets,5,0.5)) MCTrueFlavour=5;
    //       else if(ZbbUtils::isJetPartonMatched(*theBTagJet,genParticlesCollection,genJets,4,0.5) && 
    // 	      ZbbUtils::isJetPartonMatched(*theBTagJet,genParticlesCollection,genJets,5,0.5)) MCTrueFlavour=4; 
    //       else MCTrueFlavour=1;
    //       btagjets.trueflavour_bjets[k] = MCTrueFlavour;
    
    //sv mass  
    const reco::SecondaryVertexTagInfo &svTagInfo = *(theBTaggedJets[k].tagInfoSecondaryVertex("secondaryVertex"));
    
    if (svTagInfo.nVertices() > 0){
      btagjets.svmass_bjets[k]    = svTagInfo.secondaryVertex(0).p4().mass();
    }
  }  

  for (UInt_t j = 0; j < theSortedJets.size(); ++j){ 
    if ( ZbbUtils::isGoodJet(theSortedJets[j],ec.bestZcandidate_,lCuts_) && ZbbUtils::isBJet(theSortedJets[j],bTagAlgoWP) ) { //select the b-jet and be sure it's a good jet

      theBjetIntex_ = j;

      //fill bjet stuff  
      bjet.pT_bjet  = theSortedJets[j].pt();
      bjet.eta_bjet = theSortedJets[j].eta();
      bjet.phi_bjet = theSortedJets[j].phi();
      bjet_et       = theSortedJets[j].et();
      UInt_t bjet_constsize = theSortedJets[j].getPFConstituents().size();
   
      for(UInt_t i=0; i<bjet_constsize; i++){

	reco::PFCandidatePtr jetPFConst = theSortedJets[j].getPFConstituent(i);
 
	if (!(jetPFConst.isAvailable()) || (jetPFConst.isNull())){
	  continue;
	}
	if (jetPFConst->trackRef().isNull()){
	  continue;
	}
	if (jetPFConst->muonRef().isNonnull()  || jetPFConst->gsfTrackRef().isNonnull()  ){
	  continue;
	}
	
	for( reco::Vertex::trackRef_iterator tk = theZVertex.tracks_begin(); tk<theZVertex.tracks_end(); tk++ ){
	  if(tk->key() == jetPFConst->trackRef().key()) {
	    bjet_ptsum_nosigcut += jetPFConst->pt();
	  } // closes if tkref = pfcandidate 
	} // closes loop on vertex tracks
	
	Double_t distance = (jetPFConst->vz() - theZVertex.z());
	Double_t error = pow( (pow(jetPFConst->trackRef()->dzError(),2) + pow(theZVertex.zError(),2)),0.5);
	//error = vertex.zError();
	Double_t sig = distance/error;
	
	if (fabs(sig)<sigcut){
	  bjet_ptsum  += jetPFConst->pt();
	  bjet_ptsumx += jetPFConst->px();
	  bjet_ptsumy += jetPFConst->py();
	}
	bjet_ptsumall+= jetPFConst->pt();
      }
      break;
    }// if btagged and good jet
    else continue;
  }//loop over sorted jets

  bjet.frac1_bjet  = bjet_ptsum/bjet_et;
  bjet.frac2_bjet  = bjet_ptsum/bjet_ptsumall;
  bjet.frac3_bjet  = pow((pow(bjet_ptsumx,2) +pow(bjet_ptsumy,2)),0.5)/bjet_et;
  bjet.frac2b_bjet = bjet_ptsum_nosigcut/bjet_ptsumall;
  
  Double_t jet2_ptsum = 0;
  Double_t jet2_ptsum_nosigcut = 0;
  Double_t jet2_ptsumx = 0;
  Double_t jet2_ptsumy = 0;
  Double_t jet2_ptsumall = 0;
  Double_t jet2_et = 0;
  
  for (UInt_t j = 0; j < theSortedJets.size(); ++j){ 
    if ( /*fabs(theSortedJets[j].pt() - bjet.pT_bjet) > 0.01*/ 
	j!=theBjetIntex_ && VtxAssociatorsUtils::jetVertex(theZVertex,theSortedJets[j],4,sigcut,etcut)){ //ensure the jet is not the tagged one

      //fill jet2 stuff
      jet2.eta_jet2 = theSortedJets[j].eta();
      jet2.phi_jet2 = theSortedJets[j].phi();
      jet2.pT_jet2  = theSortedJets[j].pt();
      jet2_et       = theSortedJets[j].et();
      UInt_t jet2_constsize = theSortedJets[j].getPFConstituents().size();
            
      for(UInt_t i=0; i<jet2_constsize; i++){

	reco::PFCandidatePtr jetPFConst = theSortedJets[j].getPFConstituent(i);

	if (!(jetPFConst.isAvailable()) || (jetPFConst.isNull())){
	  continue;
	}
	if (jetPFConst->trackRef().isNull()){
	  continue;
	}
	if (jetPFConst->muonRef().isNonnull() || jetPFConst->gsfTrackRef().isNonnull()){
	  continue;
	}

	for( reco::Vertex::trackRef_iterator tk = theZVertex.tracks_begin(); tk<theZVertex.tracks_end(); tk++ ){
	  if(tk->key() == jetPFConst->trackRef().key()) {
	    jet2_ptsum_nosigcut += jetPFConst->pt();
	  } // closes if tkref = pfcandidate 
	} // closes loop on vertex tracks

	Double_t distance = (jetPFConst->vz() - theZVertex.z());
	Double_t error = pow( (pow(jetPFConst->trackRef()->dzError(),2) + pow(theZVertex.zError(),2)),0.5);
	//error = vertex.zError();
	Double_t sig = distance/error;
	
	if (fabs(sig)<sigcut){
	  jet2_ptsum  += jetPFConst->pt();
	  jet2_ptsumx += jetPFConst->px();
	  jet2_ptsumy += jetPFConst->py();
	}
	jet2_ptsumall+= jetPFConst->pt();
      }
      jet2.frac1_jet2  = jet2_ptsum/jet2_et;
      jet2.frac2_jet2  = jet2_ptsum/jet2_ptsumall;
      jet2.frac3_jet2  = pow((pow(jet2_ptsumx,2) +pow(jet2_ptsumy,2)),0.5)/jet2_et;
      jet2.frac2b_jet2 = jet2_ptsum_nosigcut/jet2_ptsumall;
      break;
    }//if subleading jet in the event
    else continue;
  }//loop over sorted Jets
  
  //fill MET stuff
  MET.pT_MET=(*mets)[0].pt();
  MET.px_MET=(*mets)[0].px();
  MET.py_MET=(*mets)[0].py();
  
  ZbbNtuple->Fill();

  return;

}

 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ZbbMatchedCandidate methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ZbbMatchedCandidate::book(){
  edm::Service<TFileService> fs;
  std::string DirName = &(*N_);
  TFileDirectory HistoVars = fs->mkdir(DirName);
  
  TH1F::SetDefaultSumw2(kTRUE);
 
  // all flavours
  h_zllmass_         = HistoVars.make<TH1F>(N_+"_zllmass","l^{+}l^{-} mass; M(l^{+}l^{-}) (GeV)",60,60,120);
  h_zllpt_           = HistoVars.make<TH1F>(N_+"_zllpt","l^{+}l^{-} p_{T}; p_{T}(l^{+}l^{-}) (GeV)",100,0,200);
  h_zllrapidity_     = HistoVars.make<TH1F>(N_+"_zllrapidity","l^{+}l^{-} y; y(l^{+}l^{-}) (GeV)",100,-3,3);
  h_zllCosThetaStar_ = HistoVars.make<TH1F>(N_+"_llCosThetaStar","Cosine of Z decay frame l^{+} angle;  cos(#theta*_{l^{+}})",20,-1.,1.);
  h_lept1pt_         = HistoVars.make<TH1F>(N_+"_lept1pt","leading lepton p_{T}; p^{lead lept}_{T} (GeV)",100,0,200);
  h_lept2pt_         = HistoVars.make<TH1F>(N_+"_lept2pt","subleading lepton p_{T}; p^{2nd lept}_{T} (GeV)",100,0,200);
  h_lept1eta_        = HistoVars.make<TH1F>(N_+"_lept1eta","leading lepton #eta; #eta^{lead lept}",25,-2.5,2.5);
  h_lept2eta_        = HistoVars.make<TH1F>(N_+"_lept2eta","subleading lepton #eta; #eta^{2nd lept}",25,-2.5,2.5);
  h_nb_              = HistoVars.make<TH1F>(N_+"_nb","b-tagged jet count; N_{b-jets}",10,-0.5,9.5);
  h_sumbHtOversumHt_ = HistoVars.make<TH1F>(N_+"_sumbHtOversumHt","Scalar sum b-jet p_{T}/ Scalar sum of all jet p_{T};(sum p^{b}_{T})/(sum all p^{j}_{T})",50,0.,1.);
  h_bjetmass_        = HistoVars.make<TH1F>(N_+"_bjetmass","all b-jet mass; all b-jets m_{bj} (GeV)",100,0.,100.);
  h_bjet1pt_         = HistoVars.make<TH1F>(N_+"_bjet1pt","leading bjet p_{T}; leading b-jet p_{T} (GeV)",50,15,215);
  h_bjet1eta_        = HistoVars.make<TH1F>(N_+"_bjet1eta","leading bjet #eta; leading b-jet #eta",25,-2.5,2.5);
  h_bjet2pt_         = HistoVars.make<TH1F>(N_+"_bjet2pt","subleading bjet p_{T}; subleading b-jet p_{T} (GeV)",50,15,215);
  h_bjet2eta_        = HistoVars.make<TH1F>(N_+"_bjet2eta","subleading bjet #eta; subleading b-jet #eta",15,-2.5,2.5);
  h_bbM_             = HistoVars.make<TH1F>(N_+"_bbjetM","b#bar{b} invariant mass; M(b#bar{b}) (GeV)",50,0,1000);
  h_bbDeltaR_        = HistoVars.make<TH1F>(N_+"_bbDeltaR","#Delta R between b-jets for b#bar{b} events; #DeltaR(j_{b},j_{#bar{b}})",30,0,5.); 
  
  // only muons
  h_zmumumass_         = HistoVars.make<TH1F>(N_+"_zmumumass","#mu^{+}#mu^{-} mass; M(#mu^{+}#mu^{-}) (GeV)",30,60,120);
  h_zmumupt_           = HistoVars.make<TH1F>(N_+"_zmumupt","#mu^{+}#mu^{-} p_{T}; p_{T}(#mu^{+}#mu^{-}) (GeV)",50,0,200);
  h_zmumurapidity_     = HistoVars.make<TH1F>(N_+"_zmumurapidity","#mu^{+}#mu^{-} y; y(#mu^{+}#mu^{-}) (GeV)",25,-3,3);
  h_zmumuCosThetaStar_ = HistoVars.make<TH1F>(N_+"_mumuCosThetaStar","Cosine of Z decay frame #mu^{+} angle;  cos(#theta*_{#mu^{+}})",20,-1.,1.);
  h_mu1pt_             = HistoVars.make<TH1F>(N_+"_mu1pt","leading muon p_{T}; p^{lead #mu}_{T} (GeV)",50,0,200);
  h_mu2pt_             = HistoVars.make<TH1F>(N_+"_mu2pt","subleading muon p_{T}; p^{2nd #mu}_{T} (GeV)",50,0,200);
  h_mu1eta_            = HistoVars.make<TH1F>(N_+"_mu1eta","leading muon #eta; #eta^{lead #mu}",25,-2.5,2.5);
  h_mu2eta_            = HistoVars.make<TH1F>(N_+"_mu2eta","subleading muon #eta; #eta^{2nd #mu}",25,-2.5,2.5);
  h_mmnb_              = HistoVars.make<TH1F>(N_+"_mmnb","b-tagged jet count; N_{b-jets}",10,-0.5,9.5);
  h_mmsumbHtOversumHt_ = HistoVars.make<TH1F>(N_+"_mmsumbHtOversumHt","Scalar sum b-jet p_{T}/ Scalar sum of all jet p_{T};(sum p^{b}_{T})/(sum p^{j}_{T})",25,0.,1.);
  h_mmbjetmass_        = HistoVars.make<TH1F>(N_+"_mmbjetmass","all b-jet mass; all b-jets m_{bj} (GeV)",50,0.,100.);
  h_mmbjet1pt_         = HistoVars.make<TH1F>(N_+"_mmbjet1pt","leading bjet p_{T}; leading b-jet p_{T} (GeV)",50,15,215);
  h_mmbjet1eta_        = HistoVars.make<TH1F>(N_+"_mmbjet1eta","leading bjet #eta; leading b-jet #eta",25,-2.5,2.5);
  h_mmbjet2pt_         = HistoVars.make<TH1F>(N_+"_mmbjet2pt","subleading bjet p_{T}; subleading b-jet p_{T} (GeV)",50,15,215);
  h_mmbjet2eta_        = HistoVars.make<TH1F>(N_+"_mmbjet2eta","subleading bjet #eta; subleading b-jet #eta",15,-2.5,2.5);
  h_mmbbM_             = HistoVars.make<TH1F>(N_+"_mmbbjetM","b#bar{b} invariant mass; M(b#bar{b}) (GeV)",50,0,1000);
  h_mmbbDeltaR_        = HistoVars.make<TH1F>(N_+"_mmbbDeltaR","#Delta R between b-jets for b#bar{b} events; #DeltaR(j_{b},j_{#bar{b}})",30,0,5.); 

  // only electrons
  h_zeleelemass_         = HistoVars.make<TH1F>(N_+"_zeleelemass","e^{+}e^{-} mass; M(e^{+}e^{-}) (GeV)",30,60,120);
  h_zeleelept_           = HistoVars.make<TH1F>(N_+"_zeleelept","e^{+}e^{-} p_{T}; p_{T}(e^{+}e^{-}) (GeV)",50,0,200);
  h_zeleelerapidity_     = HistoVars.make<TH1F>(N_+"_zeleelerapidity","e^{+}e^{-} y; y(e^{+}e^{-}) (GeV)",25,-3,3);
  h_zeleeleCosThetaStar_ = HistoVars.make<TH1F>(N_+"_eleeleCosThetaStar","Cosine of Z decay frame e^{+} angle;  cos(#theta*_{e^{+}})",20,-1.,1.);
  h_ele1pt_              = HistoVars.make<TH1F>(N_+"_ele1pt","leading electron p_{T}; p^{lead ele}_{T} (GeV)",50,0,200);
  h_ele2pt_              = HistoVars.make<TH1F>(N_+"_ele2pt","subleading electron p_{T}; p^{2nd ele}_{T} (GeV)",50,0,200);
  h_ele1eta_             = HistoVars.make<TH1F>(N_+"_ele1eta","leading electron #eta; #eta^{lead ele}",25,-2.5,2.5);
  h_ele2eta_             = HistoVars.make<TH1F>(N_+"_ele2eta","subleading electron #eta; #eta^{2nd ele}",25,-2.5,2.5);
  h_eenb_                = HistoVars.make<TH1F>(N_+"_eenb","b-tagged jet count; N_{b-jets}",10,-0.5,9.5);
  h_eesumbHtOversumHt_   = HistoVars.make<TH1F>(N_+"_eesumbHtOversumHt","Scalar sum b-jet p_{T}/ Scalar sum of all jet p_{T};(sum p^{b}_{T})/(sum all p^{j}_{T})",25,0.,1.);
  h_eebjetmass_          = HistoVars.make<TH1F>(N_+"_eebjetmass","all b-jet mass; all b-jets m_{bj} (GeV)",50,0.,100.);
  h_eebjet1pt_           = HistoVars.make<TH1F>(N_+"_eebjet1pt","leading bjet p_{T}; leading b-jet p_{T} (GeV)",50,15,215);
  h_eebjet1eta_          = HistoVars.make<TH1F>(N_+"_eebjet1eta","leading bjet #eta; leading b-jet #eta",25,-2.5,2.5);
  h_eebjet2pt_           = HistoVars.make<TH1F>(N_+"_eebjet2pt","subleading bjet p_{T}; subleading b-jet p_{T} (GeV)",50,15,215);
  h_eebjet2eta_          = HistoVars.make<TH1F>(N_+"_eebjet2eta","subleading bjet #eta; subleading b-jet #eta",15,-2.5,2.5);
  h_eebbM_               = HistoVars.make<TH1F>(N_+"_eebbjetM","b#bar{b} invariant mass; M(b#bar{b}) (GeV)",50,0,1000);
  h_eebbDeltaR_          = HistoVars.make<TH1F>(N_+"_eebbDeltaR","#Delta R between b-jets for b#bar{b} events; #DeltaR(j_{b},j_{#bar{b}})",30,0,5.); 
  
}

//--------------------------------- filling method for MC ----------------------------------------------------------
void ZbbMatchedCandidate::fillMC(const EventCategory& ec_mu, const EventCategory& ec_ele, edm::Handle<edm::View<pat::Jet> > jets, std::string bTagAlgoWP,
				 edm::Handle<reco::GenParticleCollection> genParticlesCollection,const reco::GenJetCollection & genJets){

  int nBjets = 0;
  
  Bool_t thecut_;   
  EventCategory ec; 
  Bool_t isMuChannel_(false);
  
  if(ec_mu.isZLL()){
    thecut_= ZbbUtils::ProduceTheCut(ec_mu,N_,"");
    ec = ec_mu;
    isMuChannel_ = true;
  } else if(ec_ele.isZLL()) {
    thecut_ = ZbbUtils::ProduceTheCut(ec_ele,N_,"");
    ec = ec_ele;
    isMuChannel_ = false;
  } else {
    return;
  }
  
  Double_t w = ZbbUtils::getTheWeight(ec,N_);

  // vector of b-tagged jets
  std::vector<pat::Jet> theBtaggedJets;
  // vector of good Jets
  std::vector<pat::Jet> theGoodJets;

  if(thecut_){
    
    if(ZbbUtils::isZLLCandidateMCMatched(ec.bestZcandidate_,genParticlesCollection)){
      
      TLorentzVector ZP4_(ec.bestZcandidate_.px(),ec.bestZcandidate_.py(),ec.bestZcandidate_.pz(),ec.bestZcandidate_.energy());
      
      h_zllmass_->Fill(ZP4_.M(),w); 
      h_zllpt_->Fill(std::min(ZP4_.Perp(),199.9),w);   
      h_zllrapidity_->Fill(ZP4_.Rapidity(),w);
      
      if(isMuChannel_){
	h_zmumumass_->Fill(ZP4_.M(),w); 
	h_zmumupt_->Fill(std::min(ZP4_.Perp(),199.9),w);   
	h_zmumurapidity_->Fill(ZP4_.Rapidity(),w);	
      } else {
	h_zeleelemass_->Fill(ZP4_.M(),w); 
	h_zeleelept_->Fill(std::min(ZP4_.Perp(),199.9),w);   
	h_zeleelerapidity_->Fill(ZP4_.Rapidity(),w);
      }

      const reco::Candidate* lepton1 = ec.bestZcandidate_.daughter(0);
      const reco::Candidate* lepton2 = ec.bestZcandidate_.daughter(1); 
      
      h_zllCosThetaStar_->Fill(ZbbUtils::CostCS(ec.bestZcandidate_),w);
      h_lept1pt_->Fill(std::min(lepton1->pt(),199.9),w); 
      h_lept2pt_->Fill(std::min(lepton2->pt(),199.9),w); 
      h_lept1eta_->Fill(lepton1->eta(),w);
      h_lept2eta_->Fill(lepton2->eta(),w);

      if(isMuChannel_){	      
	h_zmumuCosThetaStar_->Fill(ZbbUtils::CostCS(ec.bestZcandidate_),w);
	h_mu1pt_->Fill(std::min(lepton1->pt(),199.9),w); 
	h_mu2pt_->Fill(std::min(lepton2->pt(),199.9),w); 
	h_mu1eta_->Fill(lepton1->eta(),w);
	h_mu2eta_->Fill(lepton2->eta(),w);
      } else {
	h_zeleeleCosThetaStar_->Fill(ZbbUtils::CostCS(ec.bestZcandidate_),w);
	h_ele1pt_->Fill(std::min(lepton1->pt(),199.9),w); 
	h_ele2pt_->Fill(std::min(lepton2->pt(),199.9),w); 
	h_ele1eta_->Fill(lepton1->eta(),w);
	h_ele2eta_->Fill(lepton2->eta(),w);
      }

      for(edm::View<pat::Jet>::const_iterator jetcand=jets->begin(); jetcand!=jets->end(); ++jetcand){
	
	if(ZbbUtils::isJetIdOk((*jetcand),"loose") && ZbbUtils::isGoodJet((*jetcand),ec.bestZcandidate_,lCuts_)){	
	  theGoodJets.push_back(*jetcand);
	  if(ZbbUtils::isBJet((*jetcand),bTagAlgoWP)){
	    if(ZbbUtils::isBJetMCMatched((*jetcand),genParticlesCollection,genJets)){
	      theBtaggedJets.push_back(*jetcand);
	      nBjets++;
	    } // close if matched b-jet
	  } // close if b-tag jet     
	} // close if good jet
      } // close the for on jets
      
      Double_t allJetsSumPt_(0);
      Double_t allbJetsSumPt_(0);
      
      for(UInt_t i=0; i<theGoodJets.size();i++){
	allJetsSumPt_+=theGoodJets[i].pt();
      }
      
      h_nb_->Fill(nBjets,w);
      if(isMuChannel_){
	h_mmnb_->Fill(nBjets,w);
      } else{
	h_eenb_->Fill(nBjets,w);
      }
      
      if(nBjets>0){
	for(UInt_t i=0; i<theBtaggedJets.size();i++){
	  allbJetsSumPt_+=theBtaggedJets[i].pt();
	}
	
	h_sumbHtOversumHt_->Fill((allbJetsSumPt_/allJetsSumPt_),w);
	h_bjetmass_->Fill(theBtaggedJets[0].mass(),w);
	h_bjet1pt_->Fill(std::min(theBtaggedJets[0].pt(),214.99),w); 
	h_bjet1eta_->Fill(theBtaggedJets[0].eta(),w);      

	if(isMuChannel_){
	  h_mmsumbHtOversumHt_->Fill((allbJetsSumPt_/allJetsSumPt_),w);
	  h_mmbjetmass_->Fill(theBtaggedJets[0].mass(),w);
	  h_mmbjet1pt_->Fill(std::min(theBtaggedJets[0].pt(),214.99),w); 
	  h_mmbjet1eta_->Fill(theBtaggedJets[0].eta(),w);   
	} else {
	  h_eesumbHtOversumHt_->Fill((allbJetsSumPt_/allJetsSumPt_),w);
	  h_eebjetmass_->Fill(theBtaggedJets[0].mass(),w);
	  h_eebjet1pt_->Fill(std::min(theBtaggedJets[0].pt(),214.99),w); 
	  h_eebjet1eta_->Fill(theBtaggedJets[0].eta(),w);   
	}

	TLorentzVector b1P4_(theBtaggedJets[0].px(),theBtaggedJets[0].py(),theBtaggedJets[0].pz(),theBtaggedJets[0].energy());
	
	if(nBjets>1){
	  h_bjet2pt_->Fill(std::min(theBtaggedJets[1].pt(),214.99),w);        
	  h_bjet2eta_->Fill(theBtaggedJets[1].eta(),w);
	  TLorentzVector b2P4_(theBtaggedJets[1].px(),theBtaggedJets[1].py(),theBtaggedJets[1].pz(),theBtaggedJets[1].energy());
	  TLorentzVector bbP4_ = b1P4_ + b2P4_;
	  h_bbM_->Fill(bbP4_.M(),w);
	  h_bbDeltaR_->Fill(ROOT::Math::VectorUtil::DeltaR(theBtaggedJets[0].momentum(),theBtaggedJets[1].momentum()),w);

	  if(isMuChannel_){ 
	    h_mmbjet2pt_->Fill(std::min(theBtaggedJets[1].pt(),214.99),w);        
	    h_mmbjet2eta_->Fill(theBtaggedJets[1].eta(),w);
	    h_mmbbM_->Fill(bbP4_.M(),w);
	    h_mmbbDeltaR_->Fill(ROOT::Math::VectorUtil::DeltaR(theBtaggedJets[0].momentum(),theBtaggedJets[1].momentum()),w);
	  } else {
	    h_eebjet2pt_->Fill(std::min(theBtaggedJets[1].pt(),214.99),w);        
	    h_eebjet2eta_->Fill(theBtaggedJets[1].eta(),w);
	    h_eebbM_->Fill(bbP4_.M(),w);
	    h_eebbDeltaR_->Fill(ROOT::Math::VectorUtil::DeltaR(theBtaggedJets[0].momentum(),theBtaggedJets[1].momentum()),w);
	  }
	} //if at least 2 b-tagged jets
      } // if at least 1 b-tagged jets
    } // if the Z is MC matched
  } // closes if right category
}

//---------------------------------------- filling method for DATA -------------------------------------------------
void ZbbMatchedCandidate::fillDATA(const EventCategory& ec_mu, const EventCategory& ec_ele, edm::Handle<edm::View<pat::Jet> > jets, std::string bTagAlgoWP){

  int nBjets = 0;
  
  Bool_t thecut_;   
  EventCategory ec; 
  Bool_t isMuChannel_(false);
  
  if(ec_mu.isZLL()){
    thecut_= ZbbUtils::ProduceTheCut(ec_mu,N_,"");
    ec = ec_mu;
    isMuChannel_ = true;
  } else if(ec_ele.isZLL()) {
    thecut_ = ZbbUtils::ProduceTheCut(ec_ele,N_,"");
    ec = ec_ele;
    isMuChannel_ = false;
  } else {
    return;
  }
  
  Double_t w = ZbbUtils::getTheWeight(ec,N_);
  
  // vector of b-tagged jets
  std::vector<pat::Jet> theBtaggedJets;
  // vector of good Jets
  std::vector<pat::Jet> theGoodJets;
  
  if(thecut_){
    
    TLorentzVector ZP4_(ec.bestZcandidate_.px(),ec.bestZcandidate_.py(),ec.bestZcandidate_.pz(),ec.bestZcandidate_.energy());
    
    h_zllmass_->Fill(ZP4_.M(),w); 
    h_zllpt_->Fill(std::min(ZP4_.Perp(),199.9),w);   
    h_zllrapidity_->Fill(ZP4_.Rapidity(),w);
    
    if(isMuChannel_){
      h_zmumumass_->Fill(ZP4_.M(),w); 
      h_zmumupt_->Fill(std::min(ZP4_.Perp(),199.9),w);   
      h_zmumurapidity_->Fill(ZP4_.Rapidity(),w);	
    } else {
      h_zeleelemass_->Fill(ZP4_.M(),w); 
      h_zeleelept_->Fill(std::min(ZP4_.Perp(),199.9),w);   
      h_zeleelerapidity_->Fill(ZP4_.Rapidity(),w);
    }
    
    const reco::Candidate* lepton1 = ec.bestZcandidate_.daughter(0);
    const reco::Candidate* lepton2 = ec.bestZcandidate_.daughter(1); 
    
    h_zllCosThetaStar_->Fill(ZbbUtils::CostCS(ec.bestZcandidate_),w);
    h_lept1pt_->Fill(std::min(lepton1->pt(),199.9),w); 
    h_lept2pt_->Fill(std::min(lepton2->pt(),199.9),w); 
    h_lept1eta_->Fill(lepton1->eta(),w);
    h_lept2eta_->Fill(lepton2->eta(),w);
    
    if(isMuChannel_){	      
      h_zmumuCosThetaStar_->Fill(ZbbUtils::CostCS(ec.bestZcandidate_),w);
      h_mu1pt_->Fill(std::min(lepton1->pt(),199.9),w); 
      h_mu2pt_->Fill(std::min(lepton2->pt(),199.9),w); 
      h_mu1eta_->Fill(lepton1->eta(),w);
      h_mu2eta_->Fill(lepton2->eta(),w);
    } else {
      h_zeleeleCosThetaStar_->Fill(ZbbUtils::CostCS(ec.bestZcandidate_),w);
      h_ele1pt_->Fill(std::min(lepton1->pt(),199.9),w); 
      h_ele2pt_->Fill(std::min(lepton2->pt(),199.9),w); 
      h_ele1eta_->Fill(lepton1->eta(),w);
      h_ele2eta_->Fill(lepton2->eta(),w);
    }
    
    for(edm::View<pat::Jet>::const_iterator jetcand=jets->begin(); jetcand!=jets->end(); ++jetcand){
      
      if(ZbbUtils::isJetIdOk((*jetcand),"loose") && ZbbUtils::isGoodJet((*jetcand),ec.bestZcandidate_,lCuts_)){	
	theGoodJets.push_back(*jetcand);
	if(ZbbUtils::isBJet((*jetcand),bTagAlgoWP)){
	  theBtaggedJets.push_back(*jetcand);
	  nBjets++;
	} // close if b-tag jet     
      } // close if good jet
    } // close the for on jets
    
    Double_t allJetsSumPt_(0);
    Double_t allbJetsSumPt_(0);
    
    for(UInt_t i=0; i<theGoodJets.size();i++){
      allJetsSumPt_+=theGoodJets[i].pt();
    }
    
    h_nb_->Fill(nBjets,w);
    if(isMuChannel_){
      h_mmnb_->Fill(nBjets,w);
    } else{
      h_eenb_->Fill(nBjets,w);
    }
    
    if(nBjets>0){
      for(UInt_t i=0; i<theBtaggedJets.size();i++){
	allbJetsSumPt_+=theBtaggedJets[i].pt();
      }
      
      h_sumbHtOversumHt_->Fill((allbJetsSumPt_/allJetsSumPt_),w);
      h_bjetmass_->Fill(theBtaggedJets[0].mass(),w);
      h_bjet1pt_->Fill(std::min(theBtaggedJets[0].pt(),214.99),w); 
      h_bjet1eta_->Fill(theBtaggedJets[0].eta(),w);      
      
      if(isMuChannel_){
	h_mmsumbHtOversumHt_->Fill((allbJetsSumPt_/allJetsSumPt_),w);
	h_mmbjetmass_->Fill(theBtaggedJets[0].mass(),w);
	h_mmbjet1pt_->Fill(std::min(theBtaggedJets[0].pt(),214.99),w); 
	h_mmbjet1eta_->Fill(theBtaggedJets[0].eta(),w);   
      } else {
	h_eesumbHtOversumHt_->Fill((allbJetsSumPt_/allJetsSumPt_),w);
	h_eebjetmass_->Fill(theBtaggedJets[0].mass(),w);
	h_eebjet1pt_->Fill(std::min(theBtaggedJets[0].pt(),214.99),w); 
	h_eebjet1eta_->Fill(theBtaggedJets[0].eta(),w);   
      }
      
      TLorentzVector b1P4_(theBtaggedJets[0].px(),theBtaggedJets[0].py(),theBtaggedJets[0].pz(),theBtaggedJets[0].energy());
      
      if(nBjets>1){
	h_bjet2pt_->Fill(std::min(theBtaggedJets[1].pt(),214.99),w);        
	h_bjet2eta_->Fill(theBtaggedJets[1].eta(),w);
	TLorentzVector b2P4_(theBtaggedJets[1].px(),theBtaggedJets[1].py(),theBtaggedJets[1].pz(),theBtaggedJets[1].energy());
	TLorentzVector bbP4_ = b1P4_ + b2P4_;
	h_bbM_->Fill(bbP4_.M(),w);
	h_bbDeltaR_->Fill(ROOT::Math::VectorUtil::DeltaR(theBtaggedJets[0].momentum(),theBtaggedJets[1].momentum()),w);
	
	if(isMuChannel_){ 
	  h_mmbjet2pt_->Fill(std::min(theBtaggedJets[1].pt(),214.99),w);        
	  h_mmbjet2eta_->Fill(theBtaggedJets[1].eta(),w);
	  h_mmbbM_->Fill(bbP4_.M(),w);
	  h_mmbbDeltaR_->Fill(ROOT::Math::VectorUtil::DeltaR(theBtaggedJets[0].momentum(),theBtaggedJets[1].momentum()),w);
	} else {
	  h_eebjet2pt_->Fill(std::min(theBtaggedJets[1].pt(),214.99),w);        
	  h_eebjet2eta_->Fill(theBtaggedJets[1].eta(),w);
	  h_eebbM_->Fill(bbP4_.M(),w);
	  h_eebbDeltaR_->Fill(ROOT::Math::VectorUtil::DeltaR(theBtaggedJets[0].momentum(),theBtaggedJets[1].momentum()),w);
	}
      } //if at least 2 b-tagged jets
    } // if at least 1 b-tagged jets
  } // closes if right category
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ZbbCollections methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ZbbCollections::book(){
  edm::Service<TFileService> fs;
  std::string DirName = &(*N_);
  TFileDirectory HistoVars = fs->mkdir(DirName);
  
  h_eventWeight_ = HistoVars.make<TH1F>(N_+"_eventWeight","event weight; w_{event}",100,0.5,1.5);
     
  TH1F::SetDefaultSumw2(kTRUE);

  // lepton specific
  h_zllmass_   = HistoVars.make<TH1F>(N_+"_zmass","l^{+}l^{-} mass; M(l^{+}l^{-}) (GeV)",60,60,120);
  h_zllpt_     = HistoVars.make<TH1F>(N_+"_zpt","l^{+}l^{-} p_{T}; p_{T}(l^{+}l^{-}) (GeV)",50,0,150);
  h_zlleta_    = HistoVars.make<TH1F>(N_+"_zeta","l^{+}l^{-} #eta; #eta(l^{+}l^{-}) (GeV)",60,-3,3);
  h_lldeltaPhi_= HistoVars.make<TH1F>(N_+"_lldeltaPhi","#Delta #phi of di-lepton ; #Delta#phi(l^{+}l^{-}) (rad)",60,0.,TMath::Pi());  
  h_lldeltaR_  = HistoVars.make<TH1F>(N_+"_lldeltaR","#Delta R of di-lepton; #DeltaR(l^{+}l^{-})",100,0,6); 
  h_zllCosThetaCS_   = HistoVars.make<TH1F>(N_+"_llCosThetaCS","Cosine of Collins-Soper frame angle; #theta_{CS}(l^{+})",20,-1.,1.); 
  h_zllCosThetaStar_ = HistoVars.make<TH1F>(N_+"_llCosThetaStar","Cosine of Z decay frame l^{+} angle; #theta*(l^{+})",20,-1.,1.); 

  // muon specific
  h_zmumumass_ = HistoVars.make<TH1F>(N_+"_zmassMu","#mu^{+}#mu^{-} mass; M(#mu^{+}#mu^{-}) (GeV)",60,60,120);
  h_zmumupt_   = HistoVars.make<TH1F>(N_+"_zptMu","#mu^{+}#mu^{-} p_{T}; p_{T}(#mu^{+}#mu^{-}) (GeV)",50,0,100);
  h_zmumueta_   = HistoVars.make<TH1F>(N_+"_zetaMu","#mu^{+}#mu^{-} #eta; #eta(#mu^{+}#mu^{-}) (GeV)",60,-3,3);
  h_mu1pt_     = HistoVars.make<TH1F>(N_+"_mu1pt","leading muon p_{T}; p^{lead #mu}_{T} (GeV)",50,0,100);
  h_mu2pt_     = HistoVars.make<TH1F>(N_+"_mu2pt","subleading muon p_{T}; p^{2nd #mu}_{T} (GeV)",50,0,100);
  h_mu1eta_    = HistoVars.make<TH1F>(N_+"_mu1eta","leading muon #eta; #eta^{lead #mu}",50,-2.5,2.5);
  h_mu2eta_    = HistoVars.make<TH1F>(N_+"_mu2eta","subleading muon #eta; #eta^{2nd #mu}",50,-2.5,2.5);
  h_mumudeltaPhi_= HistoVars.make<TH1F>(N_+"_mumudeltaPhi","#Delta #phi of di-muon ; #Delta#phi(#mu^{+}#mu^{-}) (rad)",60,0.,TMath::Pi());  
  h_mumudeltaR_  = HistoVars.make<TH1F>(N_+"_mumudeltaR","#Delta R of di-muon; #DeltaR(#mu^{+}#mu^{-})",100,0,6); 
  h_zmumuCosThetaCS_   = HistoVars.make<TH1F>(N_+"_mumuCosThetaCS","Cosine of Collins-Soper frame angle; #theta_{CS}(#mu^{+})",20,-1.,1.); 
  h_zmumuCosThetaStar_ = HistoVars.make<TH1F>(N_+"_mumuCosThetaStar","Cosine of Z decay frame #mu^{+} angle; #theta*(#mu^{+})",20,-1.,1.); 
  
  // electron specific
  h_zeemass_  = HistoVars.make<TH1F>(N_+"_zmassEle","e^{+}e^{-} mass; M(e^{+}e^{-}) (GeV)",60,60,120);
  h_zeept_    = HistoVars.make<TH1F>(N_+"_zptEle","e^{+}e^{-} p_{T}; p_{T}(e^{+}e^{-}) (GeV)",50,0,100);
  h_zeeeta_   = HistoVars.make<TH1F>(N_+"_zetaEle","e^{+}e^{-} #eta; #eta(e^{+}e^{-}) (GeV)",60,-3,3);
  h_ele1pt_   = HistoVars.make<TH1F>(N_+"_ele1pt","leading electron p_{T};  p^{lead e}_{T} (GeV)",50,0,100);
  h_ele2pt_   = HistoVars.make<TH1F>(N_+"_ele2pt","subleading electron p_{T};  p^{2nd e}_{T} (GeV)",50,0,100);
  h_ele1eta_  = HistoVars.make<TH1F>(N_+"_ele1eta","leading electron #eta; #eta^{lead e}",50,-2.5,2.5);
  h_ele2eta_  = HistoVars.make<TH1F>(N_+"_ele2eta","subleading electron #eta; #eta^{2nd e} ",50,-2.5,2.5);
  h_eedeltaPhi_ = HistoVars.make<TH1F>(N_+"_eedeltaPhi","#Delta #phi of di-electron; #Delta#phi(e^{+}e^{-}) (rad)",60,0.,TMath::Pi());  
  h_eedeltaR_   = HistoVars.make<TH1F>(N_+"_eedeltaR","#Delta R of di-electron; #DeltaR(e^{+}e^{-})",100,0,6); 
  h_zeeCosThetaCS_   = HistoVars.make<TH1F>(N_+"_eeCosThetaCS","Cosine of Collins-Soper frame angle; #theta_{CS}(e^{+})",20,-1.,1.); 
  h_zeeCosThetaStar_ = HistoVars.make<TH1F>(N_+"_eeCosThetaStar","Cosine of Z decay frame e^{+} angle; #theta*(e^{+})",20,-1.,1.); 

  // jet specific
  h_nj_            = HistoVars.make<TH1F>(N_+"_nj","jet count; N_{jets}",15,-0.5,14.5);
  h_SSVHEdisc_     = HistoVars.make<TH1F>(N_+"_SSVHEdisc","SSVHEdisc; SSV HE discriminant",100,-2,10);
  h_SSVHPdisc_     = HistoVars.make<TH1F>(N_+"_SSVHPdisc","SSVHPdisc; SSV HP discriminant",100,-2,10);
  h_TCHEdisc_      = HistoVars.make<TH1F>(N_+"_TCHEdisc","TCHEdisc; TC HE discriminant",100,-10,10);
  h_TCHPdisc_      = HistoVars.make<TH1F>(N_+"_TCHPdisc","TCHPdisc; TC HP discriminant",100,-10,10);
  h_CSVdisc_       = HistoVars.make<TH1F>(N_+"_CSVdisc","CSVdisc; CSV discriminant",48,-0.1,1.1);
  h_JPdisc_        = HistoVars.make<TH1F>(N_+"_JPdisc","JPdisc; JP discriminant",48,0.,1.1);
  h_met_           = HistoVars.make<TH1F>(N_+"_MET","MET; #slash{E}_{T} (GeV)",100,0,200);
  h_phimet_        = HistoVars.make<TH1F>(N_+"_METphi","MET #phi; MET #phi (rad)",70,-TMath::Pi(),TMath::Pi());
  h_jetpt_         = HistoVars.make<TH1F>(N_+"_jetpt","All jets p_{T}; Jet p_{T} (GeV)",100,15,215);
  h_jeteta_        = HistoVars.make<TH1F>(N_+"_jeteta","All jets #eta; Jet #eta",25,-2.5, 2.5);
  h_jetphi_        = HistoVars.make<TH1F>(N_+"_jetphi","All jet #phi; Jet #phi (rad)",70,-TMath::Pi(),TMath::Pi());
  h_jetoverlapLept_= HistoVars.make<TH1F>(N_+"_jetoverlapLept","jets overlaps with leptons from Z",2,-0.5,1.5);

  TString jetOverlapLeptBinLabels[2] ={"no overlap","has overlap"};
  
  for(UInt_t bin=1; bin<=2; bin++){
    h_jetoverlapLept_->GetXaxis()->SetBinLabel(bin,jetOverlapLeptBinLabels[bin-1]);    
  }

  h_rawjetspt_     = HistoVars.make<TH1F>(N_+"_rawjetpt","All raw jets p_{T}; Jet p_{T} (GeV)",100,15,215);
  h_rawjetseta_    = HistoVars.make<TH1F>(N_+"_rawjeteta","All raw jets #eta; Jet #eta",25,-2.5, 2.5);

  h_jet1pt_        = HistoVars.make<TH1F>(N_+"_jet1pt","leading jet p_{T}; leading jet p_{T} (GeV)",100,15,515);
  h_jet1eta_       = HistoVars.make<TH1F>(N_+"_jet1eta","leading jet #eta; leading jet #eta",50,-2.5,2.5);
  h_jet2pt_        = HistoVars.make<TH1F>(N_+"_jet2pt","subleading jet p_{T}; subleading jet p_{T} (GeV)",100,15,515);
  h_jet2eta_       = HistoVars.make<TH1F>(N_+"_jet2eta","subleading jet #eta; subleading jet #eta",50,-2.5,2.5);
  h_nhf_           = HistoVars.make<TH1F>(N_+"_nhf","neutral hadron energy fraction; E_{n. had}/E_{tot}",101,0,1.01);
  h_nef_           = HistoVars.make<TH1F>(N_+"_nef","neutral EmEnergy fraction; E_{n. em}/E_{tot}",101,0,1.01);
  h_nconstituents_ = HistoVars.make<TH1F>(N_+"_npf","total multiplicity; jet total multiplicity",10,-0.5,9.5);
  h_chf_           = HistoVars.make<TH1F>(N_+"_chf","charged hadron energy fraction; E_{ch. had}/E_{tot}",101,0,1.01);
  h_nch_           = HistoVars.make<TH1F>(N_+"_nch","charged multiplicity; jet charged multiplicity",50,-0.5,49.5);
  h_cef_           = HistoVars.make<TH1F>(N_+"_cef","charged EmEnergy fraction; E_{ch. em}/E_{tot}",101,0,1.01);
  h_alljetDeltaR_  = HistoVars.make<TH1F>(N_+"_alljetDeltaR","#Delta R between all jet pairs; #DeltaR(j,j)",100,0,10); 
  h_alljetDeltaPhi_= HistoVars.make<TH1F>(N_+"_alljetDeltaPhi","#Delta #phi of all jet pairs; #Delta#phi(j,j) (rad)",60,0.,TMath::Pi());  
  h_jetid_         = HistoVars.make<TH1F>(N_+"_jetid","Jet Id level (none, loose, medium, tight)",4,-0.5,3.5);

  TString jetIdBinLabels[4] ={"none","loose","medium","tight"};
   
  for(UInt_t bin=1; bin<=4; bin++){
    h_jetid_->GetXaxis()->SetBinLabel(bin,jetIdBinLabels[bin-1]);    
  }

  h_jetMCflav_ = HistoVars.make<TH1F>(N_+"_MCflav","Jet MC flavour;MC Jet parton pdg ID",22,-0.5,21.5);

}

void ZbbCollections::fillMu(const EventCategory& ec_mu){
   
  Bool_t thecut_= ZbbUtils::ProduceTheCut(ec_mu,N_,"");
  Double_t w = ZbbUtils::getTheWeight(ec_mu,N_);

  if(thecut_){
    h_zllmass_->Fill(ec_mu.bestZcandidate_.mass(),w);  
    h_zllpt_->Fill(ec_mu.bestZcandidate_.pt(),w);
    h_zlleta_->Fill(ec_mu.bestZcandidate_.eta(),w);
    h_zmumueta_->Fill(ec_mu.bestZcandidate_.eta(),w);
    h_zmumumass_->Fill(ec_mu.bestZcandidate_.mass(),w);  
    h_zmumupt_->Fill(ec_mu.bestZcandidate_.pt(),w);  
    h_mu1pt_->Fill(ec_mu.bestZcandidate_.daughter(0)->pt(),w);    
    h_mu2pt_->Fill(ec_mu.bestZcandidate_.daughter(1)->pt(),w);    
    h_mu1eta_->Fill(ec_mu.bestZcandidate_.daughter(0)->eta(),w);   
    h_mu2eta_->Fill(ec_mu.bestZcandidate_.daughter(1)->eta(),w);  
    double deltaRll   = ROOT::Math::VectorUtil::DeltaR(ec_mu.bestZcandidate_.daughter(0)->momentum(),ec_mu.bestZcandidate_.daughter(1)->momentum());
    double deltaPhill = ZbbUtils::myDeltaPhi(ec_mu.bestZcandidate_.daughter(0)->phi(),ec_mu.bestZcandidate_.daughter(1)->phi());
    h_lldeltaR_->Fill(deltaRll,w);
    h_lldeltaPhi_->Fill(deltaPhill,w); 
    h_mumudeltaR_->Fill(deltaRll,w);  
    h_mumudeltaPhi_->Fill(deltaPhill,w); 
    h_zllCosThetaCS_->Fill(ZbbUtils::CostCS(ec_mu.bestZcandidate_),w);   
    h_zllCosThetaStar_->Fill(ZbbUtils::CosThetaStar(ec_mu.bestZcandidate_),w); 
    h_zmumuCosThetaCS_->Fill(ZbbUtils::CostCS(ec_mu.bestZcandidate_),w);      
    h_zmumuCosThetaStar_->Fill(ZbbUtils::CosThetaStar(ec_mu.bestZcandidate_),w); 
  }
  return;  
}

void ZbbCollections::fillEle(const EventCategory& ec_ele){

  //std::cout<<"ZbbCollections::fillEle"<<std::endl;
  Bool_t thecut_ = ZbbUtils::ProduceTheCut(ec_ele,N_,"");
  Double_t w = ZbbUtils::getTheWeight(ec_ele,N_);

  if(thecut_){
    h_zllmass_->Fill(ec_ele.bestZcandidate_.mass(),w);  
    h_zllpt_->Fill(ec_ele.bestZcandidate_.pt(),w);
    h_zeeeta_->Fill(ec_ele.bestZcandidate_.eta(),w);
    h_zeemass_->Fill(ec_ele.bestZcandidate_.mass(),w); 
    h_zeemass_->Fill(ec_ele.bestZcandidate_.mass(),w);  
    h_zeept_->Fill(ec_ele.bestZcandidate_.pt(),w);  
    h_ele1pt_->Fill(ec_ele.bestZcandidate_.daughter(0)->pt(),w);    
    h_ele2pt_->Fill(ec_ele.bestZcandidate_.daughter(1)->pt(),w);    
    h_ele1eta_->Fill(ec_ele.bestZcandidate_.daughter(0)->eta(),w);   
    h_ele2eta_->Fill(ec_ele.bestZcandidate_.daughter(1)->eta(),w); 
    double deltaRll   = ROOT::Math::VectorUtil::DeltaR(ec_ele.bestZcandidate_.daughter(0)->momentum(),ec_ele.bestZcandidate_.daughter(1)->momentum());
    double deltaPhill = ZbbUtils::myDeltaPhi(ec_ele.bestZcandidate_.daughter(0)->phi(),ec_ele.bestZcandidate_.daughter(1)->phi());
    h_lldeltaR_->Fill(deltaRll,w);
    h_lldeltaPhi_->Fill(deltaPhill,w); 
    h_eedeltaR_->Fill(deltaRll,w);  
    h_eedeltaPhi_->Fill(deltaPhill,w);
    h_zllCosThetaCS_->Fill(ZbbUtils::CostCS(ec_ele.bestZcandidate_),w);   
    h_zllCosThetaStar_->Fill(ZbbUtils::CosThetaStar(ec_ele.bestZcandidate_),w); 
    h_zeeCosThetaCS_->Fill(ZbbUtils::CostCS(ec_ele.bestZcandidate_),w);      
    h_zeeCosThetaStar_->Fill(ZbbUtils::CosThetaStar(ec_ele.bestZcandidate_),w); 
  }  
  return;
}

void ZbbCollections::fillJets(const EventCategory& ec_mu,const EventCategory& ec_ele, edm::Handle<edm::View<pat::Jet> > jets, edm::Handle<edm::View<pat::MET> > mets){

  //std::cout<<"ZbbCollections::fillJets"<<std::endl;
  Int_t nJets = 0;

  Bool_t thecut_;   
  EventCategory ec; 
  
  if(ec_mu.isLLTightAndTriggerMatched()){
    thecut_= ZbbUtils::ProduceTheCut(ec_mu,N_,"");
    ec = ec_mu;
  } else if(ec_ele.isLLTightAndTriggerMatched()){
    thecut_ = ZbbUtils::ProduceTheCut(ec_ele,N_,"");
    ec = ec_ele;
  } else {
    return;
  }
  
  Double_t w = ZbbUtils::getTheWeight(ec,N_);

  // vector of good Jets
  std::vector<pat::Jet> theGoodJets;

  if(thecut_){

    // fill the event-weight after selection
    h_eventWeight_->Fill(w);

    for(edm::View<pat::Jet>::const_iterator jetcand=jets->begin(); jetcand!=jets->end(); ++jetcand){
       
      if( ZbbUtils::hasOverlap(jetcand->eta(),jetcand->phi(),ec.bestZcandidate_.daughter(0)->eta(),ec.bestZcandidate_.daughter(0)->phi(),0.5) && 
	  ZbbUtils::hasOverlap(jetcand->eta(),jetcand->phi(),ec.bestZcandidate_.daughter(1)->eta(),ec.bestZcandidate_.daughter(1)->phi(),0.5) ){
	h_jetoverlapLept_->Fill(1.,w);
      } else {
	h_jetoverlapLept_->Fill(0.,w);	
      }

      if ( ZbbUtils::isJetIdOk((*jetcand),"tight")){
	h_jetid_->Fill(3.,w);
      } else if (ZbbUtils::isJetIdOk((*jetcand),"medium")){
	h_jetid_->Fill(2.,w);
      } else if (ZbbUtils::isJetIdOk((*jetcand),"loose")){
	h_jetid_->Fill(1.,w);
      } else{ 
	h_jetid_->Fill(0.,w);
      }

      h_rawjetspt_->Fill(jetcand->pt(),w); 
      h_rawjetseta_->Fill(jetcand->eta(),w);   
	
      if(ZbbUtils::isJetIdOk((*jetcand),"loose") && ZbbUtils::isGoodJet((*jetcand),ec.bestZcandidate_,lCuts_)){
	nJets++;
	
	theGoodJets.push_back(*jetcand);

	h_jetpt_->Fill(jetcand->pt(),w);         
	h_jeteta_->Fill(jetcand->eta(),w);        
	h_jetphi_->Fill(jetcand->phi(),w);       
	h_nhf_->Fill((jetcand->neutralHadronEnergy() + jetcand->HFHadronEnergy() ) / jetcand->energy(),w);           
	h_nef_->Fill(jetcand->neutralEmEnergyFraction(),w);           
	h_nconstituents_->Fill(jetcand->numberOfDaughters(),w); 
	h_chf_->Fill(jetcand->chargedHadronEnergyFraction(),w);          
	h_nch_->Fill(jetcand->chargedHadronEnergyFraction(),w);           
	h_cef_->Fill(jetcand->chargedEmEnergyFraction(),w);           
 
	h_SSVHEdisc_->Fill(jetcand->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"),w);     
	h_SSVHPdisc_->Fill(jetcand->bDiscriminator("simpleSecondaryVertexHighPurBJetTags"),w);     
	h_TCHEdisc_->Fill(jetcand->bDiscriminator("trackCountingHighEffBJetTags"),w);    
        h_TCHPdisc_->Fill(jetcand->bDiscriminator("trackCountingHighPurBJetTags"),w);  
	h_CSVdisc_->Fill(jetcand->bDiscriminator("combinedSecondaryVertexBJetTags"),w);
	h_JPdisc_->Fill(jetcand->bDiscriminator("jetProbabilityBJetTags"),w);
	if ( ec.isMC() ) {
	  h_jetMCflav_->Fill(abs(jetcand->partonFlavour()));
	}
      } 
    } // close the for on jets 

    h_met_->Fill((*mets)[0].et(),w);           
    h_phimet_->Fill((*mets)[0].phi(),w);  
    h_nj_->Fill(nJets); 
 
    //if good jets

    if(nJets>0){
      h_jet1pt_->Fill(theGoodJets[0].pt(),w);       
      h_jet1eta_->Fill(theGoodJets[0].eta(),w);      
      if(nJets>1){
	h_jet2pt_->Fill(theGoodJets[1].pt(),w);       
	h_jet2eta_->Fill(theGoodJets[1].eta(),w);
      }
    }

    for(std::vector<pat::Jet>::const_iterator jet1=theGoodJets.begin(); jet1!= theGoodJets.end(); ++jet1){
      for(std::vector<pat::Jet>::const_iterator jet2=jet1+1; jet2!=theGoodJets.end(); ++jet2){
	if(jet1!=jet2){
	  double deltaRjj   = ROOT::Math::VectorUtil::DeltaR(jet1->momentum(),jet2->momentum());
	  double deltaPhijj = ZbbUtils::myDeltaPhi(jet1->phi(),jet2->phi());
	  h_alljetDeltaR_->Fill(deltaRjj,w);  
	  h_alljetDeltaPhi_->Fill(deltaPhijj,w);  
	}
      }
    }
  } // if right category
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ZbbBasicComponents methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ZbbBasicComponents::book(){
  edm::Service<TFileService> fs;
  std::string DirName = &(*N_);
  TFileDirectory HistoVars = fs->mkdir(DirName);
  
  TH1F::SetDefaultSumw2(kTRUE);

  //jet stuff section
  h_electron_SIP_         = HistoVars.make<TH1F>(N_+"_electron_SIP","electron IP/#sigma_{IP};electron IP/#sigma_{IP}",60,0.,20.);  
  h_electron_iso_         = HistoVars.make<TH1F>(N_+"_electron_iso","electron iso; (sum_{trk} p_{T} + sum_{cal} E_{T})/p_{T} (combRelIso);",50,0.,0.1);
  h_electron_eop_         = HistoVars.make<TH1F>(N_+"_electron_eop","electron eop; electron E/p;",40,0.,1.);
  h_electron_clus_        = HistoVars.make<TH1F>(N_+"_electron_clus","electron clus; electron E_{1x5}/E_{5x5}",40,0.,1.);
  h_electronpair_type_    = HistoVars.make<TH1F>(N_+"_electronpairType","electron pair class; electron pair type",3,-0.5,2.5);
  h_electron_type_        = HistoVars.make<TH1F>(N_+"_electronType","electron class; electron type",2,-0.5,1.5);
  h_electron_eIDs_        = HistoVars.make<TH1F>(N_+"_electron_eIDs","electron eIDS; electron eleID;",4,0.,4.);
  h_electron_eIDcut_      = HistoVars.make<TH1F>(N_+"_Ele_eleID","Electron ID cut",8,-0.5,7.5);
  h_electron_eIDcut_->GetXaxis()->SetBinLabel(1,"Fails");
  h_electron_eIDcut_->GetXaxis()->SetBinLabel(2,"ele ID only");
  h_electron_eIDcut_->GetXaxis()->SetBinLabel(3,"ele iso only");
  h_electron_eIDcut_->GetXaxis()->SetBinLabel(4,"ele ID+iso");
  h_electron_eIDcut_->GetXaxis()->SetBinLabel(5,"conv rej");
  h_electron_eIDcut_->GetXaxis()->SetBinLabel(6,"conv rej+ID");
  h_electron_eIDcut_->GetXaxis()->SetBinLabel(7,"conv rej+iso");
  h_electron_eIDcut_->GetXaxis()->SetBinLabel(8,"all");

  TString electrontypesLabels[2]={"Barrel","Endcap"};
  TString muontypesLabels[4]={"Barrel","Overlap","Encap","CSCOnly"};
  TString leptonpairtypesLabels[3]={"Barrel-Barrel","Barrel-Endcap","Endcap-Endcap"};
  
  for(UInt_t bin=1; bin<=3; bin++){
    if(bin <3)   h_electron_type_->GetXaxis()->SetBinLabel(bin,electrontypesLabels[bin-1]);     
    h_electronpair_type_->GetXaxis()->SetBinLabel(bin,leptonpairtypesLabels[bin-1]);    
  }
  
  h_muon_SIP_             = HistoVars.make<TH1F>(N_+"_muon_SIP","muon IP/#sigma_{IP};muon IP/#sigma_{IP}",60,0.,20.); 
  h_muon_iso_             = HistoVars.make<TH1F>(N_+"_muon_iso","muon iso; (sum_{trk} p_{T} + sum_{cal} E_{T}) / p_{T} (combRelIso)",50,0.,0.1);
  h_muon_chi2_            = HistoVars.make<TH1F>(N_+"_muon_chi2","muon #chi^{2}/ndof; muon #chi^{2}/ndof", 100,0.,10.);
  h_muon_trackerhits_     = HistoVars.make<TH1F>(N_+"_muon_trackerhits","muon Trk hits;tracker hits",40,-0.5,39.5);
  h_muon_pixelhits_       = HistoVars.make<TH1F>(N_+"_muon_pixelhits","muon Pxl hits;pixel hits",10,-0.5,9.5);
  h_muon_muonhits_        = HistoVars.make<TH1F>(N_+"_muon_muonhits","muon Muon hits;muon hits",60,-0.,59.5);
  
  h_muon_muonhitsCSCOnly_ = HistoVars.make<TH1F>(N_+"_muon_muhitsCSCOnly","muon Muon hits (CSC Only);muon hits (|#eta|>2.1)",50,-0.,49.5);
  h_muon_muonhitsBarrel_  = HistoVars.make<TH1F>(N_+"_muon_muhitsBarrel","muon Muon hits (DT Only);muon hits (|#eta|<0.9)",60,-0.,59.5);
  h_muon_muonhitsOverlap_ = HistoVars.make<TH1F>(N_+"_muon_muhitsOverlap","muon Muon hits (DT+CSC+RPC);muon hits (0.9<#eta<1.2)",50,-0.,49.5);
  h_muon_muonhitsEndcaps_ = HistoVars.make<TH1F>(N_+"_muon_muhitsEndcaps","muon Muon hits (CSC+RPC);muon hits (1.2<#eta<2.1)",50,-0.,49.5);
    
  h_muon_numberOfMatches_ = HistoVars.make<TH1F>(N_+"_muon_nOfMatches","muon number of matches;n.of matches",10,-0.5,9.5);
  h_muonpair_type_        = HistoVars.make<TH1F>(N_+"_muonpairType","muon pair class; muon pair type",3,-0.5,2.5);
  h_muon_type_            = HistoVars.make<TH1F>(N_+"_muonType","muon class; muon type",4,-0.5,3.5);

  for(UInt_t bin=1; bin<=4; bin++){
    h_muon_type_->GetXaxis()->SetBinLabel(bin,muontypesLabels[bin-1]);     
    if(bin <4) h_muonpair_type_->GetXaxis()->SetBinLabel(bin,leptonpairtypesLabels[bin-1]);    
  }
}

void ZbbBasicComponents::fillEle(const EventCategory& ec_ele){
  
  //std::cout<<"ZbbBasicComponents::fillEle"<<std::endl;
  const reco::Candidate* lep0 = ec_ele.bestZcandidate_.daughter(0);
  const reco::Candidate* lep1 = ec_ele.bestZcandidate_.daughter(1);
  
  Bool_t thecut_= ZbbUtils::ProduceTheCut(ec_ele,N_,"");
  Double_t w = ZbbUtils::getTheWeight(ec_ele,N_);

  if(thecut_){

    const pat::Electron* ele0 = dynamic_cast<const pat::Electron*>(&(*lep0->masterClone()));
    const pat::Electron* ele1 = dynamic_cast<const pat::Electron*>(&(*lep1->masterClone()));

    if(ele0->isEB()&&ele1->isEB()){
      h_electronpair_type_->Fill(0.,w);
    } else if( (ele0->isEB()&&ele1->isEE())||(ele0->isEE()&&ele1->isEB()) ){
      h_electronpair_type_->Fill(1.,w);
    } else {
      h_electronpair_type_->Fill(2.,w);
    }

    std::vector< const pat::Electron* > electrons;
    electrons.push_back(ele0);
    electrons.push_back(ele1);

    for(UInt_t i=0; i<electrons.size(); i++){

      if(electrons[i]->isEB()){
	h_electron_type_->Fill(0.,w);
      } else if (electrons[i]->isEE()){
	h_electron_type_->Fill(1.,w);
      }

      h_electron_SIP_->Fill(fabs(electrons[i]->dB(pat::Electron::PV3D))/electrons[i]->edB(pat::Electron::PV3D),w); 
      h_electron_iso_->Fill((electrons[i]->trackIso()+electrons[i]->caloIso())/electrons[i]->pt(),w);        
      h_electron_eop_->Fill(electrons[i]->eSeedClusterOverP(),w);        
      h_electron_clus_->Fill(electrons[i]->e1x5()/electrons[i]->e5x5(),w);       

      //Specific electron histograms
      if( electrons[i]->electronID("eidVBTFRel95")==0 )h_electron_eIDcut_->Fill(0);
      if( electrons[i]->electronID("eidVBTFRel95")==1 )h_electron_eIDcut_->Fill(1);
      if( electrons[i]->electronID("eidVBTFRel95")==2 )h_electron_eIDcut_->Fill(2);
      if( electrons[i]->electronID("eidVBTFRel95")==3 )h_electron_eIDcut_->Fill(3);
      if( electrons[i]->electronID("eidVBTFRel95")==4 )h_electron_eIDcut_->Fill(4);
      if( electrons[i]->electronID("eidVBTFRel95")==5 )h_electron_eIDcut_->Fill(5);
      if( electrons[i]->electronID("eidVBTFRel95")==6 )h_electron_eIDcut_->Fill(6);
      if( electrons[i]->electronID("eidVBTFRel95")==7 )h_electron_eIDcut_->Fill(7);
      
      if( electrons[i]->electronID("eidVBTFCom80") > 0.5 ) h_electron_eIDs_->Fill(0.,w);
      if( electrons[i]->electronID("eidVBTFCom95") > 0.5 ) h_electron_eIDs_->Fill(1.,w);
      if( electrons[i]->electronID("eidVBTFRel80") > 0.5 ) h_electron_eIDs_->Fill(2.,w);
      if( electrons[i]->electronID("eidVBTFRel95") > 0.5 ) h_electron_eIDs_->Fill(3.,w);
    }
  }
  return;
}

void ZbbBasicComponents::fillMu(const EventCategory& ec_mu){
  
  //std::cout<<"ZbbBasicComponents::fillMu"<<std::endl;
  const reco::Candidate* lep0 = ec_mu.bestZcandidate_.daughter(0);
  const reco::Candidate* lep1 = ec_mu.bestZcandidate_.daughter(1);
  
  Bool_t thecut_ = ZbbUtils::ProduceTheCut(ec_mu,N_,"");
  Double_t w = ZbbUtils::getTheWeight(ec_mu,N_);
 
  if(thecut_){
    
    const pat::Muon* muon0 = dynamic_cast<const pat::Muon*>(&(*lep0->masterClone()));
    const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(&(*lep1->masterClone())); 
    
    if(fabs(muon0->eta())<0.9 && fabs(muon1->eta())<0.9 ){
      h_muonpair_type_->Fill(0.,w);
    } else if( (fabs(muon0->eta())>0.9 && fabs(muon1->eta())<0.9) || (fabs(muon0->eta())<0.9 && fabs(muon1->eta())>0.9)){
      h_muonpair_type_->Fill(1.,w);
    } else {
      h_muonpair_type_->Fill(2.,w);
    }

    std::vector< const pat::Muon*> muons;
    muons.push_back(muon0);
    muons.push_back(muon1);
    
    for(UInt_t i=0; i<muons.size(); i++){
      
      h_muon_iso_->Fill((muons[i]->trackIso() + muons[i]->caloIso())/muons[i]->pt(),w);             
      h_muon_SIP_->Fill(fabs(muons[i]->dB(pat::Muon::PV3D))/muons[i]->edB(pat::Muon::PV3D),w);   
          
      if(muons[i]->isGlobalMuon()==true || muons[i]->isTrackerMuon()==true){
	h_muon_chi2_->Fill(muons[i]->globalTrack()->normalizedChi2(),w);            
	h_muon_trackerhits_->Fill(muons[i]->innerTrack()->hitPattern().numberOfValidTrackerHits(),w);     
	h_muon_pixelhits_->Fill(muons[i]->innerTrack()->hitPattern().numberOfValidPixelHits(),w);       
	h_muon_muonhits_->Fill(muons[i]->globalTrack()->hitPattern().numberOfValidMuonHits(),w);       
	
	if(fabs(muons[i]->eta())<0.9){
	  h_muon_type_->Fill(0.,w);
	  h_muon_muonhitsBarrel_->Fill(muons[i]->globalTrack()->hitPattern().numberOfValidMuonHits(),w);
	} else if( muons[i]->eta() >0.9 && muons[i]->eta()<1.2){
	  h_muon_type_->Fill(1.,w);
	  h_muon_muonhitsOverlap_->Fill(muons[i]->globalTrack()->hitPattern().numberOfValidMuonHits(),w);
	} else if( muons[i]->eta() >1.2 && muons[i]->eta()<2.4) {
	  h_muon_type_->Fill(2.,w);
	  h_muon_muonhitsEndcaps_->Fill(muons[i]->globalTrack()->hitPattern().numberOfValidMuonHits(),w);
	} else {
	  h_muon_type_->Fill(3.,w);
	  h_muon_muonhitsCSCOnly_->Fill(muons[i]->globalTrack()->hitPattern().numberOfValidMuonHits(),w);
	}
	h_muon_numberOfMatches_->Fill(muons[i]->numberOfMatches(),w); 
      }
    }
  }
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ZbbMCinfo methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ZbbMCinfo::book(){
  edm::Service<TFileService> fs;
  std::string DirName = &(*N_);
  TFileDirectory HistoVars = fs->mkdir(DirName);

  Int_t nbin_pt=70;
  Double_t pt_min= 0;
  Double_t pt_max=140;
  Int_t nbin_eta=80;
  Double_t eta_min=-4;
  Double_t eta_max =4;

  TH1F::SetDefaultSumw2(kTRUE);

  // Z plots
  h_GENP_pt_Z_   = HistoVars.make<TH1F>(N_+"_GENP_pt_Z","p_{T} of Z boson (genParticles); GEN p_{T} of the Z (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_eta_Z_  = HistoVars.make<TH1F>(N_+"_GENP_eta_Z","#eta of Z boson (genParticles); #eta of the Z (GEN)",nbin_eta,eta_min,eta_max); 
  h_GENP_mass_Z_ = HistoVars.make<TH1F>(N_+"_GENP_mass_Z","mass of Z boson (genParticles); GEN M_{Z} (GeV)",120,60.,120.); 

  h_LHE_pt_Z_   = HistoVars.make<TH1F>(N_+"_LHE_pt_Z","p_{T} of Z boson (LHE event); LHE p_{T} of the Z (GeV)",nbin_pt,pt_min,pt_max);  
  h_LHE_eta_Z_  = HistoVars.make<TH1F>(N_+"_LHE_eta_Z","#eta of Z boson (LHE event); #eta of the Z (LHE)",nbin_eta,eta_min,eta_max); 
  h_LHE_mass_Z_ = HistoVars.make<TH1F>(N_+"_LHE_mass_Z","mass of Z boson (LHE event); LHE M_{Z} (GeV)",120,60.,120.); 

  h_GENP_pt_ZfromLept_   = HistoVars.make<TH1F>(N_+"_GENP_pt_ZfromLept","p_{T} of ll (genParticles); GEN p_{T}(l^{+}l^{-}) (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_eta_ZfromLept_  = HistoVars.make<TH1F>(N_+"_GENP_eta_ZfromLept","#eta of ll (genParticles); #eta(l^{+}l^{-}) (GEN)",nbin_eta,eta_min,eta_max); 
  h_GENP_mass_ZfromLept_ = HistoVars.make<TH1F>(N_+"_GENP_mass_ZfromLept","mass of ll (genParticles); GEN M_{l^{+}l^{-}} (GeV)",120,60.,120.); 

  h_GENP_pt_ZfromMuons_   = HistoVars.make<TH1F>(N_+"_GENP_pt_ZfromMuons","p_{T} of #mu#mu (genParticles); GEN p_{T}(#mu^{+}#mu^{-}) (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_eta_ZfromMuons_  = HistoVars.make<TH1F>(N_+"_GENP_eta_ZfromMuons","#eta of #mu#mu (genParticles); #eta(#mu^{+}#mu^{-}) (GEN)",nbin_eta,eta_min,eta_max); 
  h_GENP_mass_ZfromMuons_ = HistoVars.make<TH1F>(N_+"_GENP_mass_ZfromMuons","mass of #mu#mu (genParticles); GEN M_{#mu^{+}#mu^{-}} (GeV)",120,60.,120.);
  
  h_GENP_pt_ZfromElectrons_   = HistoVars.make<TH1F>(N_+"_GENP_pt_ZfromElectrons","p_{T} of #it{ee} (genParticles); GEN p_{T}(e ^{+}e ^{-}) (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_eta_ZfromElectrons_  = HistoVars.make<TH1F>(N_+"_GENP_eta_ZfromElectrons","#eta of #it{ee} (genParticles); #eta(e^{+}e^{-}) (GEN)",nbin_eta,eta_min,eta_max); 
  h_GENP_mass_ZfromElectrons_ = HistoVars.make<TH1F>(N_+"_GENP_mass_ZfromElectrons","mass of #it{ee} (genParticles); GEN M_{e^{+}e^{-}} (GeV)",120,60.,120.);

  // leading b plots
  h_ncGENP_       = HistoVars.make<TH1F>(N_+"_ncGENP","number of decaying c partons in the event (genParticles); GEN n^{lch}_{c}",6,-0.5,5.5); 
  h_nbGENP_       = HistoVars.make<TH1F>(N_+"_nbGENP","number of decaying b partons in the event (genParticles); GEN n^{lbh}_{b}",6,-0.5,5.5); 
  h_GENP_pt_b1_   = HistoVars.make<TH1F>(N_+"_GENP_pt_b1","p_{T} leading b in the event (genParticles); p_{T} leading b (GeV)",nbin_pt,pt_min,pt_max);  
  h_GENP_eta_b1_  = HistoVars.make<TH1F>(N_+"_GENP_eta_b1","#eta leading b in the event (genParticles); #eta leading b",nbin_eta,eta_min,eta_max); 
  h_GENP_mass_b1_ = HistoVars.make<TH1F>(N_+"_GENP_mass_b1","mass leading b in the event (genParticles); m_{b} (GEN) GeV",100,-0.5,5.5); 
  h_GENP_recoMatched_pt_b1_  = HistoVars.make<TH1F>(N_+"_GENP_recoMatched_pt_b1","p_{T} leading b in the event (genParticles - reco matched);GEN p^{b}_{T} (GeV) (after reco match)",nbin_pt,pt_min,pt_max); 
  h_GENP_recoMatched_eta_b1_ = HistoVars.make<TH1F>(N_+"_GENP_recoMatched_eta_b1","#eta leading b in the event (genParticles - reco matched);GEN #eta^{b} (after reco match)",nbin_eta,eta_min,eta_max); 
  h_GENP_recoBtagged_pt_b1_  = HistoVars.make<TH1F>(N_+"_GENP_recoBtagged_pt_b1","p_{T} leading b in the event (genParticles - b-tag matched);GEN p^{b}_{T} (GeV) (after b-tag match)",nbin_pt,pt_min,pt_max);
  h_GENP_recoBtagged_eta_b1_ = HistoVars.make<TH1F>(N_+"_GENP_recoBtagged_eta_b1","#eta leading b in the event (genParticles - b-tag matched);GEN #eta^{b} (after b-tag match)",nbin_eta,eta_min,eta_max); 

  // for matching studies
  h_GENPDeltaR_b1LBH_RECOJet_ = HistoVars.make<TH1F>(N_+"GENPDeltaR_bLBH_RECOJet","#DeltaR(b,reco::Jet); #DeltaR(b_{GEN},j)",100,0.,5.); 

  // histos for eff b/c studies
  float bTag_pTbins[12]={0, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 240};
  h_GENP_pt_b1_etaBrl_=HistoVars.make<TH1F>(N_+"_GENP_pt_b1_etaBrl","p_{T} leading b in the event (genParticles) barrel; p_{T} leading b (GeV)",11,bTag_pTbins);
  h_GENP_pt_b1_etaFwd_=HistoVars.make<TH1F>(N_+"_GENP_pt_b1_etaFwd","p_{T} leading b in the event (genParticles) fwd; p_{T} leading b (GeV)",11,bTag_pTbins);
  h_GENP_pt_c1_etaBrl_=HistoVars.make<TH1F>(N_+"_GENP_pt_c1_etaBrl","p_{T} leading c in the event (genParticles) barrel; p_{T} leading c (GeV)",11,bTag_pTbins);
  h_GENP_pt_c1_etaFwd_=HistoVars.make<TH1F>(N_+"_GENP_pt_c1_etaFwd","p_{T} leading c in the event (genParticles) fwd; p_{T} leading c (GeV)",11,bTag_pTbins);
  h_RECO_pt_b1_etaBrl_=HistoVars.make<TH1F>(N_+"_RECO_pt_b1_etaBrl","p_{T} leading b in the event (recoJet) barrel; p_{T} leading b (GeV)",11,bTag_pTbins);
  h_RECO_pt_b1_etaFwd_=HistoVars.make<TH1F>(N_+"_RECO_pt_b1_etaFwd","p_{T} leading b in the event (recoJet) fwd; p_{T} leading b (GeV)",11,bTag_pTbins);
  h_RECO_pt_c1_etaBrl_=HistoVars.make<TH1F>(N_+"_RECO_pt_c1_etaBrl","p_{T} leading c in the event (recoJet) barrel; p_{T} leading c (GeV)",11,bTag_pTbins);
  h_RECO_pt_c1_etaFwd_=HistoVars.make<TH1F>(N_+"_RECO_pt_c1_etaFwd","p_{T} leading c in the event (recoJet) fwd; p_{T} leading c (GeV)",11,bTag_pTbins);

  h_nbLHE_       = HistoVars.make<TH1F>(N_+"_nbLHE","number of LHE b partons in the event (genParticles); LHE n_{b}",6,-0.5,5.5); 
  h_LHE_pt_b1_   = HistoVars.make<TH1F>(N_+"_LHE_pt_b1","pt leading b in the event (LHE event); p_{T} leading b (GeV)",nbin_pt,pt_min,pt_max);  
  h_LHE_eta_b1_  = HistoVars.make<TH1F>(N_+"_LHE_eta_b1","#eta leading b in the event (LHE event); #eta leading b",nbin_eta,eta_min,eta_max); 
  h_LHE_mass_b1_ = HistoVars.make<TH1F>(N_+"_LHE_mass_b1","mass leading b in the event (genParticles); m_{b} GeV",100,-0.5,5.5); 
  h_LHE_recoMatched_pt_b1_  = HistoVars.make<TH1F>(N_+"_LHE_recoMatched_pt_b1","p_{T} leading b in the event (LHE event - reco matched);LHE p^{b}_{T} (GeV) (after reco match)",nbin_pt,pt_min,pt_max);  
  h_LHE_recoMatched_eta_b1_ = HistoVars.make<TH1F>(N_+"_LHE_recoMatched_eta_b1","#eta leading b in the event (LHE event - reco matched);LHE #eta^{b} (after reco match)",nbin_eta,eta_min,eta_max); 
  h_LHE_recoBtagged_pt_b1_  = HistoVars.make<TH1F>(N_+"_LHE_recoBtagged_pt_b1","p_{T} leading b in the event (LHE event - b-tag matched);LHE p^{b}_{T} (GeV) (after b-tag match)",nbin_pt,pt_min,pt_max);  
  h_LHE_recoBtagged_eta_b1_ = HistoVars.make<TH1F>(N_+"_LHE_recoBtagged_eta_b1","#eta leading b in the event (LHE event - b-tag matched);LHE #eta^{b} (after b-tag match)",nbin_eta,eta_min,eta_max); 
  
  // correlations
  h2_LHEpt_vs_GENpt     = HistoVars.make<TH2F>(N_+"_LHEpt_Vs_GENpt","p_{T} leading b LHE vs GEN;GEN p^{b}_{T} (GeV);LHE p^{b}_{T} (GeV)",nbin_pt,pt_min,pt_max,nbin_pt,pt_min,pt_max); ;   
  h2_LHEpt_vs_RECOJetpt = HistoVars.make<TH2F>(N_+"_LHEpt_Vs_RECOJetpt","p_{T} leading b LHE vs p_{T} jet;jet p_{T} (GeV); LHE p^{b}_{T} (GeV)",nbin_pt,pt_min,pt_max,nbin_pt,pt_min,pt_max); 
  h2_GENpt_vs_RECOJetpt = HistoVars.make<TH2F>(N_+"_GENpt_Vs_RECOJetpt","p_{T} leading b GEN vs p_{T} jet;jet p_{T} (GeV); GEN p^{b}_{T} (GeV)",nbin_pt,pt_min,pt_max,nbin_pt,pt_min,pt_max); 
}

void ZbbMCinfo::fill(const EventCategory& ec_mu, const EventCategory& ec_ele,std::string bTagAlgoWP, edm::Handle<LHEEventProduct> lh_evt,edm::Handle<reco::GenParticleCollection> genParticlesCollection, const reco::GenJetCollection& genJets ,edm::Handle<edm::View<pat::Jet> > recoJets) {

  //std::cout<<"ZbbMCinfo::fill"<<std::endl;

  //**********************************
  // LHE collection
  //**********************************

  Bool_t hasLHEbPartonJetMatched_ = false;
  Bool_t hasLHEbPartonbJetMatched_ = false;
  
  // initialize pT and eta of the leading b/Z
  math::XYZTLorentzVectorD p4LHEmax(0,0,0,0); 
  math::XYZTLorentzVectorD p4ZLHE(0,0,0,0);
  int nbLHE=0;
  int nZLHE=0;
  Double_t jetPtLHEMatched(0.);

  if(lh_evt.isValid ()){
    
    const lhef::HEPEUP hepeup_ = lh_evt->hepeup();
    const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP;
    const std::vector<int> idup_ = hepeup_.IDUP;
    const std::vector<int> istup_ = hepeup_.ISTUP;
    
    // loop on the particles in the event
    int nup = hepeup_.NUP;
    for ( int iup=0; iup<nup; iup++){
      if(idup_[iup]==23 && istup_[iup]==2){
	nZLHE++;
	math::XYZTLorentzVectorD p4ZLHEtemp((pup_[iup])[0],(pup_[iup])[1],(pup_[iup])[2],(pup_[iup])[3]); 
	p4ZLHE = p4ZLHEtemp;
      }
      
      if ( idup_[iup]==5 || idup_[iup]==-5 ){
	if ( istup_[iup]==1) {
	  nbLHE++;  
	  // dump the 4-momentum into a LorentzVector
	  math::XYZTLorentzVectorD p4LHE((pup_[iup])[0],(pup_[iup])[1],(pup_[iup])[2],(pup_[iup])[3]);      	  
	  if ( p4LHE.pt() > p4LHEmax.pt() ){ 
	    p4LHEmax=p4LHE;
	    hasLHEbPartonJetMatched_ = false;
	    hasLHEbPartonbJetMatched_ = false;
	    for(edm::View<pat::Jet>::const_iterator jet=recoJets->begin(); jet!=recoJets->end(); ++jet){
	      if(ROOT::Math::VectorUtil::DeltaR(jet->momentum(),p4LHEmax.Vect())<0.3 && jet->pt()>lCuts_.bjetPtMin_&& fabs(jet->eta())<lCuts_.bjetEtaMax_){
		hasLHEbPartonJetMatched_ = true;
		jetPtLHEMatched=jet->pt();
		if (ZbbUtils::isBJet(*jet,bTagAlgoWP)){
		  hasLHEbPartonbJetMatched_ = true;
		} // if the the good jet is b-tagged
	      } // if jet is in the acceptance
	    } // loop over reco::Jets
	  } // exceeds ptmax
	} // end if status = 1
      } // end if b or bbar
    } // end loop on particles
  } // if LHE event handle is valid

  //**********************************
  // genParticle collection
  //********************************** 
  
  Bool_t isStatus1_,isStatus2_,isStatus3_,isZ_,isb_,isc_,hasHdaughter_,hasLdaughter_;
  
  Bool_t hasbPartonJetMatched_ = false;
  Bool_t hasbPartonbJetMatched_ = false;
  Bool_t isbPartonGenJetMatched_= false;

  Bool_t hascPartonJetMatched_ = false; 
  Bool_t hascPartonbJetMatched_ = false; 
  Bool_t iscPartonGenJetMatched_= false;

  // initialize pT and eta of the leading b/Z
  std::vector<reco::GenParticle> theGenMuons;
  std::vector<reco::GenParticle> theGenElectrons;
  
  math::XYZTLorentzVectorD p4GENmax_c(0,0,0,0);    
  math::XYZTLorentzVectorD p4GENmax_b(0,0,0,0);    
  math::XYZTLorentzVectorD p4ZGEN(0,0,0,0); 
  math::XYZTLorentzVectorD p4ZGENFromLept(0,0,0,0); 
  int nbGEN=0;
  int ncGEN=0;
  int nZGEN=0;
  Double_t jetPtbGENMatched(0.);
  Double_t jetPtcGENMatched(0.);
  Double_t jetEtabGENMatched(0.);
  Double_t jetEtacGENMatched(0.);
  // loop over GEN particles

  Double_t minDeltaRbRecoJet(999.);	
  for(reco::GenParticleCollection::const_iterator genp = genParticlesCollection->begin();genp != genParticlesCollection->end(); ++ genp ) {
       
    //Set boolean to false

    isStatus1_=false;
    isStatus2_=false;
    isStatus3_=false;
    isb_=false;
    isc_=false;
    isZ_=false;
    hasHdaughter_=false;
    hasLdaughter_=false;
    
    if (genp->status()==1)isStatus1_=true;
    if (genp->status()==2)isStatus2_=true;
    if (genp->status()==3)isStatus3_=true;
    
    if (genp->pdgId()==23) isZ_=true;
    if (fabs(genp->pdgId())==4)isc_=true;
    if (fabs(genp->pdgId())==5)isb_=true;

    if( (fabs(genp->pdgId()) == 11 && isStatus1_==true) && (ZbbUtils::getParentCode(&(*genp))==23) ){
      theGenElectrons.push_back(*genp);
    } else if( (fabs(genp->pdgId()) == 13 && isStatus1_==true) && (ZbbUtils::getParentCode(&(*genp))==23) ){
      theGenMuons.push_back(*genp);
    } else if( (fabs(genp->pdgId()) == 11 || fabs(genp->pdgId()) == 13) && isStatus1_==true && (ZbbUtils::getParentCode(&(*genp))!=23) ) {
      //std::cout<<"Ptl: "<<fabs(genp->pdgId())<<" parent code:"<<ZbbUtils::getParentCode(&(*genp))<<std::endl;
    }
 
    // loops on daughter to check "last before hadronization"
    for (size_t i=0;i<genp->numberOfDaughters();i++){
      if (fabs(genp->daughter(i)->pdgId()) == 91 || fabs(genp->daughter(i)->pdgId()) == 92) hasHdaughter_=true;
      if (fabs(genp->daughter(i)->pdgId()) == 11 || fabs(genp->daughter(i)->pdgId()) == 13) hasLdaughter_=true;
    }
    
    // Z part
    if(isZ_ && isStatus2_){
      nZGEN++;
      math::XYZTLorentzVectorD p4ZGENtemp(genp->px(),genp->py(),genp->pz(),genp->energy());    
      p4ZGEN=p4ZGENtemp;
      if(theGenElectrons.size()!=2 && theGenMuons.size()!=2 ){
	//std::cout<<"Ptl:" <<genp->pdgId();
	for (size_t i=0;i<genp->numberOfDaughters();i++){
	  //std::cout<< genp->daughter(i)->pdgId();
	}
	//std::cout<<std::endl;
      }
    }

    //++++++++++++++++++++++++++++++++++++++ the c-parton +++++++++++++++++++++++++++++++++++++++++++

    if (isc_){
      if ((isStatus2_ || isStatus3_) && hasHdaughter_){
	ncGEN++;
	math::XYZTLorentzVectorD p4GEN(genp->px(),genp->py(),genp->pz(),genp->energy());   
	if (p4GEN.pt()>p4GENmax_c.pt()){
	  hascPartonJetMatched_ = false;
	  hascPartonbJetMatched_ = false;
	  p4GENmax_c=p4GEN;
	  
	  iscPartonGenJetMatched_= false;
	  std::pair<int,double> GenGenAssociation = std::make_pair(-1,9999.);
	  double minDeltaRGenGen(9999.);
	  int i(0);

	  // first match the c-parton to the genJet
	  for(std::vector<reco::GenJet>::const_iterator genjet_it=genJets.begin(); genjet_it!=genJets.end(); ++genjet_it){ 
	    if(ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4GENmax_c.Vect())<0.5 ) { 
	      iscPartonGenJetMatched_=true;
	      if(ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4GENmax_c.Vect())< minDeltaRGenGen){
		minDeltaRGenGen = ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4GENmax_c.Vect());
		GenGenAssociation.first  = i;
		GenGenAssociation.second = minDeltaRGenGen;
	      } // end if minimum distance gen jet-parton 
	    } // ends if parton matched to gen jet
	    i++;
	  }// ends loop on gen jets
	  
	  if(!iscPartonGenJetMatched_ || GenGenAssociation.first == -1 ) continue;
	
	  // then if the c is matched to the gen match the reco
	  for(edm::View<pat::Jet>::const_iterator recojet_it=recoJets->begin(); recojet_it!=recoJets->end(); ++recojet_it){  
	    if(ROOT::Math::VectorUtil::DeltaR(genJets.at(GenGenAssociation.first).momentum(),recojet_it->momentum())<0.5 
	       && (ZbbUtils::isJetIdOk((*recojet_it),"loose") && recojet_it->pt()>lCuts_.bjetPtMin_&& TMath::Abs(recojet_it->eta())<lCuts_.bjetEtaMax_) ) {
	      hascPartonJetMatched_ = true;
	      jetPtcGENMatched=recojet_it->pt();   
	      jetEtacGENMatched=recojet_it->eta();
	      if (ZbbUtils::isBJet(*recojet_it,bTagAlgoWP)){
		hascPartonbJetMatched_ = true;
	      } // if the the good jet is b-tagged
	    } // end if distance gen jet-reco::Jet < 0.5
	  } //ends loop on reco::Jets  
	} // if pt exceeds ptmax 
      } // controls status == ( 2 || 3 )
    } // control over b-partons

    //++++++++++++++++++++++++++++++++++++++ the b-parton +++++++++++++++++++++++++++++++++++++++++++

    // setting the variables
    if (isb_){
      if ((isStatus2_ || isStatus3_) && hasHdaughter_){
	nbGEN++;
	math::XYZTLorentzVectorD p4GEN(genp->px(),genp->py(),genp->pz(),genp->energy());
	if (p4GEN.pt()>p4GENmax_b.pt()){
	  hasbPartonJetMatched_ = false;
	  hasbPartonbJetMatched_ = false;
	  p4GENmax_b=p4GEN;
	  
	  isbPartonGenJetMatched_= false;
	  std::pair<int,double> GenGenAssociation = std::make_pair(-1,9999.);
	  double minDeltaRGenGen(9999.);
	  int j(0);
	  
	  // first match the b-parton to the genJet
	  for(std::vector<reco::GenJet>::const_iterator genjet_it=genJets.begin(); genjet_it!=genJets.end(); ++genjet_it){ 
	    if(ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4GENmax_b.Vect())<0.5 ) { 
	      isbPartonGenJetMatched_=true;
	      if(ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4GENmax_b.Vect())< minDeltaRGenGen){
		minDeltaRGenGen = ROOT::Math::VectorUtil::DeltaR(genjet_it->momentum(),p4GENmax_b.Vect());
		GenGenAssociation.first  = j;
		GenGenAssociation.second = minDeltaRGenGen;
	      } // end if minimum distance gen jet-parton 
	    } // ends if parton matched to gen jet
	    j++;
	  }// ends loop on gen jets
	  
	  if(!isbPartonGenJetMatched_ || GenGenAssociation.first == -1 ) continue;
	
	  // then if the c is matched to the gen match the reco
	  for(edm::View<pat::Jet>::const_iterator recojet_it=recoJets->begin(); recojet_it!=recoJets->end(); ++recojet_it){ 
	    if(ROOT::Math::VectorUtil::DeltaR(recojet_it->momentum(),p4GEN.Vect())<minDeltaRbRecoJet){
	      minDeltaRbRecoJet=ROOT::Math::VectorUtil::DeltaR(recojet_it->momentum(),p4GENmax_b.Vect());
	    }  
	    if(ROOT::Math::VectorUtil::DeltaR(genJets.at(GenGenAssociation.first).momentum(),recojet_it->momentum())<0.5 
	       && (ZbbUtils::isJetIdOk((*recojet_it),"loose") && recojet_it->pt()>lCuts_.bjetPtMin_ && TMath::Abs(recojet_it->eta())<lCuts_.bjetEtaMax_) ) { 
	      hasbPartonJetMatched_ = true;
	      jetPtbGENMatched=recojet_it->pt();
	      jetEtabGENMatched=recojet_it->eta();
	      if (ZbbUtils::isBJet(*recojet_it,bTagAlgoWP)){
		hasbPartonbJetMatched_ = true;
	      } // if the the good jet is b-tagged
	    } // end if distance gen jet-reco::Jet < 0.5
	  } //ends loop on reco::Jets  
	} // if pt exceeds ptmax 
      } // controls status == ( 2 || 3 )
    } // control over b-partons
  } // loop on genParticles

  h_GENPDeltaR_b1LBH_RECOJet_->Fill(minDeltaRbRecoJet); 
  
  //-- filling the histos

  // produces the cut
  Bool_t thecut_;   
  EventCategory ec;   

  if(ec_mu.isLLTightAndTriggerMatched()){
    thecut_= ZbbUtils::ProduceTheCut(ec_mu,N_,"");
    ec = ec_mu;
  } else {
    thecut_= ZbbUtils::ProduceTheCut(ec_ele,N_,"");
    ec = ec_ele;
  }

  Double_t w =ZbbUtils::getTheWeight(ec,N_);

  if(thecut_){ 
    
    //*******************************
    // unweighted histos so far  
    //*******************************

    if(nZLHE>0){
      h_LHE_pt_Z_->Fill(p4ZLHE.pt());   
      h_LHE_eta_Z_->Fill(p4ZLHE.eta());
      h_LHE_mass_Z_->Fill(p4ZLHE.mass());
    }

    if(nZGEN>0){
      h_GENP_pt_Z_->Fill(p4ZGEN.pt());     
      h_GENP_eta_Z_->Fill(p4ZGEN.eta()); 
      h_GENP_mass_Z_->Fill(p4ZGEN.mass());
    }

    if(theGenElectrons.size()>1){
      for (size_t i=0;i<theGenElectrons.size();i++){
	math::XYZTLorentzVectorD p4ZGENFromLepttemp(theGenElectrons[i].px(),theGenElectrons[i].py(),theGenElectrons[i].pz(),theGenElectrons[i].energy());  
	p4ZGENFromLept += p4ZGENFromLepttemp;
      }
      
      h_GENP_pt_ZfromLept_->Fill(p4ZGENFromLept.pt());  
      h_GENP_eta_ZfromLept_->Fill(p4ZGENFromLept.eta()); 
      h_GENP_mass_ZfromLept_->Fill(p4ZGENFromLept.mass());
      h_GENP_pt_ZfromElectrons_->Fill(p4ZGENFromLept.pt());  
      h_GENP_eta_ZfromElectrons_->Fill(p4ZGENFromLept.eta()); 
      h_GENP_mass_ZfromElectrons_->Fill(p4ZGENFromLept.mass());
      
    } else if(theGenMuons.size()>1){
      for (size_t i=0;i<theGenMuons.size();i++){
	math::XYZTLorentzVectorD p4ZGENFromLepttemp(theGenMuons[i].px(),theGenMuons[i].py(),theGenMuons[i].pz(),theGenMuons[i].energy());  
	p4ZGENFromLept += p4ZGENFromLepttemp;
      }
      
      h_GENP_pt_ZfromLept_->Fill(p4ZGENFromLept.pt());  
      h_GENP_eta_ZfromLept_->Fill(p4ZGENFromLept.eta()); 
      h_GENP_mass_ZfromLept_->Fill(p4ZGENFromLept.mass());
      h_GENP_pt_ZfromMuons_->Fill(p4ZGENFromLept.pt());  
      h_GENP_eta_ZfromMuons_->Fill(p4ZGENFromLept.eta()); 
      h_GENP_mass_ZfromMuons_->Fill(p4ZGENFromLept.mass());
    }
    
    h_nbLHE_->Fill(nbLHE);  
    h_ncGENP_->Fill(ncGEN);
    h_nbGENP_->Fill(nbGEN);

      
    h2_LHEpt_vs_GENpt->Fill(p4GENmax_b.pt(),p4LHEmax.pt());

    // all LHE partons
    if (nbLHE>0){
      
      h_LHE_pt_b1_->Fill(p4LHEmax.pt());    
      h_LHE_eta_b1_->Fill(p4LHEmax.eta());
      h_LHE_mass_b1_->Fill(p4LHEmax.mass());

      // if parton is associated with a reconstructed jet
      if(hasLHEbPartonJetMatched_){
	h_LHE_recoMatched_pt_b1_->Fill(p4LHEmax.pt());    
	h_LHE_recoMatched_eta_b1_->Fill(p4LHEmax.eta());
	
	h2_LHEpt_vs_RECOJetpt->Fill(jetPtLHEMatched,p4LHEmax.pt());

	// if parton is associated with a reconstructed b-jet
	if(hasLHEbPartonbJetMatched_){
	  h_LHE_recoBtagged_pt_b1_->Fill(p4LHEmax.pt());     
	  h_LHE_recoBtagged_eta_b1_->Fill(p4LHEmax.eta()); 
	} 
      }
    }
    

    // binned pT spectra for eff b/c studies
    if (ncGEN>0){
      // GEN level
      if ( fabs(p4GENmax_c.eta())<1.2 ) {
	h_GENP_pt_c1_etaBrl_->Fill(p4GENmax_c.pt());
      } else if ( fabs(p4GENmax_c.eta())>=1.2  && fabs(p4GENmax_c.eta())<lCuts_.bjetEtaMax_ ) {
	h_GENP_pt_c1_etaFwd_->Fill(p4GENmax_c.pt());
      }
      // RECO jet level
      if ( hascPartonJetMatched_ ){
	if ( fabs(jetEtacGENMatched)<1.2 ) {
	  h_RECO_pt_c1_etaBrl_->Fill(jetPtcGENMatched,w);
	} else if ( fabs(jetEtacGENMatched)>=1.2  && fabs(jetEtacGENMatched)<lCuts_.bjetEtaMax_ ) {
	  h_RECO_pt_c1_etaFwd_->Fill(jetPtcGENMatched,w);
	}
      }//end of RECO jet level

    }
    if (nbGEN>0){
      // GEN level
      if ( fabs(p4GENmax_b.eta())<1.2 ) {
	h_GENP_pt_b1_etaBrl_->Fill(p4GENmax_b.pt());
      } else if ( fabs(p4GENmax_b.eta())>=1.2  && fabs(p4GENmax_b.eta())<lCuts_.bjetEtaMax_ ) {
	h_GENP_pt_b1_etaFwd_->Fill(p4GENmax_b.pt());
      }

      // RECO jet level
      if ( hasbPartonJetMatched_ ){
	if ( fabs(jetEtabGENMatched)<1.2 ) {
	  h_RECO_pt_b1_etaBrl_->Fill(jetPtbGENMatched,w);
	} else if ( fabs(jetEtabGENMatched)>=1.2  && fabs(jetEtabGENMatched)<lCuts_.bjetEtaMax_ ) {
	  h_RECO_pt_b1_etaFwd_->Fill(jetPtbGENMatched,w);
	}
      }//end of RECO jet level
    }

    if (nbGEN>0){
      // all gen partons
      h_GENP_pt_b1_->Fill(p4GENmax_b.pt());
      h_GENP_eta_b1_->Fill(p4GENmax_b.eta());
      h_GENP_mass_b1_->Fill(p4GENmax_b.mass());

      // if parton is associated with a reconstructed jet
      if(hasbPartonJetMatched_){
	h_GENP_recoMatched_pt_b1_->Fill(p4GENmax_b.pt()); 
	h_GENP_recoMatched_eta_b1_->Fill(p4GENmax_b.eta());
	
	h2_GENpt_vs_RECOJetpt->Fill(jetPtbGENMatched,p4GENmax_b.pt());

	// if parton is associated with a reconstructed b-jet
	if(hasbPartonbJetMatched_){
	  h_GENP_recoBtagged_pt_b1_->Fill(p4GENmax_b.pt()); 
	  h_GENP_recoBtagged_eta_b1_->Fill(p4GENmax_b.eta());
	}
      }
    }
  }
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ZbbABCDMatrix methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ZbbABCDMatrix::book(){
  edm::Service<TFileService> fs;
  std::string DirName = &(*N_);
  TFileDirectory HistoVars = fs->mkdir(DirName);

  TH1F::SetDefaultSumw2(kTRUE);
  TH2F::SetDefaultSumw2(kTRUE);

  h_MllVsMET_ = HistoVars.make<TH2F>(N_+"_MllVsMET","l^{+}l^{-} Mass vs MET; #slash{E}_{T} (GeV); M(l^{+}l^{-}) (GeV)",100,0,200,150,0,300);
  h_MllVsNOMASSCUT_ = HistoVars.make<TH1F>(N_+"_MllNOMASSCUT","M(l^{+}l^{-}) - dilepton mass; M(l^{+}l^{-}) (GeV)",500,0,500);

  h_MmumuVsMET_ = HistoVars.make<TH2F>(N_+"_MmumuVsMET","#mu^{+}#mu^{-} Mass vs MET; #slash{E}_{T} (GeV); M(#mu^{+}#mu^{-}) (GeV)",100,0,200,150,0,300);
  h_MmumuVsNOMASSCUT_ = HistoVars.make<TH1F>(N_+"_MmumuNOMASSCUT","M(#mu^{+}#mu^{-}) - dilepton mass; M(#mu^{+}#mu^{-}) (GeV)",500,0,500);

  h_MeeVsMET_ = HistoVars.make<TH2F>(N_+"_MeeVsMET","e^{+}e^{-} Mass vs MET; #slash{E}_{T} (GeV); M(e^{+}e^{-}) (GeV)",100,0,200,150,0,300);
  h_MeeVsNOMASSCUT_ = HistoVars.make<TH1F>(N_+"_MeeNOMASSCUT","M(e^{+}e^{-}) - dilepton mass; M(e^{+}e^{-}) (GeV)",500,0,500);

}

void ZbbABCDMatrix::fill(const EventCategory& ec_mu,const EventCategory& ec_ele, edm::Handle<edm::View<pat::MET> > mets){

  //std::cout<<"ZbbABCDMatrix::fill"<<std::endl;

  Bool_t thecut_;   
  EventCategory ec; 
  
  if(ec_mu.isLL()){
    thecut_=ZbbUtils::ProduceTheCut(ec_mu,N_,"ZLL"); 
    ec = ec_mu;
  } else {
    thecut_=ZbbUtils::ProduceTheCut(ec_ele,N_,"ZLL"); 
    ec = ec_ele;
  }
  
  Double_t w = ZbbUtils::getTheWeight(ec,N_);

  if(thecut_){ 
    h_MllVsMET_->Fill((*mets)[0].et(),ec.bestZcandidate_.mass(),w);
    h_MllVsNOMASSCUT_->Fill(ec.bestZcandidate_.mass(),w);
  }
  return;
}

void ZbbABCDMatrix::fillMu(const EventCategory& ec_mu,edm::Handle<edm::View<pat::MET> > mets){

  //std::cout<<"ZbbABCDMatrix::fill"<<std::endl;

  Bool_t thecut_=ZbbUtils::ProduceTheCut(ec_mu,N_,"ZLL");   
  Double_t w = ZbbUtils::getTheWeight(ec_mu,N_);

  if(thecut_){ 
    h_MmumuVsMET_->Fill((*mets)[0].et(),ec_mu.bestZcandidate_.mass(),w);
    h_MmumuVsNOMASSCUT_->Fill(ec_mu.bestZcandidate_.mass(),w);
  }
  return;
}

void ZbbABCDMatrix::fillEle(const EventCategory& ec_ele,edm::Handle<edm::View<pat::MET> > mets){

  //std::cout<<"ZbbABCDMatrix::fill"<<std::endl;

  Bool_t thecut_=ZbbUtils::ProduceTheCut(ec_ele,N_,"ZLL");   
  Double_t w = ZbbUtils::getTheWeight(ec_ele,N_);

  if(thecut_){ 
    h_MeeVsMET_->Fill((*mets)[0].et(),ec_ele.bestZcandidate_.mass(),w);
    h_MeeVsNOMASSCUT_->Fill(ec_ele.bestZcandidate_.mass(),w);
  }

  return;
}
