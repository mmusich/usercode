#ifndef ZbbAnalysis_AnalysisStep_VtxHistos_h
#define ZbbAnalysis_AnalysisStep_VtxHistos_h

#include "ZbbAnalysis/AnalysisStep/interface/ZbbUtils.h"
#include "ZbbAnalysis/AnalysisStep/interface/EventCategory.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

/** \class VtxHistos
 *
 * A set of histograms for vertex association.
 *
 *  $Date: 2012/02/14 10:49:39 $
 *  $Revision: 1.6 $
 *  \author M. Musich - Torino
 */

#include <vector>
#include <iostream>
#include <TString.h>
#include <string>
class TH1F;
class TH2F;
class TH2I;
class TProfile;

using namespace ZbbUtils;

class VtxHistos{
  
  AcceptanceCuts lCuts_;
  TString N_;
    
  // lepton and vertex miscellanea
  TH1F* h_nvertices_;  
  TH1F* h_vx_;         
  TH1F* h_vy_;         
  TH1F* h_vz_;         
  TH1F* h_vxerr_;      
  TH1F* h_vyerr_;      
  TH1F* h_vzerr_;      
  TH1F* h_lepton_dz_;
  TH1F* h_l1v_dz_;     
  TH1F* h_l2v_dz_;     
  TH1F* h_distance_;   
  TH1F* h_sig_;   
  
  // jet vertex association figures of merit
  TH1F* h_ratio1_;     
  TH1F* h_ratio2_;     
  TH1F* h_ratio3_;   
  TH1F* h_ratio1b_;   
  TH1F* h_ratio2b_;   
  TH1F* h_ratio3b_;  
  TH1F* h_beta_;
  TH1F* h_betastar_;
  
  // ============== beta association =============
  // beta in bins of PU and asking or not b-tag
  TH1F* h_beta_antiBtag_; 
  TH1F* h_beta_Btag_;     
  
  TH1F* h_beta_inPUbins[25];
  TH1F* h_beta_inPUbins_antiBtag[25];
  TH1F* h_beta_inPUbins_btag[25];
 
  // if the jet fails beta association
  TH1F* h_PtJetFailsbetaAssoc_;
  TH1F* h_EtaJetFailsbetaAssoc_;
  TH1F* h_PtJetFailsbetaAssoc_inPUbins[25];
  TH1F* h_EtaJetFailsbetaAssoc_inPUbins[25];
  
  // beta correlation with jet kinematics
  TH2F *h_betaVsJetEta_;
  TH2F *h_betaVsJetPt_;
  TH2F *h_betaVsJetEta_inPUbins[25];
  TH2F *h_betaVsJetPt_inPUbins[25];
  
  // beta vs PU
  TH2F     *h_beta_vsPU_;
  TProfile *p_beta_vsPU_;

  // ============== betastar association ==========
  // betastar in bins of PU and asking or not b-tag
  TH1F* h_betastar_antiBtag_; 
  TH1F* h_betastar_Btag_;
 
  TH1F* h_betastar_inPUbins[25];
  TH1F* h_betastar_inPUbins_antiBtag[25];
  TH1F* h_betastar_inPUbins_btag[25];

  // if the jet fails betastar association
  TH1F* h_PtJetFailsbetastarAssoc_;
  TH1F* h_EtaJetFailsbetastarAssoc_;
  TH1F* h_PtJetFailsbetastarAssoc_inPUbins[25];
  TH1F* h_EtaJetFailsbetastarAssoc_inPUbins[25];
  
  // betastar correlation with jet kinematics
  TH2F *h_betastarVsJetEta_;
  TH2F *h_betastarVsJetPt_;
  TH2F *h_betastarVsJetEta_inPUbins[25];
  TH2F *h_betastarVsJetPt_inPUbins[25];

  // beta star vs PU
  TH2F     *h_betastar_vsPU_;
  TProfile *p_betastar_vsPU_;

  // just my check
  TH1F* h_Delta3DVtx0VtxSel_;     
  TH1F* h_Delta3DVtxSelVtxEvtCat_;

  TH1F* h_goodevent_;  

  public:
  
  VtxHistos(const TString& name, AcceptanceCuts theCuts) : 
    lCuts_(theCuts),
    N_(name){}

  void book();
  void fill(const EventCategory& ec_mu, const EventCategory& ec_ele, edm::View<reco::Vertex> vertices, edm::Handle<edm::View<pat::Jet> > jets,std::string bTagAlgoWP);
  ~VtxHistos(){} 

};

#endif
