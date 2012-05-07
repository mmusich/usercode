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
  TString N_;
    
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
  TH1F* h_ratio1_;     
  TH1F* h_ratio2_;     
  TH1F* h_ratio3_;   
  TH1F* h_ratio1b_;   
  TH1F* h_ratio2b_;   
  TH1F* h_ratio3b_;   

  TH1F* h_ratio2b_antiBtag_;   
  TH1F* h_ratio2b_Btag_;   
  
  TH1F* h_Delta3DVtx0VtxSel_;     
  TH1F* h_Delta3DVtxSelVtxEvtCat_;

  TH1F* h_ratio2b_inPUbins[25];
  TH1F* h_ratio2b_inPUbins_antiBtag[25];
  TH1F* h_ratio2b_inPUbins_btag[25];

  TH1F* h_goodevent_;  
 
  TH1F* h_PtJetFails2bAssoc_;
  TH1F* h_EtaJetFails2bAssoc_;
  
  TH2F *h_ratio2bVsJetEta_;
  TH2F *h_ratio2bVsJetPt_;
 
  TH1F *h_PtJetFails2bAssoc_inPUbins[25];
  TH1F *h_EtaJetFails2bAssoc_inPUbins[25];
  
  TH2F *h_ratio2bVsJetEta_inPUbins[25];
  TH2F *h_ratio2bVsJetPt_inPUbins[25];

  TH2F     *h_ratio2b_vsPU_;
  TProfile *p_ratio2b_vsPU_;
  
  public:
  
  VtxHistos(const TString& name) : 
    N_(name){}

  void book();
  void fill(const EventCategory& ec_mu, const EventCategory& ec_ele, edm::View<reco::Vertex> vertices, edm::Handle<edm::View<pat::Jet> > jets,std::string bTagAlgoWP);
  ~VtxHistos(){} 

};

#endif
