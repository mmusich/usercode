#include "ZbbAnalysis/AnalysisStep/interface/EventCategory.h"
#include "ZbbAnalysis/AnalysisStep/interface/VtxHistos.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbUtils.h"
#include "ZbbAnalysis/AnalysisStep/interface/VtxAssociatorsUtils.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include "TProfile.h"
#include "TCut.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TString.h"
#include "TLorentzVector.h"

#include <iostream>
#include <algorithm>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Vxt associator methods
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void VtxHistos::book(){
  edm::Service<TFileService> fs;
  std::string DirName = &(*N_);
  TFileDirectory HistoVars = fs->mkdir(DirName);

  TH1F::SetDefaultSumw2(kTRUE);

  h_nvertices_ = HistoVars.make<TH1F>(N_+"_nvertices","Nb. of valid vertices; N_{VTX}",30,-0.5,29.5);			    
  h_vx_        = HistoVars.make<TH1F>(N_+"_vx","vx; x position of PV (cm)",100,-0.5,0.5);					    
  h_vy_        = HistoVars.make<TH1F>(N_+"_vy","vy; y position of PV (cm)",100,-0.5,0.5);					    
  h_vz_        = HistoVars.make<TH1F>(N_+"_vz","vz; z position of PV (cm)",100,-25,25);		
  h_vxerr_     = HistoVars.make<TH1F>(N_+"_vxerr","vxerr; error on x position of PV (cm)",100,0,0.01);					    
  h_vyerr_     = HistoVars.make<TH1F>(N_+"_vyerr","vyerr; error on y position of PV (cm)",100,0,0.01);					    
  h_vzerr_     = HistoVars.make<TH1F>(N_+"_vzerr","vzerr; error on z position of PV (cm)",100,0,0.02);					    
  h_lepton_dz_ = HistoVars.make<TH1F>(N_+"_lepton_dz","z distance between the two Z leptons; |#Deltaz(ll)| (cm)",100,0,0.1);
  h_l1v_dz_    = HistoVars.make<TH1F>(N_+"_l1v_dz","z distance between lepton and vertex; |#Deltaz(l_{1},Z)| (cm)",100,0,0.1);	    
  h_l2v_dz_    = HistoVars.make<TH1F>(N_+"_l2v_dz","z distance between lepton and vertex; |#Deltaz(l_{2},Z)| (cm)",100,0,0.1);	    
  h_distance_  = HistoVars.make<TH1F>(N_+"_tkvtx_distance","vertex/track distance in z; |#Deltaz(tk,PV)| (cm)",100,-10,10);	    
  h_sig_       = HistoVars.make<TH1F>(N_+"_tkvtx_sig","vertex/track significance in z;|#Deltaz(tk,PV)|/#sigma_{#Deltaz}",100,0,10);		    
  h_ratio1_    = HistoVars.make<TH1F>(N_+"_jetvtx_ratio1","jet/vertex association ratio v1; jet-vertex association ratio",100,0,1);	    
  h_ratio2_    = HistoVars.make<TH1F>(N_+"_jetvtx_ratio2","jet/vertex association ratio v2; jet-vertex association ratio",100,0,1);	    
  h_ratio3_    = HistoVars.make<TH1F>(N_+"_jetvtx_ratio3","jet/vertex association ratio v3; jet-vertex association ratio",100,0,1);	    
  h_ratio1b_   = HistoVars.make<TH1F>(N_+"_jetvtx_ratio1b","jet/vertex association ratio v1 using vertexing; jet-vertex association ratio",100,0,1);
  h_ratio2b_   = HistoVars.make<TH1F>(N_+"_jetvtx_ratio2b","jet/vertex association ratio v2 using vertexing; jet-vertex association ratio",100,0,1);
  h_ratio3b_   = HistoVars.make<TH1F>(N_+"_jetvtx_ratio3b","jet/vertex association ratio v3 using vertexing; jet-vertex association ratio",100,0,1);
  
  h_beta_      = HistoVars.make<TH1F>(N_+"_beta","jet #beta; jet-vertex association #beta (%)",200,0.,1.);
  h_betastar_  = HistoVars.make<TH1F>(N_+"_betastar","jet #beta; jet-vertex association #beta* (%)",200,0.,1.);

  // specific for beta association
  h_beta_antiBtag_ = HistoVars.make<TH1F>(N_+"_betaAntiBtag","jet #beta; jet-vertex association #beta (%)",200,0.,1.);
  h_beta_Btag_     = HistoVars.make<TH1F>(N_+"_betaBtag","jet #beta; jet-vertex association #beta (%)",200,0.,1.);

  h_beta_vsPU_ = HistoVars.make<TH2F>(N_+"_betaVsPU","jet #beta vs PU;N_{vtx} number of vertices; jet-vertex association #beta (%)",25,-0.5,24.5,200,0.,1.);
  p_beta_vsPU_ = HistoVars.make<TProfile>(N_+"_betaVsPU_prof","jet #beta vs PU;N_{vtx} number of vertices; jet-vertex association #beta (%)",25,-0.5,24.5);

  h_PtJetFailsbetaAssoc_  = HistoVars.make<TH1F>(N_+"_PtJetFailsbetaAssoc","p_{T} of jet failing #beta association cut; p_{T} of jet failing #beta cut (GeV)",100.,0,100.);
  h_EtaJetFailsbetaAssoc_ = HistoVars.make<TH1F>(N_+"_EtaJetFailsbetaAssoc","#eta of jet failing #beta association cut; #eta of jet failing #beta cut",100.,-3.,3.);
  
  h_betaVsJetPt_  = HistoVars.make<TH2F>(N_+"_betaVsJetPt","jet #beta vs p_{T} of jet; p_{T} of jet (GeV);jet #beta (%)",100.,0,100.,200,0.,1.);
  h_betaVsJetEta_ = HistoVars.make<TH2F>(N_+"_betaVsJetEta","jet #beta vs #eta of jet; #eta of jet; jet #beta (%)",100.,-3.,3.,200,0.,1.);

  // specific for betastar association
  h_betastar_antiBtag_ = HistoVars.make<TH1F>(N_+"_betastarAntiBtag","jet #beta*; jet-vertex association #beta* (%)",200,0.,1.);
  h_betastar_Btag_     = HistoVars.make<TH1F>(N_+"_betastarBtag","jet #beta*; jet-vertex association #beta* (%)",200,0.,1.);

  h_betastar_vsPU_ = HistoVars.make<TH2F>(N_+"_betastarVsPU","jet #beta* vs PU;N_{vtx} number of vertices; jet-vertex association #beta* (%)",25,-0.5,24.5,200,0.,1.);
  p_betastar_vsPU_ = HistoVars.make<TProfile>(N_+"_betastarVsPU_prof","jet #beta* vs PU;N_{vtx} number of vertices; jet-vertex association #beta* (%)",25,-0.5,24.5);

  h_PtJetFailsbetastarAssoc_  = HistoVars.make<TH1F>(N_+"_PtJetFailsbetastarAssoc","p_{T} of jet failing #beta* association cut; p_{T} of jet failing #beta* cut (GeV)",100.,0,100.);
  h_EtaJetFailsbetastarAssoc_ = HistoVars.make<TH1F>(N_+"_EtaJetFailsbetastarAssoc","#eta of jet failing #beta* association cut; #eta of jet failing #beta* cut",100.,-3.,3.);

  h_betastarVsJetPt_  = HistoVars.make<TH2F>(N_+"_betastarVsJetPt","jet #beta* vs p_{T} of jet; p_{T} of jet (GeV);jet #beta* (%)",100.,0,100.,200,0.,1.);
  h_betastarVsJetEta_ = HistoVars.make<TH2F>(N_+"_betastarVsJetEta","jet #beta* vs #eta of jet; #eta of jet; jet #beta* (%)",100.,-3.,3.,200,0.,1.);

  // histograms for PU dependence of track-jet association
  TFileDirectory BetaHistoVars = fs->mkdir(DirName+"/BetaAssociation");
  TFileDirectory BetaStarHistoVars = fs->mkdir(DirName+"/BetaStarAssociation");

  for(UInt_t nPU=0;nPU<25;nPU++){
      
    // beta association part
    
    char hname1beta[128], hname2beta[128], hname3beta[128];
    char htitle1beta[128],htitle2beta[128],htitle3beta[128];
    
    sprintf(hname1beta,"_jetvtx_beta_%i_PUvtx",nPU);
    sprintf(htitle1beta,"jet #beta - %i PU; jet-vertex association #beta",nPU);
    h_beta_inPUbins[nPU] = BetaHistoVars.make<TH1F>(N_+hname1beta,htitle1beta,200,0.,1.); 

    sprintf(hname2beta,"_jetvtx_beta_%i_PUvtxAntiBtag",nPU);
    sprintf(htitle2beta,"anti btagged jet #beta - %i PU; jet-vertex association #beta",nPU);
    h_beta_inPUbins_antiBtag[nPU] = BetaHistoVars.make<TH1F>(N_+hname2beta,htitle2beta,200,0.,1.); 
    
    sprintf(hname3beta,"_jetvtx_beta_%i_PUvtxBtag",nPU);
    sprintf(htitle3beta,"btagged jet #beta - %i PU; jet-vertex association #beta",nPU);
    h_beta_inPUbins_btag[nPU] = BetaHistoVars.make<TH1F>(N_+hname3beta,htitle3beta,200,0.,1.); 

    char hname4[128], hname5[128], hname6[128],hname7[128];
    char htitle4[128],htitle5[128],htitle6[128],htitle7[128];
    
    sprintf(hname4,"_PtJetFailsbetaAssoc_%i_PUvtx",nPU);
    sprintf(htitle4,"p_{T} of jet failing #beta association cut- %i PU; p_{T} of jet failing #beta cut",nPU);
    h_PtJetFailsbetaAssoc_inPUbins[nPU] = BetaHistoVars.make<TH1F>(N_+hname4,htitle4,100,0.,100.); 
    
    sprintf(hname5,"_EtaJetFailsbetaAssoc_%i_PUvtx",nPU);
    sprintf(htitle5,"#eta of jet failing #beta association - %i PU; #eta of jet failing #beta cut",nPU);
    h_EtaJetFailsbetaAssoc_inPUbins[nPU] = BetaHistoVars.make<TH1F>(N_+hname5,htitle5,100,-3,3); 
    
    sprintf(hname6,"_betaVsJetPt_%i_PUvtx",nPU);
    sprintf(htitle6,"jet #beta vs p_{T} of jet - %i PU;p_{T} of jet (GeV); jet #beta",nPU);
    h_betaVsJetPt_inPUbins[nPU] = BetaHistoVars.make<TH2F>(N_+hname6,htitle6,100.,0,100.,200,0.,1.); 
    
    sprintf(hname7,"_betaVsJetEta_%i_PUvtx",nPU);
    sprintf(htitle7,"jet #beta vs #eta of jet - %i PU;#eta of jet; jet #beta",nPU);
    h_betaVsJetEta_inPUbins[nPU] = BetaHistoVars.make<TH2F>(N_+hname7,htitle7,100.,-3.,3.,200,0.,1.); 

    // betastar association part

    char hname1betastar[128], hname2betastar[128], hname3betastar[128];
    char htitle1betastar[128],htitle2betastar[128],htitle3betastar[128];
    
    sprintf(hname1betastar,"_jetvtx_betastar_%i_PUvtx",nPU);
    sprintf(htitle1betastar,"jet #beta* - %i PU; jet-vertex association #beta*",nPU);
    h_betastar_inPUbins[nPU] = BetaStarHistoVars.make<TH1F>(N_+hname1betastar,htitle1betastar,200,0.,1.); 

    sprintf(hname2betastar,"_jetvtx_betastar_%i_PUvtxAntiBtag",nPU);
    sprintf(htitle2betastar,"anti btagged jet #beta* - %i PU; jet-vertex association #beta*",nPU);
    h_betastar_inPUbins_antiBtag[nPU] = BetaStarHistoVars.make<TH1F>(N_+hname2betastar,htitle2betastar,200,0.,1.); 
    
    sprintf(hname3betastar,"_jetvtx_betastar_%i_PUvtxBtag",nPU);
    sprintf(htitle3betastar,"btagged jet #beta* - %i PU; jet-vertex association #beta*",nPU);
    h_betastar_inPUbins_btag[nPU] = BetaStarHistoVars.make<TH1F>(N_+hname3betastar,htitle3betastar,200,0.,1.);

    char hname8[128], hname9[128], hname10[128],hname11[128];
    char htitle8[128],htitle9[128],htitle10[128],htitle11[128];
    
    sprintf(hname8,"_PtJetFailsbetastarAssoc_%i_PUvtx",nPU);
    sprintf(htitle8,"p_{T} of jet failing #beta* association cut- %i PU; p_{T} of jet failing #beta* cut",nPU);
    h_PtJetFailsbetastarAssoc_inPUbins[nPU] = BetaStarHistoVars.make<TH1F>(N_+hname8,htitle8,100,0.,100.); 
    
    sprintf(hname9,"_EtaJetFailsbetastarAssoc_%i_PUvtx",nPU);
    sprintf(htitle9,"#eta of jet failing #beta* association - %i PU; #eta of jet failing #beta* cut",nPU);
    h_EtaJetFailsbetastarAssoc_inPUbins[nPU] = BetaStarHistoVars.make<TH1F>(N_+hname9,htitle9,100,-3,3); 
    
    sprintf(hname10,"_betastarVsJetPt_%i_PUvtx",nPU);
    sprintf(htitle10,"jet #beta* vs p_{T} of jet - %i PU;p_{T} of jet (GeV); jet #beta*",nPU);
    h_betastarVsJetPt_inPUbins[nPU] = BetaStarHistoVars.make<TH2F>(N_+hname10,htitle10,100.,0,100.,200,0.,1.); 
    
    sprintf(hname11,"_betastarVsJetEta_%i_PUvtx",nPU);
    sprintf(htitle11,"jet #beta* vs #eta of jet - %i PU;#eta of jet; jet #beta*",nPU);
    h_betastarVsJetEta_inPUbins[nPU] = BetaStarHistoVars.make<TH2F>(N_+hname11,htitle11,100.,-3.,3.,200,0.,1.); 
 
  }

  // decision if jet is good for beta cut
  h_goodevent_ = HistoVars.make<TH1F>(N_+"_jetvtx_goodevent","pass or not Z+jet to vertex association",2,0,2);
  TString goodeventBinLabels[2] ={"Z+Jet associated","Z+Jet not associated"};
  for(UInt_t bin=1; bin<=2; bin++){
    h_goodevent_->GetXaxis()->SetBinLabel(bin,goodeventBinLabels[bin-1]);    
  }

  // just my checks
  h_Delta3DVtx0VtxSel_      = HistoVars.make<TH1F>(N_+"_Delta3DVtx0VtxSel","#Delta3D Vtx 0 w.r.t Vtx Selected; #Delta3D(PV_{0},PV_{sel}) (cm)",100,0,30);	 
  h_Delta3DVtxSelVtxEvtCat_ = HistoVars.make<TH1F>(N_+"_Delta3DVtxSelVtxEvtCat","#Delta3D Vtx Sel w.r.t Vtx Evt Category; #Delta3D(PV_{sel},PV_{evt cat}) (cm)",100,0,30);

}

void VtxHistos::fill(const EventCategory& ec_mu, const EventCategory& ec_ele, edm::View<reco::Vertex> vertices, edm::Handle<edm::View<pat::Jet> > jets,std::string bTagAlgoWP){
    
  Bool_t thecut_;   
  EventCategory ec; 
  
  if(ec_mu.isZLL()){
    thecut_= ZbbUtils::ProduceTheCut(ec_mu,N_,"");
    ec = ec_mu;
  } else {
    thecut_ = ZbbUtils::ProduceTheCut(ec_ele,N_,"");
    ec = ec_ele;
  }

  Double_t w = ZbbUtils::getTheWeight(ec,N_);

  if(thecut_){
    
    unsigned int vertexCollectionSize = vertices.size();  
    unsigned int nvvertex = 0;
    for (unsigned int i=0; i<vertexCollectionSize; i++) {
      const reco::Vertex& vertex = vertices.at(i);
      if (vertex.isValid()) nvvertex++;
    }
    
    // lists of track refs necessary for PU/PV beta/betastar
    std::list<int> trackrefs_PV = VtxAssociatorsUtils::buildTrackRefs(vertices,false);
    std::list<int> trackrefs_PU = VtxAssociatorsUtils::buildTrackRefs(vertices,true);
    
    h_nvertices_->Fill(nvvertex,w);
    // h_nvertices_->Fill(vertices.size());
      
    // select the vertex associated to the Z  
    reco::Vertex vertex = VtxAssociatorsUtils::findPrimaryVertex(ec.bestZcandidate_,vertices);
    reco::Vertex vertexFromEvtCat = ec.theZvertex_;

    const math::XYZPoint zvtx(vertex.position().x(),vertex.position().y(),vertex.position().z());
    const math::XYZPoint zvtxEvCat(vertexFromEvtCat.position().x(),vertexFromEvtCat.position().y(),vertexFromEvtCat.position().z());
    const math::XYZPoint pvtx(vertices.at(0).position().x(),vertices.at(0).position().y(),vertices.at(0).position().z());
    
    h_Delta3DVtx0VtxSel_->Fill(TMath::Sqrt( pow(zvtx.x()-pvtx.x(),2) + pow(zvtx.y()-pvtx.y(),2) + pow(zvtx.z()-pvtx.z(),2) ));     
    h_Delta3DVtxSelVtxEvtCat_->Fill(TMath::Sqrt( pow(zvtx.x()-zvtxEvCat.x(),2) + pow(zvtx.y()-zvtxEvCat.y(),2) + pow(zvtx.z()-zvtxEvCat.z(),2) ) );
    
    h_vx_->Fill(vertex.x(),w); 
    h_vy_->Fill(vertex.y(),w); 
    h_vz_->Fill(vertex.z(),w);
    h_vxerr_->Fill(vertex.xError(),w);
    h_vyerr_->Fill(vertex.yError(),w);
    h_vzerr_->Fill(vertex.zError(),w);
    
    // relevant quantities to monitor: Z vs primary vertex
    const reco::Candidate* lepton1 = ec.bestZcandidate_.daughter(0);
    const reco::Candidate* lepton2 = ec.bestZcandidate_.daughter(1);
    
    h_lepton_dz_->Fill(abs(lepton1->vz()-lepton2->vz()),w);  
    h_l1v_dz_->Fill(abs(lepton1->vz()-vertex.z()),w);       
    h_l2v_dz_->Fill(abs(lepton2->vz()-vertex.z()),w);       
    
    // relevant quantities to monitor: jets vs primary vertex
    for(unsigned i=0; i<jets->size(); ++i){
     
      if((*jets)[i].isPFJet()){

	UInt_t jet_constsize = (*jets)[i].getPFConstituents().size();

	for(UInt_t j=0; j< jet_constsize; j++){

	  reco::PFCandidatePtr jetPFConst = (*jets)[i].getPFConstituent(j);

	  //make sure the object is usable
	  //the last condition is a fix if we miss muons and electrons in the file, for rare occurences... 
	  //apparently something in the vz() calculation.
	  
          if (!(jetPFConst.isAvailable()) || (jetPFConst.isNull())){
            continue;
	  }
          if (jetPFConst->trackRef().isNull()){
	    continue;
	  }
          if (jetPFConst->muonRef().isNonnull()  || jetPFConst->gsfTrackRef().isNonnull()  ){
	    continue;
	  }
         
	  Double_t distance = (jetPFConst->vz() - vertex.z());
	  h_distance_->Fill(distance,w);
	  
	  Double_t error = pow( (pow((jetPFConst->trackRef()->dzError()),2) + pow(vertex.zError(),2)),0.5);			
	  Double_t sig = distance/error;
	  
	  h_sig_->Fill(sig,w); 
	
	}
	
	// calculate the map of jet-vertex associations
	std::map<TString,Double_t> theJVFmap = VtxAssociatorsUtils::calculateJetVertexAssociation((*jets)[i],vertex);

	// method using jet-vtx distance association
    	h_ratio1_->Fill(theJVFmap["1"],w);   
    	h_ratio2_->Fill(theJVFmap["2"],w);
    	h_ratio3_->Fill(theJVFmap["3"],w);
	
	// method using vertexing	
	h_ratio1b_->Fill(theJVFmap["1b"],w); 
    	h_ratio2b_->Fill(theJVFmap["2b"],w);
	h_ratio3b_->Fill(theJVFmap["3b"],w);
	
	float mybeta = VtxAssociatorsUtils::beta((*jets)[i],trackrefs_PV);
	float mybetastar = VtxAssociatorsUtils::betaStar((*jets)[i],trackrefs_PU);
	
	h_beta_->Fill(mybeta,w);
	h_betastar_->Fill(mybetastar,w); 

	h_betaVsJetEta_->Fill((*jets)[i].eta(),mybeta,w); 
	h_betaVsJetPt_->Fill((*jets)[i].pt(),mybeta,w);

	h_betastarVsJetEta_->Fill((*jets)[i].eta(),mybetastar,w); 
	h_betastarVsJetPt_->Fill((*jets)[i].pt(),mybetastar,w);

	// faling beta association
	if(mybeta<lCuts_.betaCut_){
	  h_PtJetFailsbetaAssoc_->Fill((*jets)[i].pt(),w);  
	  h_EtaJetFailsbetaAssoc_->Fill((*jets)[i].eta(),w); 
	}

	// failing beta star association
	if(mybetastar<lCuts_.betaStarCut_){
	  h_PtJetFailsbetastarAssoc_->Fill((*jets)[i].pt(),w);  
	  h_EtaJetFailsbetastarAssoc_->Fill((*jets)[i].eta(),w); 
	}

	for(UInt_t nPU=0;nPU<25;nPU++){
	  if((nvvertex-1)==nPU){
	    
	    // fill for all jets
	    h_beta_inPUbins[nPU]->Fill(mybeta,w);
	    h_betastar_inPUbins[nPU]->Fill(mybetastar,w);

	    if(!ZbbUtils::isBJet((*jets)[i],bTagAlgoWP) ) {         // fill for anti b-tagged jets
	      h_beta_inPUbins_antiBtag[nPU]->Fill(mybeta,w);  
	      h_betastar_inPUbins_antiBtag[nPU]->Fill(mybetastar,w);  
	      h_beta_antiBtag_->Fill(mybeta,w);
	      h_betastar_antiBtag_->Fill(mybetastar,w);
	    } else {                                                // fill for b-tagged jets
	      h_beta_inPUbins_btag[nPU]->Fill(mybeta,w);
	      h_betastar_inPUbins_btag[nPU]->Fill(mybetastar,w);
	      h_beta_Btag_->Fill(mybeta,w);    
	      h_betastar_Btag_->Fill(mybetastar,w); 
	    }

	    // failing beta association
	    if(mybeta<lCuts_.betaCut_){
	      h_PtJetFailsbetaAssoc_inPUbins[nPU]->Fill((*jets)[i].pt(),w);
	      h_EtaJetFailsbetaAssoc_inPUbins[nPU]->Fill((*jets)[i].eta(),w);
	    }

	    // failing betastar association
	    if(mybetastar<lCuts_.betaStarCut_){
	      h_PtJetFailsbetastarAssoc_inPUbins[nPU]->Fill((*jets)[i].pt(),w);
	      h_EtaJetFailsbetastarAssoc_inPUbins[nPU]->Fill((*jets)[i].eta(),w);
	    }
	    
	    h_betaVsJetPt_inPUbins[nPU]->Fill((*jets)[i].pt(),mybeta,w);
	    h_betaVsJetEta_inPUbins[nPU]->Fill((*jets)[i].eta(),mybeta,w);
	    
	    h_betastarVsJetPt_inPUbins[nPU]->Fill((*jets)[i].pt(),mybetastar,w);
	    h_betastarVsJetEta_inPUbins[nPU]->Fill((*jets)[i].eta(),mybetastar,w);
	    
	  }
	}
	
	h_beta_vsPU_->Fill(nvvertex-1,mybeta,w);
	p_beta_vsPU_->Fill(nvvertex-1,mybeta,w);

	h_betastar_vsPU_->Fill(nvvertex-1,mybetastar,w);
	p_betastar_vsPU_->Fill(nvvertex-1,mybetastar,w);

      } else {
    	continue;
      }
    }
    h_goodevent_->Fill(VtxAssociatorsUtils::checkVertexAssociation(ec.bestZcandidate_,jets,vertices),w); 
  }
}
