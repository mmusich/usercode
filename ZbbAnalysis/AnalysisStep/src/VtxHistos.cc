#include "ZbbAnalysis/AnalysisStep/interface/EventCategory.h"
#include "ZbbAnalysis/AnalysisStep/interface/VtxHistos.h"
#include "ZbbAnalysis/AnalysisStep/interface/ZbbUtils.h"
#include "ZbbAnalysis/AnalysisStep/interface/VtxAssociatorsUtils.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include "TCut.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TProfile.h"
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

  h_ratio2b_antiBtag_ = HistoVars.make<TH1F>(N_+"_jetvtx_ratio2bAntiBtag","jet/vertex association ratio v2 using vertexing; jet-vertex association ratio",100,0,1);
  h_ratio2b_Btag_     = HistoVars.make<TH1F>(N_+"_jetvtx_ratio2bBtag","jet/vertex association ratio v2 using vertexing; jet-vertex association ratio",100,0,1);

  h_ratio2b_vsPU_ = HistoVars.make<TH2F>(N_+"_jetvtx_ratio2bVsPU","jet/vertex association ratio v2 using vertexing vs PU;N_{vtx} number of vertices; jet-vertex association ratio (#beta)",25,-0.5,24.5,100,0.,1.);
  p_ratio2b_vsPU_ = HistoVars.make<TProfile>(N_+"_jetvtx_ratio2bVsPU_prof","jet/vertex association ratio v2 using vertexing vs PU;N_{vtx} number of vertices; jet-vertex association ratio (#beta)",25,-0.5,24.5);

  // specific for 2b association
  h_PtJetFails2bAssoc_  = HistoVars.make<TH1F>(N_+"_PtJetFails2bAssoc","p_{T} of jet failing 2b jet/vertexc association; p_{T} of jet failing #beta cut",100.,0,100.);
  h_EtaJetFails2bAssoc_ = HistoVars.make<TH1F>(N_+"_EtaJetFails2bAssoc","#eta of jet failing 2b jet/vertexc association; #eta of jet failing #beta cut",100.,-3.,3.);

  h_ratio2bVsJetPt_  = HistoVars.make<TH2F>(N_+"_ratio2bVsJetPt","jet/vertex association 2b vs p_{T} of jet; p_{T} of jet (GeV);jet #beta",100.,0,100.,100,0,1);
  h_ratio2bVsJetEta_ = HistoVars.make<TH2F>(N_+"_ratio2bVsJetEta","jet/vertex association 2b vs #eta of jet; #eta of jet; jet #beta",100.,-3.,3.,100,0,1);
  
  h_Delta3DVtx0VtxSel_      = HistoVars.make<TH1F>(N_+"_Delta3DVtx0VtxSel","#Delta3D Vtx 0 w.r.t Vtx Selected; #Delta3D(PV_{0},PV_{sel}) (cm)",100,0,30);	 
  h_Delta3DVtxSelVtxEvtCat_ = HistoVars.make<TH1F>(N_+"_Delta3DVtxSelVtxEvtCat","#Delta3D Vtx Sel w.r.t Vtx Evt Category; #Delta3D(PV_{sel},PV_{evt cat}) (cm)",100,0,30);

  // histograms for PU dependence of track-jet association
  
  TFileDirectory BetaHistoVars = fs->mkdir(DirName+"/BetaAssociation");

  for(UInt_t nPU=0;nPU<25;nPU++){
    
    char histoname1[128], histoname2[128], histoname3[128];
    char histotitle1[128],histotitle2[128],histotitle3[128];
    
    sprintf(histoname1,"_jetvtx_ratio2b_%i_PUvtx",nPU);
    sprintf(histotitle1,"jet/vertex association ratio v2 using vertexing - %i PU; jet-vertex association ratio",nPU);
    h_ratio2b_inPUbins[nPU] = BetaHistoVars.make<TH1F>(N_+histoname1,histotitle1,100,0,1); 

    sprintf(histoname2,"_jetvtx_ratio2b_%i_PUvtxAntiBtag",nPU);
    sprintf(histotitle2,"anti btagged jet/vertex association ratio v2 using vertexing - %i PU; jet-vertex association ratio",nPU);
    h_ratio2b_inPUbins_antiBtag[nPU] = BetaHistoVars.make<TH1F>(N_+histoname2,histotitle2,100,0,1); 
    
    sprintf(histoname3,"_jetvtx_ratio2b_%i_PUvtxBtag",nPU);
    sprintf(histotitle3,"btagged jet/vertex association ratio v2 using vertexing - %i PU; jet-vertex association ratio",nPU);
    h_ratio2b_inPUbins_btag[nPU] = BetaHistoVars.make<TH1F>(N_+histoname3,histotitle3,100,0,1); 

    char histoname4[128], histoname5[128], histoname6[128],histoname7[128];
    char histotitle4[128],histotitle5[128],histotitle6[128],histotitle7[128];
    
    sprintf(histoname4,"_PtJetFails2bAssoc_%i_PUvtx",nPU);
    sprintf(histotitle4,"p_{T} of jet failing 2b jet/vertexc association - %i PU; p_{T} of jet failing #beta cut",nPU);
    h_PtJetFails2bAssoc_inPUbins[nPU] = BetaHistoVars.make<TH1F>(N_+histoname4,histotitle4,100,0.,100.); 
    
    sprintf(histoname5,"_EtaJetFails2bAssoc_%i_PUvtx",nPU);
    sprintf(histotitle5,"#eta of jet failing 2b jet/vertexc association - %i PU; #eta of jet failing #beta cut",nPU);
    h_EtaJetFails2bAssoc_inPUbins[nPU] = BetaHistoVars.make<TH1F>(N_+histoname5,histotitle5,100,-3,3); 
    
    sprintf(histoname6,"_ratio2bVsJetPt_%i_PUvtx",nPU);
    sprintf(histotitle6,"jet/vertex association 2b vs p_{T} of jet - %i PU;p_{T} of jet (GeV); jet #beta",nPU);
    h_ratio2bVsJetPt_inPUbins[nPU] = BetaHistoVars.make<TH2F>(N_+histoname6,histotitle6,100.,0,100.,100,0,1); 
    
    sprintf(histoname7,"_ratio2bVsJetEta_%i_PUvtx",nPU);
    sprintf(histotitle7,"jet/vertex association 2b vs #eta of jet - %i PU;#eta of jet; jet #beta",nPU);
    h_ratio2bVsJetEta_inPUbins[nPU] = BetaHistoVars.make<TH2F>(N_+histoname7,histotitle7,100.,-3.,3.,100,0,1); 

  }

  h_goodevent_ = HistoVars.make<TH1F>(N_+"_jetvtx_goodevent","pass or not Z+jet to vertex association",2,0,2);

  TString goodeventBinLabels[2] ={"Z+Jet associated","Z+Jet not associated"};
   
  for(UInt_t bin=1; bin<=2; bin++){
    h_goodevent_->GetXaxis()->SetBinLabel(bin,goodeventBinLabels[bin-1]);    
  }

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
    
    h_nvertices_->Fill(nvvertex,w);
    // h_nvertices_->Fill(vertices.size());
  
    // cut on the signficance
    Double_t sigcut = 5;
    
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
      
      Double_t ptsum = 0.;
      Double_t ptsumx = 0.;
      Double_t ptsumy = 0.;
      Double_t ptsumall = 0.;
      Double_t ratio = 0;

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
	  if( fabs(sig)<sigcut ){
	    ptsum  += jetPFConst->pt();
	    ptsumx += jetPFConst->px();
	    ptsumy += jetPFConst->py();
	  }
	  ptsumall += jetPFConst->pt();
	}
	
    	h_ratio1_->Fill(ptsum/(*jets)[i].et(),w); 
    
    	if(ptsumall>0) {
    	  ratio = ptsum/ptsumall;
    	} else{
    	  ratio = -1;
    	}
    
    	h_ratio2_->Fill(ratio,w);
    	h_ratio3_->Fill( pow((pow(ptsumx,2)+pow(ptsumy,2)),0.5) /(*jets)[i].et(),w);
	
	// method using vertexing
	
	ptsum = 0.;
	ptsumx = 0.;
	ptsumy = 0.;
	ptsumall = 0.;
	ratio = 0;
	
	for(UInt_t j=0; j< jet_constsize; j++){
	  
	  reco::PFCandidatePtr jetPFConst = (*jets)[i].getPFConstituent(j);

	  if(jetPFConst->trackRef().isNull()){
	    continue;
	  }
	  
	  for( reco::Vertex::trackRef_iterator tk = vertex.tracks_begin(); tk<vertex.tracks_end(); tk++ ){
	    if(tk->key() == jetPFConst->trackRef().key()) {
	      ptsum    += jetPFConst->pt();
	      ptsumx   += jetPFConst->px();
	      ptsumy   += jetPFConst->py();
	    } // closes if tkref = pfcandidate 
	  } // closes loop on vertex tracks
	  ptsumall += jetPFConst->pt();
	} // closes loop on pfconstituents

	h_ratio1b_->Fill(ptsum/(*jets)[i].et(),w); 
	
    	if(ptsumall>0) {
    	  ratio = ptsum/ptsumall;
    	} else{
    	  ratio = -1;
    	}
	
    	h_ratio2b_->Fill(ratio,w);

	h_ratio2bVsJetEta_->Fill((*jets)[i].eta(),ratio,w); 
	h_ratio2bVsJetPt_->Fill((*jets)[i].pt(),ratio,w);
  
	if(ratio<0.15){
	  h_PtJetFails2bAssoc_->Fill((*jets)[i].pt(),w);  
	  h_EtaJetFails2bAssoc_->Fill((*jets)[i].eta(),w); 
	}

	for(UInt_t nPU=0;nPU<25;nPU++){
	  if((nvvertex-1)==nPU){
	    // fill for all jets
	    h_ratio2b_inPUbins[nPU]->Fill(ratio,w);
	    if(!ZbbUtils::isBJet((*jets)[i],bTagAlgoWP) ) {         // fill for anti b-tagged jets
	      h_ratio2b_inPUbins_antiBtag[nPU]->Fill(ratio,w); 
	      h_ratio2b_antiBtag_->Fill(ratio,w); 
	    } else {                                                // fill for b-tagged jets
	      h_ratio2b_inPUbins_btag[nPU]->Fill(ratio,w);
	      h_ratio2b_Btag_->Fill(ratio,w);    
	    }

	    if(ratio<0.15){
	      h_PtJetFails2bAssoc_inPUbins[nPU]->Fill((*jets)[i].pt(),w);
	      h_EtaJetFails2bAssoc_inPUbins[nPU]->Fill((*jets)[i].eta(),w);
	    }
	    
	    h_ratio2bVsJetPt_inPUbins[nPU]->Fill((*jets)[i].pt(),ratio,w);
	    h_ratio2bVsJetEta_inPUbins[nPU]->Fill((*jets)[i].eta(),ratio,w);
	    
	  }
	}

	h_ratio2b_vsPU_->Fill(nvvertex-1,ratio,w);
	p_ratio2b_vsPU_->Fill(nvvertex-1,ratio,w);
	  
    	h_ratio3b_->Fill( pow((pow(ptsumx,2)+pow(ptsumy,2)),0.5) /(*jets)[i].et(),w);

      } else {
    	continue;
      }
    }
    h_goodevent_->Fill(VtxAssociatorsUtils::checkVertexAssociation(ec.bestZcandidate_,jets,vertices),w); 
  }
}
