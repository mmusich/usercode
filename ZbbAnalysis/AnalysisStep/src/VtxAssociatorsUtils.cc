#include "ZbbAnalysis/AnalysisStep/interface/VtxAssociatorsUtils.h"

using namespace VtxAssociatorsUtils;

///////////////////////////////////////////////////////////////
Bool_t 
VtxAssociatorsUtils::checkVertexAssociation(const reco::CompositeCandidate& ZCand, edm::Handle<edm::View<pat::Jet> > jets, edm::View<reco::Vertex> vertices){
  // condition 1: both leptons are from the same vertex
  if(!loosezVertex(ZCand, 0.1)){
    return false;
  }
  //condition 2: strict lepton to vertex association
  const reco::Vertex& vertex = findPrimaryVertex(ZCand,vertices);
  if(!tightzVertex(ZCand, 0.1,vertex)){
    return false;
  }
  //condition 3: at least one of the jets is associated to the z vertex.
  for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
    if(jetVertex(vertex,*jet,2,2,0.5)){
      return true;
    }
  }
  return false;
}

///////////////////////////////////////////////////////////////
reco::Vertex
VtxAssociatorsUtils::findPrimaryVertex(const reco::CompositeCandidate& ZCand,edm::View<reco::Vertex> vertices){
  
  //look for closest vertex
  //vertices are already filtered... no need to check quality
  
  const reco::Candidate* lepton1 = ZCand.daughter(0);
  const reco::Candidate* lepton2 = ZCand.daughter(1);
  
  Double_t vz = (lepton1->vz()+lepton2->vz())/2.;
  Double_t mindz = 1000;
  reco::Vertex closestVertex = vertices.at(0);
  unsigned int vertexCollectionSize = vertices.size();  
  for (unsigned int i=0; i<vertexCollectionSize; i++) {
    reco::Vertex vertex = vertices.at(i);
    if( fabs(vertex.z()-vz)<mindz ){
      mindz = fabs(vertex.z()-vz);
      closestVertex = vertex;
    }
  }
  return closestVertex;
}

///////////////////////////////////////////////////////////////
Bool_t 
VtxAssociatorsUtils::loosezVertex(const reco::CompositeCandidate& ZCand, Double_t cut){
 
  const reco::Candidate* lepton1 = ZCand.daughter(0);
  const reco::Candidate* lepton2 = ZCand.daughter(1);
 
  //loose criteria: both leptons close one of each other
  return fabs(lepton1->vz()-lepton2->vz())<cut;
}

///////////////////////////////////////////////////////////////
//Candidate leptons version
Bool_t 
VtxAssociatorsUtils::tightzVertex(const reco::CompositeCandidate& ZCand, Double_t cut,const reco::Vertex& vertex){

  const reco::Candidate* lepton1 = ZCand.daughter(0);
  const reco::Candidate* lepton2 = ZCand.daughter(1);

  //strict criteria: both leptons are close to the (same) vertex
  return (fabs(lepton1->vz()-vertex.z())<cut and fabs(lepton2->vz()-vertex.z())<cut);
}

///////////////////////////////////////////////////////////////
//Other leptons version
Bool_t 
VtxAssociatorsUtils::tightzVertexMu(const pat::Muon& muon, Double_t cut,const reco::Vertex& vertex){

  //strict criteria: lepton is close to the vertex
  reco::TrackRef theTrack = muon.track();
  const math::XYZPoint myVertex(vertex.position().x(),vertex.position().y(),vertex.position().z());
  return (fabs(theTrack->dz(myVertex))<cut);
}

Bool_t 
VtxAssociatorsUtils::tightzVertexEle(const pat::Electron& electron, Double_t cut,const reco::Vertex& vertex){

  //strict criteria: lepton is close to the vertex
  reco::GsfTrackRef theTrack = electron.gsfTrack();
  const math::XYZPoint myVertex(vertex.position().x(),vertex.position().y(),vertex.position().z());
  return (fabs(theTrack->dz(myVertex))<cut);
}


///////////////////////////////////////////////////////////////
Bool_t 
VtxAssociatorsUtils::jetVertex(const reco::Vertex& vertex, const pat::Jet& jet, Int_t algo, Double_t sigmaCut, Double_t fraction){

  Bool_t decision=false;

  if (algo==1){ 
    decision = jetVertex_1(vertex,jet,sigmaCut,fraction);
  } else if(algo==2){ 
    decision = jetVertex_2(vertex,jet,sigmaCut,fraction);
  } else if(algo==3){ 
    decision = jetVertex_3(vertex,jet,sigmaCut,fraction);
  } else if(algo==4){
    decision = jetVertex_2b(vertex,jet,fraction);
  } else {
    std::cout<<"Error: unrecognized jet-vertex matching algorithm"<<std::endl<<std::cout<<"Please type a number between 1-4"<<std::endl;
  }
  
  return decision;
}

///////////////////////////////////////////////////////////////
Bool_t 
VtxAssociatorsUtils::jetVertex_1(const reco::Vertex& vertex, const pat::Jet& jet, Double_t sigcut, Double_t etcut){
  
  Double_t ptsum = 0;
  UInt_t jetconstsize = jet.getPFConstituents().size();

  for(UInt_t i=0; i<jetconstsize; i++){

    reco::PFCandidatePtr jetPFConst = jet.getPFConstituent(i);

    if (!(jetPFConst.isAvailable()) || (jetPFConst.isNull())){
      continue;
    }
    if(jetPFConst->trackRef().isNull()){
      continue;
    }
    if (jetPFConst->muonRef().isNonnull() || jetPFConst->gsfTrackRef().isNonnull()  ){
      continue;
    }
    
    Double_t distance = (jetPFConst->vz() - vertex.z());
    Double_t error = pow( (pow(jetPFConst->trackRef()->dzError(),2) + pow(vertex.zError(),2)),0.5);
    //error = vertex.zError();
    Double_t sig = distance/error;
    if (fabs(sig)<sigcut){
      ptsum += jetPFConst->pt();
    }
  }

  if (ptsum/jet.et() > etcut){
      return true;
  } else{
    return false;
  }
}

///////////////////////////////////////////////////////////////
Bool_t 
VtxAssociatorsUtils::jetVertex_2(const reco::Vertex& vertex, const pat::Jet& jet, Double_t sigcut, Double_t ptcut){
  
  Double_t ptsum = 0;
  Double_t ptsumall = 0;
  UInt_t jetconstsize = jet.getPFConstituents().size();

  for(UInt_t i=0; i<jetconstsize; i++){

    reco::PFCandidatePtr jetPFConst = jet.getPFConstituent(i);

    if (!(jetPFConst.isAvailable()) || (jetPFConst.isNull())){
      continue;
    }
    if (jetPFConst->trackRef().isNull()){
      continue;
    }
    if (jetPFConst->muonRef().isNonnull() || jetPFConst->gsfTrackRef().isNonnull()  ){
      continue;
    }

    Double_t distance = (jetPFConst->vz() - vertex.z());
    Double_t error = pow( (pow(jetPFConst->trackRef()->dzError(),2) + pow(vertex.zError(),2)),0.5);
    //error = vertex.zError();
    Double_t sig = distance/error;
    if (fabs(sig)<sigcut){
      ptsum+= jetPFConst->pt();
    }
    ptsumall+= jetPFConst->pt();
  }
  
  if( ptsumall==0.){
    return false;
  } else {
    if (ptsum/ptsumall > ptcut){
      return true;
    } else {
      return false;
    }
  }
}

///////////////////////////////////////////////////////////////
Bool_t 
VtxAssociatorsUtils::jetVertex_3(const reco::Vertex& vertex, const pat::Jet& jet, Double_t sigcut, Double_t etcut){

 Double_t ptsumx = 0;
 Double_t ptsumy = 0;
 UInt_t jetconstsize = jet.getPFConstituents().size();

  for(UInt_t i=0; i<jetconstsize; i++){

    reco::PFCandidatePtr jetPFConst = jet.getPFConstituent(i);

    if (!(jetPFConst.isAvailable()) || (jetPFConst.isNull())){
      continue;
    }
    if (jetPFConst->trackRef().isNull()){
      continue;
    }
    if (jetPFConst->muonRef().isNonnull() || jetPFConst->gsfTrackRef().isNonnull()  ){
      continue;
    }

    Double_t distance = (jetPFConst->vz() - vertex.z());
    Double_t error = pow( (pow(jetPFConst->trackRef()->dzError(),2) + pow(vertex.zError(),2)),0.5);
    //error = vertex.zError();
    Double_t sig = distance/error;
    if( fabs(sig)<sigcut ){
      ptsumx += jetPFConst->px();
      ptsumy += jetPFConst->py();
    }
  }

  if( pow((pow(ptsumx,2) +pow(ptsumy,2)),0.5)/jet.et()  > etcut ){
    return true;
  } else {
    return false;
  }
}


Bool_t
VtxAssociatorsUtils::jetVertex_2b(const reco::Vertex& vertex,const pat::Jet& jet,Double_t ptcut){

  Double_t ptsum = 0.;
  Double_t ptsumall = 0.;
  UInt_t jetconstsize = jet.getPFConstituents().size();

  for(UInt_t i=0; i<jetconstsize; i++){

    reco::PFCandidatePtr jetPFConst = jet.getPFConstituent(i);
     
    if (!(jetPFConst.isAvailable()) || (jetPFConst.isNull())){
      continue;
    }
    if (jetPFConst->trackRef().isNull()){
      continue;
    }
    if (jetPFConst->muonRef().isNonnull() || jetPFConst->gsfTrackRef().isNonnull()  ){
      continue;
    }
    
    for( reco::Vertex::trackRef_iterator tk = vertex.tracks_begin(); tk<vertex.tracks_end(); tk++ ){
      if(tk->key() == jetPFConst->trackRef().key()) {
	ptsum += jetPFConst->pt();
      } // closes if tkref = pfcandidate 
    } // closes loop on vertex tracks
    ptsumall += jetPFConst->pt();
  } // closes loop on pfconstituents
  
  if( ptsumall==0){
    return false;
  } 
  
  if(ptsum/ptsumall > ptcut){
    return true;
  } else{
    return false;
  }  
}
