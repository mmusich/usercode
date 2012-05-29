#include "ZbbAnalysis/AnalysisStep/interface/VtxAssociatorsUtils.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <list>


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

///////////////////////////////////////////////////////////////
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

std::map<TString,Double_t> 
VtxAssociatorsUtils::calculateJetVertexAssociation(const pat::Jet& jet, reco::Vertex vertex){

  std::map<TString,Double_t> theJVFmap;

  Double_t ptsum = 0.;
  Double_t ptsumx = 0.;
  Double_t ptsumy = 0.;
  Double_t ptsumall = 0.;
  Double_t ratio1 = 0;
  Double_t ratio2 = 0;
  Double_t ratio3 = 0;         
  
  Double_t ratio1b = 0;	      
  Double_t ratio2b = 0;	      
  Double_t ratio3b = 0;        

  Double_t sigcut = 5;

  UInt_t jet_constsize = jet.getPFConstituents().size();

  for(UInt_t j=0; j< jet_constsize; j++){
      
    reco::PFCandidatePtr jetPFConst = jet.getPFConstituent(j);
    
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
    Double_t error = pow( (pow((jetPFConst->trackRef()->dzError()),2) + pow(vertex.zError(),2)),0.5);			
    Double_t sig = distance/error;
    
    if( fabs(sig)<sigcut ){
      ptsum  += jetPFConst->pt();
      ptsumx += jetPFConst->px();
      ptsumy += jetPFConst->py();
    }
    ptsumall += jetPFConst->pt();
  }  
  
  // calculate the Jet Vertex fractions
  ratio1 = ptsum/jet.et();    
 
  if(ptsumall>0) {
    ratio2 = ptsum/ptsumall;
  } else{
    ratio2 = -1;
  }
 
  ratio3 =  pow((pow(ptsumx,2)+pow(ptsumy,2)),0.5)/jet.et();

  // fill the map
  theJVFmap["1"]=ratio1;
  theJVFmap["2"]=ratio2;
  theJVFmap["3"]=ratio3;
  
  // reset the counters to zero
  ptsum = 0.;	  
  ptsumx = 0.;	  
  ptsumy = 0.;	  
  ptsumall = 0.; 

  for(UInt_t j=0; j< jet_constsize; j++){
    
    reco::PFCandidatePtr jetPFConst = jet.getPFConstituent(j);
    
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

  // calculate the Jet Vertex fractions
  ratio1b = ptsum/jet.et();    
   
  if(ptsumall>0) {
    ratio2b = ptsum/ptsumall;
  } else{
    ratio2b = -1;
  }
 
  ratio3b =  pow((pow(ptsumx,2)+pow(ptsumy,2)),0.5)/jet.et();
  
  // fill the map
  theJVFmap["1b"]=ratio1b;
  theJVFmap["2b"]=ratio2b;
  theJVFmap["3b"]=ratio3b;
  
  return theJVFmap;

}


// ------------ method called to build the reference to track refs -------------
std::list<int>
VtxAssociatorsUtils::buildTrackRefs(edm::View<reco::Vertex> vertices,bool isPU){
  
  if(!isPU){
    std::list<int> _trackrefs_PV;
    // loop over the tracks making the PV, and store the keys
    _trackrefs_PV.clear();
    if(vertices.size()>0) {
      const reco::Vertex pv = vertices.operator[](0);
      for( reco::Vertex::trackRef_iterator tk = pv.tracks_begin(); tk < pv.tracks_end(); ++tk) {
	_trackrefs_PV.push_back(tk->key());
      }
    }
    return _trackrefs_PV;
  } else {
    std::list<int> _trackrefs_PU;
    // loop over the tracks making the PU vertices, and store the keys
    _trackrefs_PU.clear();
    if(vertices.size()>1) {
      for(edm::View<reco::Vertex>::const_iterator PUvertex = vertices.begin()+1; PUvertex<vertices.end(); ++PUvertex) {
	for( reco::Vertex::trackRef_iterator tk = PUvertex->tracks_begin(); tk < PUvertex->tracks_end(); ++tk) {
	  _trackrefs_PU.push_back(tk->key());
	}
      }
    }
    return _trackrefs_PU;
  }
}

// ------------ method called to produce the beta  ------------
float
VtxAssociatorsUtils::beta(pat::Jet const& jet,std::list<int> _trackrefs_PV){
  // definition of beta: ratio of charged pT from first vertex over the sum of all the charged pT in jet. 
  // by definition, 0 if there is no primary vertex.
  if(_trackrefs_PV.empty()) return 0.;
  // now loop over the jet charged constituents, and count pt
  float ptsum = 0.;
  float ptsumall = 0.;
  const std::vector< reco::PFCandidatePtr > constituents = jet.getPFConstituents();
  for(std::vector< reco::PFCandidatePtr >::const_iterator jetconst = constituents.begin(); jetconst < constituents.end(); ++jetconst) {
    if((*jetconst)->trackRef().isNull()) continue;
    std::list<int>::const_iterator found = find(_trackrefs_PV.begin(), _trackrefs_PV.end(), (*jetconst)->trackRef().key());
    if(found!=_trackrefs_PV.end()) {
      ptsum += (*jetconst)->pt();
    }
    ptsumall += (*jetconst)->pt();
   
  }
  if(ptsumall<0.001) // non-null: 0.001 is much lower than any pt cut at reco level.
    return -1.;
  return ptsum/ptsumall;
}

// ------------ method called to produce the beta*  -----------
float
VtxAssociatorsUtils::betaStar(pat::Jet const& jet, std::list<int> _trackrefs_PU){
  // defined as the ratio of charged pT coming from good PU vertices over the sum of all charged pT in jet. 
  // by definition, 0 if there is no PU vertex.
  if(_trackrefs_PU.empty()) return 0.;
  // now loop over the jet charged constituents, and count pt
  float ptsum = 0.;
  float ptsumall = 0.;
  const std::vector< reco::PFCandidatePtr > constituents = jet.getPFConstituents();
  for(std::vector< reco::PFCandidatePtr >::const_iterator jetconst = constituents.begin(); jetconst < constituents.end(); ++jetconst) {
    if((*jetconst)->trackRef().isNull()) continue;
    std::list<int>::const_iterator found = find(_trackrefs_PU.begin(), _trackrefs_PU.end(), (*jetconst)->trackRef().key());   
    if(found!=_trackrefs_PU.end()) ptsum += (*jetconst)->pt();
    ptsumall += (*jetconst)->pt();
  }
  if(ptsumall<0.001) // non-null: 0.001 is much lower than any pt cut at reco level.
    return -1.;
  return ptsum/ptsumall;
}
