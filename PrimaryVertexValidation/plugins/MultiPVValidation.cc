// -*- C++ -*-
//
// Package:    MultiPVValidation
// Class:      MultiPVValidation
// 
/**\class MultiPVValidation MultiPVValidation.cc Alignment/OfflineValidation/plugins/MultiPVValidation.cc

 Description: Validate alignment constants using unbiased vertex residuals

 Implementation:
 <Notes on implementation>
*/
//
// Original Author:  Marco Musich
//         Created:  Tue Mar 02 10:39:34 CDT 2010
//

// system include files
#include <memory>


// user include files
#include "Alignment/OfflineValidation/plugins/MultiPVValidation.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "FWCore/Framework/interface/Event.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h>
#include <SimDataFormats/TrackingHit/interface/PSimHit.h>

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include "TNtuple.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include <Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h>
#include <DataFormats/GeometrySurface/interface/Surface.h>
#include <DataFormats/GeometrySurface/interface/GloballyPositioned.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include "MagneticField/Engine/interface/MagneticField.h" 
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/PrimaryVertexProducer/interface/GapClusterizerInZ.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ.h"

 const int kBPIX = PixelSubdetector::PixelBarrel;
 const int kFPIX = PixelSubdetector::PixelEndcap;

// Constructor

MultiPVValidation::MultiPVValidation(const edm::ParameterSet& iConfig)
  : theConfig(iConfig), 
    theTrackFilter_(iConfig.getParameter<edm::ParameterSet>("TkFilterParameters"))
{
  //now do what ever initialization is needed
  debug_    = iConfig.getParameter<bool>       ("Debug");  
  TrackCollectionTag_      = iConfig.getParameter<edm::InputTag>("TrackCollectionTag");  
  filename_ = iConfig.getParameter<std::string>("OutputFileName"); 
  // theTrackClusterizer_ = new GapClusterizerInZ(iConfig.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkGapClusParameters"));
  
  // select and configure the track clusterizer  
  std::string clusteringAlgorithm=iConfig.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<std::string>("algorithm");
  if (clusteringAlgorithm=="gap"){
    theTrackClusterizer_ = new GapClusterizerInZ(iConfig.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkGapClusParameters"));
  }else if(clusteringAlgorithm=="DA"){
    theTrackClusterizer_ = new DAClusterizerInZ(iConfig.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
  }else{
    throw VertexException("PrimaryVertexProducerAlgorithm: unknown clustering algorithm: " + clusteringAlgorithm);  
  }
}
   
// Destructor
MultiPVValidation::~MultiPVValidation()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
MultiPVValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace std;
  
  Nevt_++;  

  //=======================================================
  // Initialize Root-tuple variables
  //=======================================================

  SetVarToZero();
 
  //=======================================================
  // Retrieve the Magnetic Field information
  //=======================================================

  edm::ESHandle<MagneticField> theMGField;
  iSetup.get<IdealMagneticFieldRecord>().get( theMGField );

  //=======================================================
  // Retrieve the Tracking Geometry information
  //=======================================================

  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get( theTrackingGeometry );

  //=======================================================
  // Retrieve the Track information
  //=======================================================
  
  edm::Handle<reco::TrackCollection>  trackCollectionHandle;
  iEvent.getByLabel(TrackCollectionTag_, trackCollectionHandle);
  
  //=======================================================
  // Retrieve offline vartex information (only for reco)
  //=======================================================
 
  edm::Handle<reco::VertexCollection> vertices;
  try {
    iEvent.getByLabel("offlinePrimaryVertices", vertices);
  } catch (...) {
    if(debug_)
      cout << "No offlinePrimaryVertices found!" << endl;
  }
  if ( vertices.isValid() ) {
    xOfflineVertex_ = (*vertices)[0].x();
    yOfflineVertex_ = (*vertices)[0].y();
    zOfflineVertex_ = (*vertices)[0].z();
  }

  nOfflineVertices_ = vertices.product()->size();
 
  //=======================================================
  // Retrieve Beamspot information
  //=======================================================

  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
    
  if ( beamSpotHandle.isValid() )
    {
      beamSpot = *beamSpotHandle;
      BSx0_ = beamSpot.x0();
      BSy0_ = beamSpot.y0();
      BSz0_ = beamSpot.z0();
      Beamsigmaz_ = beamSpot.sigmaZ();    
      Beamdxdz_ = beamSpot.dxdz();	     
      BeamWidthX_ = beamSpot.BeamWidthX();
      BeamWidthY_ = beamSpot.BeamWidthY();
    } else
    {
      if(debug_)
	cout << "No BeamSpot found!" << endl;
    }
  
  if(debug_)
    std::cout<<"Beamspot x:"<<BSx0_<<" y:"<<BSy0_<<" z:"<<BSz0_<<std::endl; 
  
  //double sigmaz = beamSpot.sigmaZ();
  //double dxdz = beamSpot.dxdz();
  //double BeamWidth = beamSpot.BeamWidth();
  
  //=======================================================
  // Starts here ananlysis
  //=======================================================

  RunNumber_ = iEvent.id().run();
  LuminosityBlockNumber_ = iEvent.id().luminosityBlock();

  if(debug_)
    std::cout<<"MultiPVValidation::analyze() looping over "<<trackCollectionHandle->size()<< "tracks." <<std::endl;       

  //======================================================
  // Interface RECO tracks to vertex reconstruction
  //======================================================

  //edm::ESHandle<TransientTrackBuilder> theB;
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  //vector<reco::TransientTrack> t_tks = (*theB).build(trackCollectionHandle,beamSpot);
  
  std::vector<reco::TransientTrack> t_tks;
  unsigned int k = 0;   
  for(reco::TrackCollection::const_iterator track = trackCollectionHandle->begin(); track!= trackCollectionHandle->end(); ++track, ++k){
  
    reco::TrackRef trackref(trackCollectionHandle,k);
    reco::TransientTrack theTTRef = reco::TransientTrack(trackref, &*theMGField, theTrackingGeometry );
    t_tks.push_back(theTTRef);
  
  }
  
  if(debug_) {cout << "PrimaryVertexValidation"
		  << "Found: " << t_tks.size() << " reconstructed tracks" << "\n";
  }

  //======================================================
  // clusterize tracks in Z
  //======================================================

  vector< vector<reco::TransientTrack> > clusters = theTrackClusterizer_->clusterize(t_tks);
  
  if (debug_){
    cout <<  " clustering returned  "<< clusters.size() << " clusters  from " << t_tks.size() << " selected tracks" <<endl;
  }
  
  nClus_=clusters.size();  

  //======================================================
  // Starts loop on clusters 
  //======================================================

  for (vector< vector<reco::TransientTrack> >::const_iterator iclus = clusters.begin(); iclus != clusters.end(); iclus++) {

    nTracksPerClus_=0;

    unsigned int i = 0;   
    for(vector<reco::TransientTrack>::const_iterator theTrack = iclus->begin(); theTrack!= iclus->end(); ++theTrack, ++i)
      {
	if ( nTracks_ >= nMaxtracks_ ) {
	  std::cout << " MultiPVValidation::analyze() : Warning - Number of tracks: " << nTracks_ << " , greater than " << nMaxtracks_ << std::endl;
	  continue;
	}
	
	pt_[nTracks_]       = theTrack->track().pt();
	p_[nTracks_]        = theTrack->track().p();
	nhits_[nTracks_]    = theTrack->track().numberOfValidHits();
	eta_[nTracks_]      = theTrack->track().eta();
	theta_[nTracks_]    = theTrack->track().theta();
	phi_[nTracks_]      = theTrack->track().phi();
	chi2_[nTracks_]     = theTrack->track().chi2();
	chi2ndof_[nTracks_] = theTrack->track().normalizedChi2();
	charge_[nTracks_]   = theTrack->track().charge();
	qoverp_[nTracks_]   = theTrack->track().qoverp();
	dz_[nTracks_]       = theTrack->track().dz();
	dxy_[nTracks_]      = theTrack->track().dxy();
	
	math::XYZPoint point(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());
	dxyBs_[nTracks_]    = theTrack->track().dxy(point);
	dzBs_[nTracks_]     = theTrack->track().dz(point);

	xPCA_[nTracks_]     = theTrack->track().vertex().x();
	yPCA_[nTracks_]     = theTrack->track().vertex().y();
	zPCA_[nTracks_]     = theTrack->track().vertex().z(); 
	
	//=======================================================
	// Retrieve rechit information
	//=======================================================  
	
	int nRecHit1D =0;
	int nRecHit2D =0;
	int nhitinTIB =0; 
	int nhitinTOB =0; 
	int nhitinTID =0; 
	int nhitinTEC =0; 
	int nhitinBPIX=0;
	int nhitinFPIX=0; 
	
	for (trackingRecHit_iterator iHit = theTrack->recHitsBegin(); iHit != theTrack->recHitsEnd(); ++iHit) {
	  if((*iHit)->isValid()) {	
	    
	    if (this->isHit2D(**iHit)) {++nRecHit2D;}
	    else {++nRecHit1D; }
	    
	    int type =(*iHit)->geographicalId().subdetId();
	    
	    if(type==int(StripSubdetector::TIB)){++nhitinTIB;}
	    if(type==int(StripSubdetector::TOB)){++nhitinTOB;}
	    if(type==int(StripSubdetector::TID)){++nhitinTID;}
	    if(type==int(StripSubdetector::TEC)){++nhitinTEC;}
	    if(type==int(                kBPIX)){++nhitinBPIX;}
	    if(type==int(                kFPIX)){++nhitinFPIX;}
	  }
	}      

	nhits1D_[nTracks_]     =nRecHit1D;
	nhits2D_[nTracks_]     =nRecHit2D;
	nhitsBPIX_[nTracks_]   =nhitinBPIX;
	nhitsFPIX_[nTracks_]   =nhitinFPIX;
	nhitsTIB_[nTracks_]    =nhitinTIB;
	nhitsTID_[nTracks_]    =nhitinTID;
	nhitsTOB_[nTracks_]    =nhitinTOB;
	nhitsTEC_[nTracks_]    =nhitinTEC;
	
	//=======================================================
	// Good tracks for vertexing selection
	//=======================================================  

	bool hasTheProbeFirstPixelLayerHit = false;
	hasTheProbeFirstPixelLayerHit = this->hasFirstLayerPixelHits((*theTrack));
	if (theTrackFilter_((*theTrack))&&hasTheProbeFirstPixelLayerHit){
	  isGoodTrack_[nTracks_]=1;
	}
      
	//=======================================================
	// Fit unbiased vertex
	//=======================================================
	
	vector<reco::TransientTrack> theFinalTracks;
	theFinalTracks.clear();

	for(vector<reco::TransientTrack>::const_iterator tk = iclus->begin(); tk!= iclus->end(); ++tk){
	  
	  hasTheProbeFirstPixelLayerHit = this->hasFirstLayerPixelHits((*tk));
	  if (theTrackFilter_((*tk))&&hasTheProbeFirstPixelLayerHit){
	    if( tk == theTrack ) continue;
	    else {
	      theFinalTracks.push_back((*tk));
	    }
	  }
	}
	
	if(theFinalTracks.size() > 2){
	    
	  if(debug_)
	    std::cout <<"MultiPVValidation::analyze() :Transient Track Collection size: "<<theFinalTracks.size()<<std::endl;
	  
	  try{
	      
	    VertexFitter<5>* fitter = new AdaptiveVertexFitter;
	    TransientVertex theFittedVertex = fitter->vertex(theFinalTracks);
	    
	    if(theFittedVertex.isValid ()){
	      
	      if(theFittedVertex.hasTrackWeight()){
		for(size_t rtracks= 0; rtracks < theFinalTracks.size(); rtracks++){
		  sumOfWeightsUnbiasedVertex_[nTracks_] += theFittedVertex.trackWeight(theFinalTracks[rtracks]);
		}
	      }
	      
	      const math::XYZPoint myVertex(theFittedVertex.position().x(),theFittedVertex.position().y(),theFittedVertex.position().z());
	      hasRecVertex_[nTracks_]    = 1;
	      xUnbiasedVertex_[nTracks_] = theFittedVertex.position().x();
	      yUnbiasedVertex_[nTracks_] = theFittedVertex.position().y();
	      zUnbiasedVertex_[nTracks_] = theFittedVertex.position().z();
	      chi2normUnbiasedVertex_[nTracks_] = theFittedVertex.normalisedChiSquared();
	      chi2UnbiasedVertex_[nTracks_] = theFittedVertex.totalChiSquared();
	      DOFUnbiasedVertex_[nTracks_] = theFittedVertex.degreesOfFreedom();   
	      tracksUsedForVertexing_[nTracks_] = theFinalTracks.size();
	      dxyFromMyVertex_[nTracks_] = theTrack->track().dxy(myVertex);
	      dzFromMyVertex_[nTracks_]  = theTrack->track().dz(myVertex);
	      dszFromMyVertex_[nTracks_] = theTrack->track().dsz(myVertex);
	      
	      if(debug_){
		std::cout<<"MultiPVValidation::analyze() :myVertex.x()= "<<myVertex.x()<<" myVertex.y()= "<<myVertex.y()<<" theFittedVertex.z()= "<<myVertex.z()<<std::endl;	  
		std::cout<<"MultiPVValidation::analyze() : theTrack->track().dz(myVertex)= "<<theTrack->track().dz(myVertex)<<std::endl;
		std::cout<<"MultiPVValidation::analyze() : zPCA -myVertex.z() = "<<(theTrack->track().vertex().z() -myVertex.z() )<<std::endl; 
	      }// ends if debug_
	    } // ends if the fitted vertex is Valid
	    
	  }  catch ( cms::Exception& er ) {
	    LogTrace("MultiPVValidation::analyze RECO")<<"caught std::exception "<<er.what()<<std::endl;
	  }
	} //ends if theFinalTracks.size() > 2
	
	else {
	  if(debug_)
	      std::cout << "MultiPVValidation::analyze() :Not enough tracks to make a vertex.  Returns no vertex info" << std::endl;
	}
	  
	++nTracks_;  
	++nTracksPerClus_;

	if(debug_)
	  cout<< "Track "<<i<<" : pT = "<<theTrack->track().pt()<<endl;
	
      }// for loop on tracks

  } // for loop on track clusters
  
  rootTree_->Fill();
  
} 

// ------------ method called to discriminate 1D from 2D hits  ------------
bool MultiPVValidation::isHit2D(const TrackingRecHit &hit) const
{
  if (hit.dimension() < 2) {
    return false; // some (muon...) stuff really has RecHit1D
  } else {
    const DetId detId(hit.geographicalId());
    if (detId.det() == DetId::Tracker) {
      if (detId.subdetId() == kBPIX || detId.subdetId() == kFPIX) {
        return true; // pixel is always 2D
      } else { // should be SiStrip now
        if (dynamic_cast<const SiStripRecHit2D*>(&hit)) return false; // normal hit
        else if (dynamic_cast<const SiStripMatchedRecHit2D*>(&hit)) return true; // matched is 2D
        else if (dynamic_cast<const ProjectedSiStripRecHit2D*>(&hit)) return false; // crazy hit...
        else {
          edm::LogError("UnkownType") << "@SUB=AlignmentTrackSelector::isHit2D"
                                      << "Tracker hit not in pixel and neither SiStripRecHit2D nor "
                                      << "SiStripMatchedRecHit2D nor ProjectedSiStripRecHit2D.";
          return false;
        }
      }
    } else { // not tracker??
      edm::LogWarning("DetectorMismatch") << "@SUB=AlignmentTrackSelector::isHit2D"
                                          << "Hit not in tracker with 'official' dimension >=2.";
      return true; // dimension() >= 2 so accept that...
    }
  }
  // never reached...
}

// ------------ method to check the presence of pixel hits  ------------
bool MultiPVValidation::hasFirstLayerPixelHits(const reco::TransientTrack track)
{
  bool accepted = false;
  // hit pattern of the track
  const reco::HitPattern& p = track.hitPattern();      
  for (int i=0; i<p.numberOfHits(); i++) {
    uint32_t pattern = p.getHitPattern(i);   
    if (p.pixelBarrelHitFilter(pattern) || p.pixelEndcapHitFilter(pattern) ) {
      if (p.getLayer(pattern) == 1) {
	if (p.validHitFilter(pattern)) {
	  accepted = true;
	}
      }
    }
  }
  return accepted;
} 

// ------------ method called once each job before begining the event loop  ------------
void MultiPVValidation::beginJob()
{
  edm::LogInfo("beginJob") << "Begin Job" << std::endl;
  // Define TTree for output

  Nevt_    = 0;
  
  rootFile_ = new TFile(filename_.c_str(),"recreate");
  rootTree_ = new TTree("tree","PV Validation tree");
  
  // Track Paramters 
  rootTree_->Branch("nTracks",&nTracks_,"nTracks/I");
  rootTree_->Branch("nTracksPerClus",&nTracksPerClus_,"nTracksPerClus/I");
  rootTree_->Branch("nClus",&nClus_,"nClus/I");
  rootTree_->Branch("nOfflineVertices",&nOfflineVertices_,"nOfflineVertices/I");
  rootTree_->Branch("nTracksPerOfflineVertex",&nTracksPerOfflineVertex_,"nTracksPerOfflineVertex/I");
  rootTree_->Branch("RunNumber",&RunNumber_,"RunNumber/i");
  rootTree_->Branch("LuminosityBlockNumber",&LuminosityBlockNumber_,"LuminosityBlockNumber/i");
  rootTree_->Branch("xOfflineVertex",&xOfflineVertex_,"xOfflineVertex/D");
  rootTree_->Branch("yOfflineVertex",&yOfflineVertex_,"yOfflineVertex/D");
  rootTree_->Branch("zOfflineVertex",&zOfflineVertex_,"zOfflineVertex/D");
  rootTree_->Branch("BSx0",&BSx0_,"BSx0/D");
  rootTree_->Branch("BSy0",&BSy0_,"BSy0/D");
  rootTree_->Branch("BSz0",&BSz0_,"BSz0/D");
  rootTree_->Branch("Beamsigmaz",&Beamsigmaz_,"Beamsigmaz/D");
  rootTree_->Branch("Beamdxdz",&Beamdxdz_,"Beamdxdz/D");
  rootTree_->Branch("BeamWidthX",&BeamWidthX_,"BeamWidthX/D");
  rootTree_->Branch("BeamWidthY",&BeamWidthY_,"BeamWidthY/D");
  rootTree_->Branch("pt",&pt_,"pt[nTracks]/D");
  rootTree_->Branch("p",&p_,"p[nTracks]/D");
  rootTree_->Branch("nhits",&nhits_,"nhits[nTracks]/I");
  rootTree_->Branch("nhits1D",&nhits1D_,"nhits1D[nTracks]/I");
  rootTree_->Branch("nhits2D",&nhits2D_,"nhits2D[nTracks]/I");
  rootTree_->Branch("nhitsBPIX",&nhitsBPIX_,"nhitsBPIX[nTracks]/I");
  rootTree_->Branch("nhitsFPIX",&nhitsFPIX_,"nhitsFPIX[nTracks]/I");
  rootTree_->Branch("nhitsTIB",&nhitsTIB_,"nhitsTIB[nTracks]/I");
  rootTree_->Branch("nhitsTID",&nhitsTID_,"nhitsTID[nTracks]/I");
  rootTree_->Branch("nhitsTOB",&nhitsTOB_,"nhitsTOB[nTracks]/I");
  rootTree_->Branch("nhitsTEC",&nhitsTEC_,"nhitsTEC[nTracks]/I");
  rootTree_->Branch("eta",&eta_,"eta[nTracks]/D");
  rootTree_->Branch("theta",&theta_,"theta[nTracks]/D");
  rootTree_->Branch("phi",&phi_,"phi[nTracks]/D");
  rootTree_->Branch("chi2",&chi2_,"chi2[nTracks]/D");
  rootTree_->Branch("chi2ndof",&chi2ndof_,"chi2ndof[nTracks]/D");
  rootTree_->Branch("charge",&charge_,"charge[nTracks]/I");
  rootTree_->Branch("qoverp",&qoverp_,"qoverp[nTracks]/D");
  rootTree_->Branch("dz",&dz_,"dz[nTracks]/D");
  rootTree_->Branch("dxy",&dxy_,"dxy[nTracks]/D");
  rootTree_->Branch("dzBs",&dzBs_,"dzBs[nTracks]/D");
  rootTree_->Branch("dxyBs",&dxyBs_,"dxyBs[nTracks]/D");
  rootTree_->Branch("xPCA",&xPCA_,"xPCA[nTracks]/D");
  rootTree_->Branch("yPCA",&yPCA_,"yPCA[nTracks]/D");
  rootTree_->Branch("zPCA",&zPCA_,"zPCA[nTracks]/D");
  rootTree_->Branch("xUnbiasedVertex",&xUnbiasedVertex_,"xUnbiasedVertex[nTracks]/D");
  rootTree_->Branch("yUnbiasedVertex",&yUnbiasedVertex_,"yUnbiasedVertex[nTracks]/D");
  rootTree_->Branch("zUnbiasedVertex",&zUnbiasedVertex_,"zUnbiasedVertex[nTracks]/D");
  rootTree_->Branch("chi2normUnbiasedVertex",&chi2normUnbiasedVertex_,"chi2normUnbiasedVertex[nTracks]/F");
  rootTree_->Branch("chi2UnbiasedVertex",&chi2UnbiasedVertex_,"chi2UnbiasedVertex[nTracks]/F");
  rootTree_->Branch("DOFUnbiasedVertex",&DOFUnbiasedVertex_," DOFUnbiasedVertex[nTracks]/F");
  rootTree_->Branch("sumOfWeightsUnbiasedVertex",&sumOfWeightsUnbiasedVertex_,"sumOfWeightsUnbiasedVertex[nTracks]/F");
  rootTree_->Branch("tracksUsedForVertexing",&tracksUsedForVertexing_,"tracksUsedForVertexing[nTracks]/I");
  rootTree_->Branch("dxyFromMyVertex",&dxyFromMyVertex_,"dxyFromMyVertex[nTracks]/D");
  rootTree_->Branch("dzFromMyVertex",&dzFromMyVertex_,"dzFromMyVertex[nTracks]/D");
  rootTree_->Branch("dszFromMyVertex",&dszFromMyVertex_,"dszFromMyVertex[nTracks]/D");
  rootTree_->Branch("hasRecVertex",&hasRecVertex_,"hasRecVertex[nTracks]/I");
  rootTree_->Branch("isGoodTrack",&isGoodTrack_,"isGoodTrack[nTracks]/I");
}

// ------------ method called once each job just after ending the event loop  ------------
void MultiPVValidation::endJob() 
{

  std::cout<<"######################################"<<std::endl;
  std::cout<<"Number of analyzed events: "<<Nevt_<<std::endl;
  std::cout<<"######################################"<<std::endl;
  
   if ( rootFile_ ) {
     rootFile_->Write();
     rootFile_->Close();
   }
}

void MultiPVValidation::SetVarToZero() {
  
  nTracks_ = 0;
  nClus_ = 0;
  nOfflineVertices_=0;
  RunNumber_ =0;
  LuminosityBlockNumber_=0;
  xOfflineVertex_ =-999.;
  yOfflineVertex_ =-999.;
  zOfflineVertex_ =-999.;
  BSx0_ = -999.;
  BSy0_ = -999.;
  BSz0_ = -999.;
  Beamsigmaz_=-999.;
  Beamdxdz_=-999.;   
  BeamWidthX_=-999.;
  BeamWidthY_=-999.;

  for ( int i=0; i<nMaxtracks_; ++i ) {
    pt_[i]        = 0;
    p_[i]         = 0;
    nhits_[i]     = 0;
    nhits1D_[i]   = 0;
    nhits2D_[i]   = 0;
    nhitsBPIX_[i]  = 0;
    nhitsFPIX_[i]  = 0;
    nhitsTIB_[i]   = 0;
    nhitsTID_[i]   = 0;
    nhitsTOB_[i]   = 0;
    nhitsTEC_[i]   = 0;
    eta_[i]       = 0;
    theta_[i]       = 0;
    phi_[i]       = 0;
    chi2_[i]      = 0;
    chi2ndof_[i]  = 0;
    charge_[i]    = 0;
    qoverp_[i]    = 0;
    dz_[i]        = 0;
    dxy_[i]       = 0;
    dzBs_[i]      = 0;
    dxyBs_[i]     = 0;
    xPCA_[i]      = 0;
    yPCA_[i]      = 0;
    zPCA_[i]      = 0;
    xUnbiasedVertex_[i] =0;    
    yUnbiasedVertex_[i] =0;
    zUnbiasedVertex_[i] =0;
    chi2normUnbiasedVertex_[i]=0;
    chi2UnbiasedVertex_[i]=0;
    DOFUnbiasedVertex_[i]=0;
    sumOfWeightsUnbiasedVertex_[i]=0;
    tracksUsedForVertexing_[i]=0;
    dxyFromMyVertex_[i]=0;
    dzFromMyVertex_[i]=0;
    dszFromMyVertex_[i]=0;
    hasRecVertex_[i] = 0;
    isGoodTrack_[i]  = 0;
  } 
}

//define this as a plug-in
DEFINE_FWK_MODULE(MultiPVValidation);
