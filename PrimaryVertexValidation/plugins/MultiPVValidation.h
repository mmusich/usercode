#ifndef MultiPVValidation_h
#define MultiPVValidation_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TTree.h"

#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/RoadSearchCloud/interface/RoadSearchCloud.h"
#include "DataFormats/RoadSearchCloud/interface/RoadSearchCloudCollection.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"

#include "DataFormats/SiStripDigi/interface/SiStripRawDigi.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "DataFormats/Common/interface/EDProduct.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "RecoVertex/VertexPrimitives/interface/VertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFindingBase.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackClusterizerInZ.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerAlgorithm.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

// system include files
#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <map>

//
// class decleration
//

class MultiPVValidation : public edm::EDAnalyzer {

 public:
  explicit MultiPVValidation(const edm::ParameterSet&);
  ~MultiPVValidation();

 private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  bool isHit2D(const TrackingRecHit &hit) const;
  bool hasFirstLayerPixelHits(const reco::TransientTrack track);

  // ----------member data ---------------------------
  edm::ParameterSet theConfig;
  edm::InputTag  TrackCollectionTag_;
  bool debug_;
  int Nevt_;
  //  TrackFilterForPVFinding theTrackFilter_;
  TrackFilterForPVFindingBase* theTrackFilter_; 
  TrackClusterizerInZ* theTrackClusterizer_;

  // Output
  bool lightNtupleSwitch_;   // switch to keep only info for daily validation
  std::string filename_;     // output file name 
  TFile* rootFile_;          
  TTree* rootTree_;
  
  // Root-Tuple variables :
  //=======================
  void SetVarToZero();  

  static const int nMaxtracks_ = 1000;

  // event-related quantities
  int nTracks_;
  int nTracksPerClus_;
  int nClus_;
  int nOfflineVertices_;
  unsigned int RunNumber_;
  unsigned int LuminosityBlockNumber_;
  double  xOfflineVertex_;
  double  yOfflineVertex_;
  double  zOfflineVertex_;
  double BSx0_;
  double BSy0_;
  double BSz0_;
  double Beamsigmaz_;
  double Beamdxdz_;   
  double BeamWidthX_;
  double BeamWidthY_;

  // track-related quantities
  double pt_[nMaxtracks_];   
  double p_[nMaxtracks_];    
  int nhits_[nMaxtracks_];
  int nhits1D_[nMaxtracks_];
  int nhits2D_[nMaxtracks_];
  int nhitsBPIX_[nMaxtracks_]; 
  int nhitsFPIX_[nMaxtracks_];
  int nhitsTIB_[nMaxtracks_];
  int nhitsTID_[nMaxtracks_];
  int nhitsTOB_[nMaxtracks_];
  int nhitsTEC_[nMaxtracks_];       
  double eta_[nMaxtracks_];
  double theta_[nMaxtracks_];
  double phi_[nMaxtracks_];
  double chi2_[nMaxtracks_];
  double chi2ndof_[nMaxtracks_];
  int    charge_[nMaxtracks_];
  double qoverp_[nMaxtracks_];
  double dz_[nMaxtracks_];
  double dxy_[nMaxtracks_];
  double dxyBs_[nMaxtracks_]; 
  double dzBs_[nMaxtracks_];  
  double xPCA_[nMaxtracks_];
  double yPCA_[nMaxtracks_];
  double zPCA_[nMaxtracks_];
  double xUnbiasedVertex_[nMaxtracks_];
  double yUnbiasedVertex_[nMaxtracks_];
  double zUnbiasedVertex_[nMaxtracks_];
  float  chi2normUnbiasedVertex_[nMaxtracks_]; 
  float  chi2UnbiasedVertex_[nMaxtracks_];
  float  DOFUnbiasedVertex_[nMaxtracks_];
  float  sumOfWeightsUnbiasedVertex_[nMaxtracks_];
  int    tracksUsedForVertexing_[nMaxtracks_];
  double dxyFromMyVertex_[nMaxtracks_];
  double dzFromMyVertex_[nMaxtracks_];
  double dszFromMyVertex_[nMaxtracks_];
  int   hasRecVertex_[nMaxtracks_];
  int   isGoodTrack_[nMaxtracks_];

};

#endif
