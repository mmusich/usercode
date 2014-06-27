/***********************************************************
*
*  Class PixelStdHitNtuplizer
*
*  - v0.1: Obtained merging PixelTree and StdHitNtuplizer
*
*  Authors: E. Migliore, M. Musich 
*
************************************************************/

// DataFormats
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h" 
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h" 
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Geometry
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"

// Structures needed
#include "AuxCode/SLHCSimPhase2/interface/RecHit.h"
#include "AuxCode/SLHCSimPhase2/interface/evt.h"

// For ROOT
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>

// STD
#include <memory>
#include <string>
#include <iostream>

using namespace std;
using namespace edm;
using namespace reco;
  
class PixelStdHitNtuplizer : public edm::EDAnalyzer
{
public:
   
  explicit PixelStdHitNtuplizer(const edm::ParameterSet& conf);
  virtual ~PixelStdHitNtuplizer();
  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const edm::Event& e, const edm::EventSetup& es);
 
protected:
 
  void fillEvt(const edm::Event& );
  //void fillPRecHit(const int subid, SiPixelRecHitCollection::const_iterator pixeliter,
  //                 const GeomDet* PixGeom);
  void fillPRecHit(const int detid_db, const int subid, 
		   const int layer_num,const int ladder_num,const int module_num,
		   const int disk_num,const int blade_num,const int panel_num,const int side_num,
		   SiPixelRecHitCollection::DetSet::const_iterator pixeliter,
		   const int num_simhit,
		   std::vector<PSimHit>::const_iterator closest_simhit,
		   const GeomDet* PixGeom);
  void fillPRecHit(const int detid_db, const int subid,
		   const int layer_num,const int ladder_num,const int module_num,
		   const int disk_num,const int blade_num,const int panel_num,const int side_num,
		   trackingRecHit_iterator pixeliter,
		   const int num_simhit,
		   std::vector<PSimHit>::const_iterator closest_simhit,
		   const GeomDet* PixGeom);
 
private:
  edm::ParameterSet conf_;
  edm::InputTag src_;
  edm::InputTag rphiRecHits_;
  edm::InputTag stereoRecHits_;
  edm::InputTag matchedRecHits_;
 
  static const int DIGIMAX = 1000000;
  static const int CLUSTERMAX = 100000;
  static const int DGPERCLMAX = 100;  
     
  struct RecHit recHit_;
  struct evt evt_;

  void init();

  TFile * tfile_;
  TTree * pixeltree_;
};
 
PixelStdHitNtuplizer::PixelStdHitNtuplizer(edm::ParameterSet const& conf) : 
  conf_(conf), 
  src_( conf.getParameter<edm::InputTag>( "src" ) ),
  rphiRecHits_( conf.getParameter<edm::InputTag>("rphiRecHits") ),
  stereoRecHits_( conf.getParameter<edm::InputTag>("stereoRecHits") ),
  matchedRecHits_( conf.getParameter<edm::InputTag>("matchedRecHits") ),
  tfile_(0), 
  pixeltree_(0)
{
}
 
 
PixelStdHitNtuplizer::~PixelStdHitNtuplizer() { }  
 
void PixelStdHitNtuplizer::endJob() 
{
  std::cout << " PixelStdHitNtuplizer::endJob" << std::endl;
  tfile_->Write();
  tfile_->Close();
}
 
void PixelStdHitNtuplizer::beginJob()
{
  std::cout << " PixelStdHitNtuplizer::beginJob" << std::endl;
  std::string outputFile = conf_.getParameter<std::string>("OutputFile");
  
  tfile_ = new TFile ( outputFile.c_str() , "RECREATE" );
  pixeltree_ = new TTree("PixelNtuple","Pixel hit analyzer ntuple");
 
  int bufsize = 64000;
  //Common Branch
  pixeltree_->Branch("evt",&evt_,"run/I:evtnum/I", bufsize);

  //Common Branch
  pixeltree_->Branch("pixel_recHit", &recHit_, 
		     "pdgid/I:process:q/F:x:y:xx:xy:yy:row:col:gx:gy:gz:subid/I:module:layer:ladder:disk:blade:panel:side:nsimhit:spreadx:spready:hx/F:hy:tx:ty:tz:theta:phi:DgN/I:DgRow[DgN]/I:DgCol[DgN]/I:DgDetId[DgN]/I:DgAdc[DgN]/F:DgCharge[DgN]/F:DgClI[DgN]/I:ClN/I:ClDgN[ClN]/I:ClDgI[ClN][100]/I]", bufsize);

  // Several branches
  //   pixeltree_->Branch("evt",     &evt_ ,          "run/I:evtnum/I", bufsize);
  //   pixeltree_->Branch("pdgid",   &recHit_.pdgid,  "pdgid/I"  );
  //   pixeltree_->Branch("process", &recHit_.process,"process/I");
  //   pixeltree_->Branch("q",       &recHit_.q       ,"q/F"	   );	   
  //   pixeltree_->Branch("x",       &recHit_.x       ,"x/F"	   );	   
  //   pixeltree_->Branch("y",	     &recHit_.y       ,"y/F"	   );	   
  //   pixeltree_->Branch("xx",	     &recHit_.xx      ,"xx/F"	   );	   
  //   pixeltree_->Branch("xy",	     &recHit_.xy      ,"xy/F"	   );	   
  //   pixeltree_->Branch("yy",	     &recHit_.yy      ,"yy/F"	   );	   
  //   pixeltree_->Branch("row",     &recHit_.row     ,"row/F"     );		   
  //   pixeltree_->Branch("col",     &recHit_.col     ,"col/F"     );		   
  //   pixeltree_->Branch("gx",	     &recHit_.gx      ,"gx/F"	   );	   
  //   pixeltree_->Branch("gy",	     &recHit_.gy      ,"gy/F"	   );	   
  //   pixeltree_->Branch("gz",	     &recHit_.gz      ,"gz/F"	   );	   
  //   pixeltree_->Branch("subid",   &recHit_.subid   ,"subid/I"   );	   
  //   pixeltree_->Branch("module",  &recHit_.module  ,"module/I"  );	   
  //   pixeltree_->Branch("layer",   &recHit_.layer    ,"layer/I"  );	   
  //   pixeltree_->Branch("ladder",  &recHit_.ladder   ,"ladder/I" );	   
  //   pixeltree_->Branch("disk",    &recHit_.disk     ,"disk/I"   );	   
  //   pixeltree_->Branch("blade",   &recHit_.blade    ,"blade/I"  );	   
  //   pixeltree_->Branch("panel",   &recHit_.panel    ,"panel/I"  );	   
  //   pixeltree_->Branch("side",    &recHit_.side     ,"side/I"   );	   
  //   pixeltree_->Branch("nsimhit", &recHit_.nsimhit  ,"nsimhit/I");	   
  //   pixeltree_->Branch("spreadx", &recHit_.spreadx  ,"spreadx/I");	   
  //   pixeltree_->Branch("spready", &recHit_.spready  ,"spready/I");	   
  //   pixeltree_->Branch("hx",	     &recHit_.hx       ,"hx/F");	   
  //   pixeltree_->Branch("hy",	     &recHit_.hy       ,"hy/F");	   
  //   pixeltree_->Branch("tx",      &recHit_.tx       ,"tx/F");	   
  //   pixeltree_->Branch("ty",	     &recHit_.ty       ,"ty/F");	   
  //   pixeltree_->Branch("tz",	     &recHit_.tz       ,"tz/F");	   
  //   pixeltree_->Branch("theta",   &recHit_.theta    ,"theta/F" );	   
  //   pixeltree_->Branch("phi",     &recHit_.phi      ,"phi/F"    );		   
  //   pixeltree_->Branch("DgN",     &recHit_.fDgN     ,"DgN/I"    );		   
  //   pixeltree_->Branch("DgRow",   recHit_.fDgRow    ,"DgRow[DgN]/I"  );	   
  //   pixeltree_->Branch("DgCol",   recHit_.fDgCol    ,"DgCol[DgN]/I"   ); 
  //   pixeltree_->Branch("DgDetId", recHit_.fDgDetId  ,"DgDetId[DgN]/I" );  
  //   pixeltree_->Branch("DgAdc",   recHit_.fDgAdc    ,"DgAdc[DgN]/F"   ); 
  //   pixeltree_->Branch("DgCharge",recHit_.fDgCharge ,"DgCharge[DgN]/F");  
  //   pixeltree_->Branch("DgClI",   recHit_.fDgClI    ,"DgClI[DgN]/I"   ); 
  //   pixeltree_->Branch("ClN",     &recHit_.fClN     ,"ClN/I"	  );	   
  //   pixeltree_->Branch("ClDgN",   recHit_.fClDgN    ,"ClDgN[ClN]/I"  ); 
  //   pixeltree_->Branch("ClDgI",   recHit_.fClDgI    ,"ClDgI[ClN][100]/I");

}
 
// Functions that gets called by framework every event
void PixelStdHitNtuplizer::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  //Retrieve tracker topology from geometry
  edm::ESHandle<TrackerTopology> tTopoHandle;
  es.get<IdealGeometryRecord>().get(tTopoHandle);
  const TrackerTopology* const tTopo = tTopoHandle.product();
 
  // geometry setup
  edm::ESHandle<TrackerGeometry>        geometry;
 
  es.get<TrackerDigiGeometryRecord>().get(geometry);
  const TrackerGeometry*  theGeometry = &(*geometry);
 
  // fastsim rechits
  //edm::Handle<SiTrackerGSRecHit2DCollection> theGSRecHits;
  //edm::InputTag hitProducer;
  //hitProducer = conf_.getParameter<edm::InputTag>("HitProducer");
  //e.getByLabel(hitProducer, theGSRecHits);
 
  std::vector<PSimHit> matched;
  std::vector<PSimHit>::const_iterator closest_simhit;

  edm::Handle<SiPixelRecHitCollection> recHitColl;
  e.getByLabel( src_, recHitColl);
 
  // for finding matched simhit
  TrackerHitAssociator associate( e, conf_ );
 
  //  std::cout << " Step A: Standard RecHits found " << (recHitColl.product())->dataSize() << std::endl;
  if((recHitColl.product())->dataSize() > 0) {
    SiPixelRecHitCollection::const_iterator recHitIdIterator      = (recHitColl.product())->begin();
    SiPixelRecHitCollection::const_iterator recHitIdIteratorEnd   = (recHitColl.product())->end();
 
    std::string detname ;

    // Loop over Detector IDs
    for ( ; recHitIdIterator != recHitIdIteratorEnd; recHitIdIterator++) {
      SiPixelRecHitCollection::DetSet detset = *recHitIdIterator;
 
      if( detset.empty() ) continue;
      DetId detId = DetId(detset.detId()); // Get the Detid object
 
      const GeomDet* geomDet( theGeometry->idToDet(detId) );
       
      // Loop over rechits for this detid
      SiPixelRecHitCollection::DetSet::const_iterator rechitRangeIteratorBegin = detset.begin();
      SiPixelRecHitCollection::DetSet::const_iterator rechitRangeIteratorEnd   = detset.end();
      SiPixelRecHitCollection::DetSet::const_iterator iterRecHit;
      for ( iterRecHit = rechitRangeIteratorBegin; 
            iterRecHit != rechitRangeIteratorEnd; ++iterRecHit) {
	// get matched simhit
	matched.clear();
	matched = associate.associateHit(*iterRecHit);
	if ( !matched.empty() ) {
	  float closest = 9999.9;
	  std::vector<PSimHit>::const_iterator closestit = matched.begin();
	  LocalPoint lp = iterRecHit->localPosition();
	  float rechit_x = lp.x();
	  float rechit_y = lp.y();
	  //loop over simhits and find closest
	  for (std::vector<PSimHit>::const_iterator m = matched.begin(); m<matched.end(); m++) 
	    {
	      float sim_x1 = (*m).entryPoint().x();
	      float sim_x2 = (*m).exitPoint().x();
	      float sim_xpos = 0.5*(sim_x1+sim_x2);
	      float sim_y1 = (*m).entryPoint().y();
	      float sim_y2 = (*m).exitPoint().y();
	      float sim_ypos = 0.5*(sim_y1+sim_y2);
             
	      float x_res = sim_xpos - rechit_x;
	      float y_res = sim_ypos - rechit_y;
	      float dist = sqrt(x_res*x_res + y_res*y_res);
	      if ( dist < closest ) {
		closest = dist;
		closestit = m;
	      }
	    } // end of simhit loop
	  closest_simhit = closestit;
	} // end matched emtpy
	unsigned int subid = detId.subdetId();
	int detid_db = detId.rawId();
	int layer_num = -99,ladder_num=-99,module_num=-99,disk_num=-99,blade_num=-99,panel_num=-99,side_num=-99;       
	if ( ( subid == PixelSubdetector::PixelBarrel ) || ( subid == PixelSubdetector::PixelEndcap ) ) {
	  // 1 = PXB, 2 = PXF
	  if ( subid ==  PixelSubdetector::PixelBarrel ) {
	    layer_num   = tTopo->pxbLayer(detId.rawId());
	    ladder_num  = tTopo->pxbLadder(detId.rawId());
	    module_num  = tTopo->pxbModule(detId.rawId());
	    //	    std::cout <<"\ndetId = "<<subid<<" : "<<tTopo->pxbLayer(detId.rawId())<<" , "<<tTopo->pxbLadder(detId.rawId())<<" , "<< tTopo->pxbModule(detId.rawId());
	  } else if ( subid ==  PixelSubdetector::PixelEndcap ) {
	    module_num  = tTopo->pxfModule(detId());
	    disk_num    = tTopo->pxfDisk(detId());
	    blade_num   = tTopo->pxfBlade(detId());
	    panel_num   = tTopo->pxfPanel(detId());
	    side_num    = tTopo->pxfSide(detId());
	  }
	  int num_simhit = matched.size();
	  fillPRecHit(  detid_db, subid,layer_num,ladder_num,module_num,disk_num,blade_num,panel_num,side_num,
			iterRecHit, num_simhit, closest_simhit, geomDet );
	  fillEvt(e);
	  pixeltree_->Fill();
	  init();
	}
      } // end of rechit loop
    } // end of detid loop
  } // end of loop test on recHitColl size
 
  // Now loop over recotracks
 
  edm::Handle<View<reco::Track> >  trackCollection;
  edm::InputTag trackProducer;
  trackProducer = conf_.getParameter<edm::InputTag>("trackProducer");
  e.getByLabel(trackProducer, trackCollection);
 
  /*
    std::cout << " num of reco::Tracks with "
    << trackProducer.process()<<":"
    << trackProducer.label()<<":"
    << trackProducer.instance()
    << ": " << trackCollection->size() << "\n";
  */
  int rT = 0;
  for(View<reco::Track>::size_type i=0; i<trackCollection->size(); ++i){
    ++rT;
    RefToBase<reco::Track> track(trackCollection, i);
    int iT = 0;
    //std::cout << " num of hits for track " << rT << " = " << track->recHitsSize() << std::endl;
    for(trackingRecHit_iterator ih=track->recHitsBegin(); ih != track->recHitsEnd(); ++ih) {
      ++iT;
      //std::cout<<" analyzing hit n. "<<iT<<std::endl;
      TrackingRecHit * hit = (*ih)->clone();
      const DetId& detId =  hit->geographicalId();
      const GeomDet* geomDet( theGeometry->idToDet(detId) );

      const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>(hit);
         
      if(pixhit){
	if(pixhit->isValid() ) {
	
	  // get matched simhit
	  matched.clear();
	  matched = associate.associateHit(*pixhit);
	  	  
	  if ( !matched.empty() ) {
	    float closest = 9999.9;
	    std::vector<PSimHit>::const_iterator closestit = matched.begin();
	    LocalPoint lp = pixhit->localPosition();
	    float rechit_x = lp.x();
	    float rechit_y = lp.y();

	    //loop over simhits and find closest	   	    
	    for (std::vector<PSimHit>::const_iterator m = matched.begin(); m<matched.end(); m++) 
	      {
		float sim_x1 = (*m).entryPoint().x();
		float sim_x2 = (*m).exitPoint().x();
		float sim_xpos = 0.5*(sim_x1+sim_x2);
		float sim_y1 = (*m).entryPoint().y();
		float sim_y2 = (*m).exitPoint().y();
		float sim_ypos = 0.5*(sim_y1+sim_y2);
		
		float x_res = sim_xpos - rechit_x;
		float y_res = sim_ypos - rechit_y;
		float dist = sqrt(x_res*x_res + y_res*y_res);
		if ( dist < closest ) {
		  closest = dist;
		  closestit = m;
		}
	      } // end of simhit loop
	    closest_simhit = closestit;
	  } // end matched emtpy
	  
	  int num_simhit = matched.size();
	  
	  int layer_num = -99,ladder_num=-99,module_num=-99,disk_num=-99,blade_num=-99,panel_num=-99,side_num=-99;
	  
	  unsigned int subid = detId.subdetId();
	  int detid_db = detId.rawId();
	  if ( ( subid == PixelSubdetector::PixelBarrel ) || ( subid == PixelSubdetector::PixelEndcap ) ) {
	    // 1 = PXB, 2 = PXF
	    if ( subid ==  PixelSubdetector::PixelBarrel ) {
	      layer_num   = tTopo->pxbLayer(detId.rawId());
	      ladder_num  = tTopo->pxbLadder(detId.rawId());
	      module_num  = tTopo->pxbModule(detId.rawId());
	      //std::cout <<"\ndetId = "<<subid<<" : "<<tTopo->pxbLayer(detId.rawId())<<" , "<<tTopo->pxbLadder(detId.rawId())<<" , "<< tTopo->pxbModule(detId.rawId()) <<std::endl;
	    } else if ( subid ==  PixelSubdetector::PixelEndcap ) {
	      module_num  = tTopo->pxfModule(detId());
	      disk_num    = tTopo->pxfDisk(detId());
	      blade_num   = tTopo->pxfBlade(detId());
	      panel_num   = tTopo->pxfPanel(detId());
	      side_num    = tTopo->pxfSide(detId());
	    }
	    
	    fillPRecHit(detid_db, subid, layer_num,ladder_num,module_num,disk_num,blade_num,panel_num,side_num, 
			ih, num_simhit, closest_simhit, geomDet );
	    	    
	    fillEvt(e);	    
	    //	    pixeltree2_->Fill();
	    init();
	    
	    /*
	      TrackingRecHit * hit = (*ih)->clone();
	      LocalPoint lp = hit->localPosition();
	      LocalError le = hit->localPositionError();
	      //            std::cout << "   lp x,y = " << lp.x() << " " << lp.y() << " lpe xx,xy,yy = "
	      //                  << le.xx() << " " << le.xy() << " " << le.yy() << std::endl;
	      std::cout << "Found RecHit in " << detname << " from detid " << detId.rawId()
	      << " subdet = " << subdetId
	      << " layer = " << layerNumber
	      << "global x/y/z/r = "
	      << geomDet->surface().toGlobal(lp).x() << " " 
	      << geomDet->surface().toGlobal(lp).y() << " " 
	      << geomDet->surface().toGlobal(lp).z() << " " 
	      << geomDet->surface().toGlobal(lp).perp() 
	      << " err x/y = " << sqrt(le.xx()) << " " << sqrt(le.yy()) << std::endl;
	    */
	  } // if ( (subid==1)||(subid==2) ) 
	} // if SiPixelHit is valid
      } // if cast is possible to SiPixelHit
      delete pixhit;
    } //end of loop on tracking rechits
  } // end of loop on recotracks
             
} // end analyze function
 
// Function for filling in all the rechits
// I know it is lazy to pass everything, but I'm doing it anyway. -EB
void PixelStdHitNtuplizer::fillPRecHit(const int detid_db, const int subid, 
				  const int layer_num,const int ladder_num,const int module_num,
				  const int disk_num,const int blade_num,const int panel_num,const int side_num,
				  SiPixelRecHitCollection::DetSet::const_iterator pixeliter,
				  const int num_simhit,
				  std::vector<PSimHit>::const_iterator closest_simhit,
				  const GeomDet* PixGeom)
{
  LocalPoint lp = pixeliter->localPosition();
  LocalError le = pixeliter->localPositionError();
 
  recHit_.x = lp.x();
  recHit_.y = lp.y();
  recHit_.xx = le.xx();
  recHit_.xy = le.xy();
  recHit_.yy = le.yy();
  //MeasurementPoint mp = topol->measurementPosition(LocalPoint(recHit_.x, recHit_.y));
  //recHit_.row = mp.x();
  //recHit_.col = mp.y();
  GlobalPoint GP = PixGeom->surface().toGlobal(pixeliter->localPosition());
  recHit_.gx = GP.x();
  recHit_.gy = GP.y();
  recHit_.gz = GP.z();
  GlobalPoint GP0 = PixGeom->surface().toGlobal(LocalPoint(0,0,0));
  recHit_.theta = GP0.theta(); 
  recHit_.phi = GP0.phi(); 

  SiPixelRecHit::ClusterRef const& Cluster =  pixeliter->cluster();
  recHit_.q = Cluster->charge();
  recHit_.spreadx = Cluster->sizeX();
  recHit_.spready = Cluster->sizeY();
 
  recHit_.subid = subid;
  recHit_.nsimhit = num_simhit;
  
  recHit_.layer = layer_num;
  recHit_.ladder= ladder_num;
  recHit_.module= module_num;
  recHit_.module= module_num;
  recHit_.disk  = disk_num;
  recHit_.blade = blade_num;
  recHit_.panel = panel_num;
  recHit_.side  = side_num;

  /*-- --*/
  const PixelGeomDetUnit *theGeomDet = dynamic_cast<const PixelGeomDetUnit*> (PixGeom );
  const PixelTopology * topol = &(theGeomDet->specificTopology());
  

  if ( subid == 1 ) {
    std::cout << "BPIX Layer "<< layer_num << " RectangularPixelTopology " << "nrows " << topol->nrows() << ", ncols " << topol->ncolumns() << ", pitchx " << topol->pitch().first << ", pitchy " << topol->pitch().second << std::endl;
  }


  if ( Cluster.isNonnull() ) { 
            // -- Get digis of this cluster
      const std::vector<SiPixelCluster::Pixel>& pixvector = Cluster->pixels();
      //      std::cout << "  Found " << pixvector.size() << " pixels for this cluster " << std::endl;
      for (unsigned int i = 0; i < pixvector.size(); ++i) {
	if (recHit_.fDgN > DIGIMAX - 1) break;
	SiPixelCluster::Pixel holdpix = pixvector[i];
	
	recHit_.fDgRow[recHit_.fDgN]    = holdpix.x;
	recHit_.fDgCol[recHit_.fDgN]    = holdpix.y;
	//	std::cout << "holdpix " << holdpix.x << " " <<  holdpix.y << std::endl;
	recHit_.fDgDetId[recHit_.fDgN]  = detid_db;
	recHit_.fDgAdc[recHit_.fDgN]    = -99.;
	recHit_.fDgCharge[recHit_.fDgN] = holdpix.adc/1000.;
	
	recHit_.fDgClI[recHit_.fDgN] = recHit_.fClN;
	
	// -- fill pointer to this digi in cluster digi index array
	if ((signed)i < DGPERCLMAX) {
	  	  recHit_.fClDgI[recHit_.fClN][i] = recHit_.fDgN;
	  	  recHit_.fClDgN[recHit_.fClN] += 1;
	} else {
	  	  recHit_.fClDgI[recHit_.fClN][DGPERCLMAX-1] = -98;
	  	  recHit_.fClDgN[recHit_.fClN] = 99;
	}

	++recHit_.fDgN;
	
      }
      ++recHit_.fClN;
      
    } // if ( Cluster.isNonnull() )
  else 
    {
      std::cout << "Pixel rechits with no associated cluster ?!?!?!!!!!!!!!!!!!!!!!!!!!!! " << std::endl;
	// for (std::vector<TrajectoryMeasurement>::const_iterator tmeasIt = tmeasColl.begin(); ...
    }
	/*-- --*/

  //std::cout << "num_simhit = " << num_simhit << std::endl;
  if(num_simhit > 0) {

    recHit_.pdgid = (*closest_simhit).particleType();
    recHit_.process= (*closest_simhit).processType();

    float sim_x1 = (*closest_simhit).entryPoint().x();
    float sim_x2 = (*closest_simhit).exitPoint().x();
    recHit_.hx = 0.5*(sim_x1+sim_x2);
    float sim_y1 = (*closest_simhit).entryPoint().y();
    float sim_y2 = (*closest_simhit).exitPoint().y();
    recHit_.hy = 0.5*(sim_y1+sim_y2);

    recHit_.tx = (*closest_simhit).localDirection().x();
    recHit_.ty = (*closest_simhit).localDirection().y();
    recHit_.tz = (*closest_simhit).localDirection().z();
   // alpha: angle with respect to local x axis in local (x,z) plane
   // float cotalpha = sim_xdir/sim_zdir;
   // beta: angle with respect to local y axis in local (y,z) plane
   // float cotbeta = sim_ydir/sim_zdir;

    //std::cout << "num_simhit x, y = " << 0.5*(sim_x1+sim_x2) << " " << 0.5*(sim_y1+sim_y2) << std::endl;
  }
  /*
    std::cout << "Found RecHit in " << subid
    << " global x/y/z : "
    << PixGeom->surface().toGlobal(pixeliter->localPosition()).x() << " " 
    << PixGeom->surface().toGlobal(pixeliter->localPosition()).y() << " " 
    << PixGeom->surface().toGlobal(pixeliter->localPosition()).z() << std::endl;
  */
}

// Function for filling in on track rechits
void PixelStdHitNtuplizer::fillPRecHit(const int detid_db, const int subid, 
				     const int layer_num,const int ladder_num,const int module_num,
				     const int disk_num,const int blade_num,const int panel_num,const int side_num,
				     trackingRecHit_iterator ih,
				     const int num_simhit,
				     std::vector<PSimHit>::const_iterator closest_simhit,
				     const GeomDet* PixGeom)
{
  TrackingRecHit * pixeliter = (*ih)->clone(); 
  LocalPoint lp = pixeliter->localPosition();
  LocalError le = pixeliter->localPositionError();
 
  recHit_.x = lp.x();
  recHit_.y = lp.y();
  recHit_.xx = le.xx();
  recHit_.xy = le.xy();
  recHit_.yy = le.yy();
  GlobalPoint GP = PixGeom->surface().toGlobal(pixeliter->localPosition());
  recHit_.gx = GP.x();
  recHit_.gy = GP.y();
  recHit_.gz = GP.z();
  GlobalPoint GP0 = PixGeom->surface().toGlobal(LocalPoint(0,0,0));
  recHit_.theta = GP0.theta(); 
  recHit_.phi = GP0.phi(); 
  recHit_.subid = subid;

  //std::cout<<"before getting the cluster"<<std::endl;

  SiPixelRecHit::ClusterRef const& Cluster =  dynamic_cast<const SiPixelRecHit*>(pixeliter)->cluster();
  recHit_.q = Cluster->charge();
  recHit_.spreadx = Cluster->sizeX();
  recHit_.spready = Cluster->sizeY();

  recHit_.nsimhit = num_simhit;

  recHit_.layer = layer_num;
  recHit_.ladder= ladder_num;
  recHit_.module= module_num;
  recHit_.module= module_num;
  recHit_.disk  = disk_num;
  recHit_.blade = blade_num;
  recHit_.panel = panel_num;
  recHit_.side  = side_num;

  //std::cout << "num_simhit = " << num_simhit << std::endl;
  if(num_simhit > 0) {

    recHit_.pdgid = (*closest_simhit).particleType();
    recHit_.process= (*closest_simhit).processType();

    float sim_x1 = (*closest_simhit).entryPoint().x();
    float sim_x2 = (*closest_simhit).exitPoint().x();
    recHit_.hx = 0.5*(sim_x1+sim_x2);
    float sim_y1 = (*closest_simhit).entryPoint().y();
    float sim_y2 = (*closest_simhit).exitPoint().y();
    recHit_.hy = 0.5*(sim_y1+sim_y2);

    recHit_.tx = (*closest_simhit).localDirection().x();
    recHit_.ty = (*closest_simhit).localDirection().y();
    recHit_.tz = (*closest_simhit).localDirection().z();
   // alpha: angle with respect to local x axis in local (x,z) plane
   // float cotalpha = sim_xdir/sim_zdir;
   // beta: angle with respect to local y axis in local (y,z) plane
   // float cotbeta = sim_ydir/sim_zdir;

    //std::cout << "num_simhit x, y = " << 0.5*(sim_x1+sim_x2) << " " << 0.5*(sim_y1+sim_y2) << std::endl;
  }

  delete pixeliter;

}
 
void
PixelStdHitNtuplizer::fillEvt(const edm::Event& E)
{
  evt_.run = E.id().run();
  evt_.evtnum = E.id().event();
}
 
void PixelStdHitNtuplizer::init()
{
  evt_.init();
  recHit_.init();
}
 
//define this as a plug-in
DEFINE_FWK_MODULE(PixelStdHitNtuplizer);
 
