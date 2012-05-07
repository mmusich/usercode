#include "ZbbAnalysis/Tools/plugins/NumberOfVerticesFilter.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
 
using namespace std;
using namespace edm;

NumberOfVerticesFilter::NumberOfVerticesFilter (const ParameterSet &pset)
{
  m_theVerticesLabel = pset.getUntrackedParameter<string>("VertexCollectionLabel");
  m_MINnumberOfVertices = pset.getUntrackedParameter<int>("MINnumberOfVertices");
  m_MAXnumberOfVertices = pset.getUntrackedParameter<int>("MAXnumberOfVertices");
  m_verbose = pset.getUntrackedParameter<bool>("verbose");
  m_numberOfEvents = 0;
}

NumberOfVerticesFilter::~NumberOfVerticesFilter (void) 
{
}

bool 
NumberOfVerticesFilter::filter (Event &iEvent, const EventSetup &eventSetup) 
{

  if(m_verbose) cout<<"inside filter"<<endl;
  nOfflineVertices_=0;
  
  //=======================================================
  // Retrieve offline vartex information
  //=======================================================
  isGoodEvent=false;

  edm::Handle<reco::VertexCollection> vertices;
  try {
    iEvent.getByLabel(m_theVerticesLabel, vertices);
  } catch (...) {
      cout << "No offlinePrimaryVertices found!" << endl;
  }
 
  nOfflineVertices_ = vertices.product()->size();
  
  if (nOfflineVertices_>= m_MINnumberOfVertices && nOfflineVertices_<= m_MAXnumberOfVertices) 
    {
      if(m_verbose){  cout << "Event no. " << m_numberOfEvents << " has " 
			   << nOfflineVertices_ << " vertices! " << endl;
      }
      
      isGoodEvent = true;
    }
  
  return isGoodEvent;
}

//define this as a plug-in
DEFINE_FWK_MODULE(NumberOfVerticesFilter);
