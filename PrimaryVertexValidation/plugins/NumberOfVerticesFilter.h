#ifndef NumberOfVerticesFilter_h
#define NumberOfVerticesFilter_h

# include "FWCore/Framework/interface/EDFilter.h"

class NumberOfVerticesFilter: public edm::EDFilter {
  
public:
  // Constructor
  explicit NumberOfVerticesFilter (const edm::ParameterSet &pset);
  
  /// Destructor
  virtual ~NumberOfVerticesFilter (void);
  
  // Operations
  virtual bool filter (edm::Event &iEvent, const edm::EventSetup &eventSetup);
  
private:
  
  // from cfg
  std::string  m_theVerticesLabel;
  int m_MAXnumberOfVertices;
  int m_MINnumberOfVertices; 
  bool m_verbose;
  // Counters
  int m_numberOfEvents;
  int nOfflineVertices_;
  bool isGoodEvent;

};

#endif
 
