#include "TH1.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


class LHEAnalyzer : public edm::EDAnalyzer {

public:
  /// default constructor
  explicit LHEAnalyzer(const edm::ParameterSet&);
  /// default destructor
  virtual ~LHEAnalyzer();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  TH1F* h_nb_;
  TH1F* h_statusb_;
  TH1F* h_massb_;
  TH1F* h_ptb1_;
  TH1F* h_etab1_;
  TH1F* h_ptball_;
  TH1F* h_etaball_;

};


LHEAnalyzer::LHEAnalyzer(const edm::ParameterSet& iConfig)
{
} 

LHEAnalyzer::~LHEAnalyzer()
{
} 



void LHEAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// based on http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/Validation/EventGenerator/plugins/DuplicationChecker.cc
// http://arxiv.org/pdf/hep-ph/0609017v1
{
  edm::Handle<LHEEventProduct> evt;
  iEvent.getByType( evt );

  const lhef::HEPEUP hepeup_ = evt->hepeup();
  const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP;
  const std::vector<int> idup_ = hepeup_.IDUP;
  const std::vector<int> istup_ = hepeup_.ISTUP;
  
  // number of b in the event
  int nb = 0; 
  // initialize pT and eta of the leading b
  double ptb_max = -1.;
  double etab_max;

  // loop on the particles in the event
  int nup = hepeup_.NUP;
  for ( int iup=0; iup<nup; iup++) {
    if ( idup_[iup]==5 || idup_[iup]==-5 ) {
      nb++;
      h_statusb_->Fill(istup_[iup]);

	if (istup_[iup]==1 ) {
      // dump the 4-momentum into a LorentzVector
      math::XYZTLorentzVectorD p4((pup_[iup])[0],(pup_[iup])[1],(pup_[iup])[2],(pup_[iup])[3]);      
      double mb=p4.mass();
      double ptb = p4.pt();
      double etab = p4.eta();

      if ( ptb > ptb_max ) {
	ptb_max = ptb; 
	etab_max = etab;
      }
      h_massb_->Fill(mb);
      h_ptball_->Fill(ptb);
      h_etaball_->Fill(etab);
    }
}
  }

  if ( ptb_max>0 ) {
    h_ptb1_->Fill(ptb_max);
    h_etab1_->Fill(etab_max);
  }
  h_nb_->Fill(nb);

}


void LHEAnalyzer::beginJob()
{
  // register to the TFileService
  edm::Service<TFileService> fs;
  
  // book histograms:
  h_nb_	     = fs->make<TH1F>("HistNb","number of b",7,-0.5,6.5);
  h_massb_   = fs->make<TH1F>("HistMassb","mass of b",120,-0.1,5.9);
  h_statusb_ = fs->make<TH1F>("HistStatusb","status of b",13,-9.5,3.5);
  h_ptb1_    = fs->make<TH1F>("HistPtLeadingb","pT leading b",150,0.,150.);
  h_etab1_   = fs->make<TH1F>("HistEtaLeadingb","eta leading b",25,-5.,5.);
  h_ptball_  = fs->make<TH1F>("HistPtAllb","pT all b",150,0.,150.);
  h_etaball_ = fs->make<TH1F>("HistEtaAllb","eta all b",25,-5.,5.);

}

void LHEAnalyzer::endJob(){
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LHEAnalyzer);
