#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include <TH1.h>
#include <TH2.h>
#include <Math/VectorUtil.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"

#include "ZbbAnalysis/AnalysisStep/interface/ZbbUtils.h"

class PatBJetTagAnalyzer : public edm::EDAnalyzer  {
    public: 
	/// constructor and destructor
	PatBJetTagAnalyzer(const edm::ParameterSet &params);
	~PatBJetTagAnalyzer();

	// virtual methods called from base class EDAnalyzer
	virtual void beginJob();
	virtual void analyze(const edm::Event &event, const edm::EventSetup &es);
        double btagSSVHE( const reco::SecondaryVertexTagInfo* svTagInfo );
        double btagSSVHP( const reco::SecondaryVertexTagInfo* svTagInfo );
        double btagNegativeSSVHE( const reco::SecondaryVertexTagInfo* svTagInfo );
        double btagNegativeSSVHP( const reco::SecondaryVertexTagInfo* svTagInfo );

    private:
	// configuration input parameters
	edm::InputTag jets_;
	edm::InputTag tracks_;
	edm::InputTag beamSpot_;
	edm::InputTag primaryVertices_;

       	// configuration cuts
        bool ismc_;
        bool redoBTag_;
        std::string bTagAlgoWP_;
  	double jetPtCut_;		// minimum (uncorrected) jet energy
	double jetEtaCut_;		// maximum |eta| for jet
	double maxDeltaR_;		// angle between jet and tracks
	double minPt_;			// track pt quality cut
	unsigned int minPixelHits_;	// minimum number of pixel hits
	unsigned int minTotalHits_;	// minimum number of total hits
	unsigned int nThTrack_;		// n-th hightest track to choose

        // one group of plots per jet flavour;
        struct Plots {
	  TH1 *discrTCHE, *discrTCHP,  *discrSSV, *discrSSVHP,  *discrSSVHE, *discrCSV, *discrJP, *jetPt;
	  TH1 *discrTCHE_fromTags, *discrTCHP_fromTags, *discrJP_fromTags, *jetPt_fromTags;      
	  TH1 *allIP, *allIPErr, *allIPSig;
	  TH1 *trackIP, *trackIPErr, *trackIPSig;
	  TH1 *negativeIP, *negativeIPErr, *negativeIPSig;
	  TH1 *nTracks, *allDeltaR;
	  TH1 *SVnVertices, *SVJetdeltaR, *SVmass, *SVmass2, *SVdist, *SVdistErr, *SVdistSig, *SVnTracks,*SVchi2;
	  TH1 *SVdist3D, *SVdistErr3D, *SVdistSig3D, *DSSVHE, *DSSVHP;        
	  TH1 *trackIPfromTagInfo, *allIPErrfromTagInfo, *allIPSigfromTagInfo, *trackptfromTagInfo, *trackpixelhitsFromTagInfo; 
	  TH2 *MCjetsVstags;
	  TH1 *MCsecvtxMass_b, *MCsecvtxMass_c, *MCsecvtxMass_q, *MCsecvtxMass_X;
	  TH1 *MCsecvtxMass_bb,*MCsecvtxMass_cc,*MCsecvtxMass_qq,*MCsecvtxMass_bc,*MCsecvtxMass_bq,*MCsecvtxMass_cq,*MCsecvtxMass_XX, *MCsecvtxMass_All; 
	} plots_;
};

PatBJetTagAnalyzer::PatBJetTagAnalyzer(const edm::ParameterSet &params) :
	jets_(params.getParameter<edm::InputTag>("jets")),
	tracks_(params.getParameter<edm::InputTag>("tracks")),
	beamSpot_(params.getParameter<edm::InputTag>("beamSpot")),
	primaryVertices_(params.getParameter<edm::InputTag>("primaryVertices")),	
	ismc_(params.getParameter<bool>("isMC")),
	redoBTag_(params.getParameter<bool>("redoBTag")),
	bTagAlgoWP_(params.getParameter<std::string>("bTagAlgoWP")),
	jetPtCut_(params.getParameter<double>("jetPtCut")),
	jetEtaCut_(params.getParameter<double>("jetEtaCut")),
	maxDeltaR_(params.getParameter<double>("maxDeltaR")),
	minPt_(params.getParameter<double>("minPt")),
	minPixelHits_(params.getParameter<unsigned int>("minPixelHits")),
	minTotalHits_(params.getParameter<unsigned int>("minTotalHits")),
	nThTrack_(params.getParameter<unsigned int>("nThTrack"))
{
}

PatBJetTagAnalyzer::~PatBJetTagAnalyzer()
{
}

void PatBJetTagAnalyzer::beginJob()
{
  // retrieve handle to auxiliary service
  //  used for storing histograms into ROOT file
  edm::Service<TFileService> fs;
  
  // book histograms for all jet flavours
  
  plots_.discrTCHE     = fs->make<TH1F>("discrTCHE","track counting (\"high efficiency\"); TCHE Discriminator; events",150,-100,50);
  plots_.discrTCHP     = fs->make<TH1F>("discrTCHP","track counting (\"high purity\"); TCHP Discriminator; events",150,-100,50);
  plots_.discrSSV      = fs->make<TH1F>("discrSSV","simple secondary vertex; SSV Discriminator; events",150,-5, 10);
  plots_.discrSSVHP    = fs->make<TH1F>("discrSSVHP","simple secondary vertex (HP) discriminator; SSV (HP) Discriminator; events",150,-5, 10);
  plots_.discrSSVHE    = fs->make<TH1F>("discrSSVHE","simple secondary vertex (HE) discriminator; SSV (HE) Discriminator; events",150,-5, 10);
  plots_.discrCSV      = fs->make<TH1F>("discrCSV","combined secondary vertex; CSV Discriminator; events",100,0,1);
  plots_.discrJP       = fs->make<TH1F>("discrJP","Jet probability discriminator;JP Discriminator; events",100,0.,3.);
  plots_.jetPt         = fs->make<TH1F>("jetPt","Jet transverse momentum; jet p_{T} (GeV); events",100,0.,200.);
  
  plots_.discrTCHE_fromTags = fs->make<TH1F>("discrTCHE_fromTags","track counting (\"high efficiency\"); TCHE Discriminator; events",150,-100,50);
  plots_.discrTCHP_fromTags = fs->make<TH1F>("discrTCHP_fromTags","track counting (\"high purity\"); TCHP Discriminator; events",150,-100,50);
  plots_.discrJP_fromTags   = fs->make<TH1F>("discrJP_fromTags","Jet probability discriminator;JP Discriminator; events",100,0.,3.);
  plots_.jetPt_fromTags     = fs->make<TH1F>("jetPt_fromTags","Jet transverse momentum;jet p_{T} (GeV); events",100,0.,200.);

  plots_.allIP         = fs->make<TH1F>("allIP","signed IP for all tracks; all tracks IP (cm); tracks",100, -0.2,0.2);
  plots_.allIPErr      = fs->make<TH1F>("allIPErr","error of signed IP for all tracks; all tracks #sigma_{IP} (cm); tracks",100,0,0.05);
  plots_.allIPSig      = fs->make<TH1F>("allIPSig","signed IP significance for all tracks; all tracks IP/#sigma_{IP}; tracks",200,-20,20);
  plots_.trackIP       = fs->make<TH1F>("trackIP","signed IP for selected positive track; IP (positive) (cm); tracks",200,-0.2,0.2);
  plots_.trackIPErr    = fs->make<TH1F>("trackIPErr","error of signed IP for selected positive track; #sigma_{IP} (positive) (cm); tracks",100, 0,0.05);
  plots_.trackIPSig    = fs->make<TH1F>("trackIPSig","signed IP significance for selected positive track; IP/#sigma_{IP} (positive); tracks",200,-20,20);
  plots_.negativeIP    = fs->make<TH1F>("negativeIP","signed IP for selected negative track; IP (negative) (cm); tracks",200,-0.2,0.2);
  plots_.negativeIPErr = fs->make<TH1F>("negativeIPErr","error of signed IP for selected negative track; #sigma_{IP} (negative) (cm); tracks",100,0,0.05);
  plots_.negativeIPSig = fs->make<TH1F>("negativeIPSig","signed IP significance for selected negative track; IP/#sigma_{IP} (negative); tracks",200,-20,20);
  plots_.nTracks       = fs->make<TH1F>("nTracks","number of usable tracks; n_{tracks} (number of tracks); events",30,-0.5,29.5);
  plots_.allDeltaR     = fs->make<TH1F>("allDeltaR","#DeltaR between track and jet; #DeltaR(tk,jet) for all tracks; tracks",100,0,10);
  plots_.SVnVertices   = fs->make<TH1F>("SVnVertices","number of secondary vertices; N_{SV} (number of secondary vertices); events",5,-0.5,4.5);
  plots_.SVJetdeltaR   = fs->make<TH1F>("SVdeltaR","#DeltaR between vertex direction and jet direction; #DeltaR(#vec{PV -SV},#vec{p}_{jet}); jets",100,0.,1);
  plots_.SVmass        = fs->make<TH1F>("SVmass","secondary vertex mass; secondary vertex mass M_{SV} (Gev/c^{2}); events",100,0.,10.);
  plots_.SVmass2       = fs->make<TH1F>("SVmass2","secondary vertex mass; secondary vertex mass M_{SV} (Gev/c^{2}); events",100,0.,10.);
  plots_.SVdist        = fs->make<TH1F>("SVdist","transverse distance between PV and SV; d_{xy}(PV,SV) (cm); events",100,-1.,5.);
  plots_.SVdistErr     = fs->make<TH1F>("SVdistErr","transverse distance error between PV and SV; #sigma_{d(PV,SV)} (cm); events",100,0.,1);
  plots_.SVdistSig     = fs->make<TH1F>("SVdistSig","transverse distance significance between PV and SV; |d_{xy}(PV,SV)|/#sigma_{d_{xy}}; events",100,0., 100.);
  plots_.SVdist3D      = fs->make<TH1F>("SVdist3D","3D distance between PV and SV; d_{3D}(PV,SV) (cm); events",100,-1.,5.);
  plots_.SVdistErr3D   = fs->make<TH1F>("SVdistErr3D","3D distance error between PV and SV; #sigma_{d(PV,SV)} (cm); events",100,0.,1);
  plots_.SVdistSig3D   = fs->make<TH1F>("SVdistSig3D","3D distance significance between PV and SV; |d_{3D}(PV,SV)|/#sigma_{d_{xy}}; events",100,0., 100.);
  plots_.DSSVHE        = fs->make<TH1F>("DSSVHE","signed discriminator D_{SSVHE}; D_{SSVHE} = sign(S)*log(1 + abs(S)); events",100,-10, 10.);
  plots_.DSSVHP        = fs->make<TH1F>("DSSVHP","signed discriminator D_{SSVHP}; D_{SSVHP} = sign(S)*log(1 + abs(S)); events",100,-10, 10.); 

  plots_.SVnTracks     = fs->make<TH1F>("SVnTracks","number of tracks at secondary vertex; number of track in SV N_{tracks}^{SV}; events",20,0., 20.);
  plots_.SVchi2        = fs->make<TH1F>("SVchi2","secondary vertex fit #chi^{2}; #chi^{2} of SV #chi^{2}(SV)",100,0., 100.);
  plots_.trackIPfromTagInfo        = fs->make<TH1F>("trackIPfromTagInfo","signed IP for all tracks; all tracks IP (cm); tracks",100, -0.2,0.2);
  plots_.allIPErrfromTagInfo       = fs->make<TH1F>("allIPErrfromTagInfo","error of signed IP for all tracks; all tracks #sigma_{IP} (cm); tracks",100,0,0.05);
  plots_.allIPSigfromTagInfo       = fs->make<TH1F>("allIPSigfromTagInfo","signed IP significance for all tracks; all tracks IP/#sigma_{IP}; tracks",200,-20,20);
  plots_.trackptfromTagInfo        = fs->make<TH1F>("trackptfromTagInfo","pt of all tracks; p_{T} of tracks in jet (GeV); tracks",100,0.,100.);
  plots_.trackpixelhitsFromTagInfo = fs->make<TH1F>("trackpixelhitsFromTagInfo","pixel hist of all tracks; N^{PXL}_{hits}; tracks",10,-0.5,9.5);

  // MC only variables
  
  plots_.MCjetsVstags = fs->make<TH2F>("MCjetsVstags","Nb. of tags vs nb of jets; N_{b-tags}; N_{jets}",10,-0.5,9.5,10,-0.5,9.5);

  plots_.MCsecvtxMass_b  = fs->make<TH1F>("MCSVmass_b","secondary vertex mass (Z+b); secondary vertex mass M_{SV} (Gev/c^{2}); events",100,0.,10.);	
  plots_.MCsecvtxMass_c  = fs->make<TH1F>("MCSVmass_c","secondary vertex mass (Z+c); secondary vertex mass M_{SV} (Gev/c^{2}); events",100,0.,10.);	
  plots_.MCsecvtxMass_q  = fs->make<TH1F>("MCSVmass_q","secondary vertex mass (Z+q); secondary vertex mass M_{SV} (Gev/c^{2}); events",100,0.,10.);	
  plots_.MCsecvtxMass_X  = fs->make<TH1F>("MCSVmass_X","secondary vertex mass (Z+nonId); secondary vertex mass M_{SV} (Gev/c^{2}); events",100,0.,10.);
					     										
  plots_.MCsecvtxMass_bb = fs->make<TH1F>("MCSVmass_bb","secondary vertex mass (Z+bb); secondary vertex mass M_{SV} (Gev/c^{2}); events",100,0.,10.);
  plots_.MCsecvtxMass_cc = fs->make<TH1F>("MCSVmass_cc","secondary vertex mass (Z+cc); secondary vertex mass M_{SV} (Gev/c^{2}); events",100,0.,10.);
  plots_.MCsecvtxMass_qq = fs->make<TH1F>("MCSVmass_qq","secondary vertex mass (Z+qq); secondary vertex mass M_{SV} (Gev/c^{2}); events",100,0.,10.);
  plots_.MCsecvtxMass_bc = fs->make<TH1F>("MCSVmass_bc","secondary vertex mass (Z+bc); secondary vertex mass M_{SV} (Gev/c^{2}); events",100,0.,10.);
  plots_.MCsecvtxMass_bq = fs->make<TH1F>("MCSVmass_bq","secondary vertex mass (Z+bq); secondary vertex mass M_{SV} (Gev/c^{2}); events",100,0.,10.);
  plots_.MCsecvtxMass_cq = fs->make<TH1F>("MCSVmass_cq","secondary vertex mass (Z+cq); secondary vertex mass M_{SV} (Gev/c^{2}); events",100,0.,10.);
  plots_.MCsecvtxMass_XX = fs->make<TH1F>("MCSVmass_XX","secondary vertex mass (Z+XX); secondary vertex mass M_{SV} (Gev/c^{2}); events",100,0.,10.);
  
  plots_.MCsecvtxMass_All = fs->make<TH1F>("MCSVmass_All","secondary vertex mass (Z+X); secondary vertex mass M_{SV} (Gev/c^{2}); events",100,0.,10.);

}

// helper function to sort the tracks by impact parameter significance
static bool significanceHigher(const Measurement1D &meas1,const Measurement1D &meas2){ return meas1.significance() > meas2.significance(); }

void PatBJetTagAnalyzer::analyze(const edm::Event &event, const edm::EventSetup &es)
{  
  // handle to the jets collection
  edm::Handle<pat::JetCollection> jetsHandle;
  event.getByLabel(jets_, jetsHandle);

  // handle to the primary vertex collection
  edm::Handle<reco::VertexCollection> pvHandle;
  event.getByLabel(primaryVertices_, pvHandle);
  
  // handle to the tracks collection
  edm::Handle<reco::TrackCollection> tracksHandle;
  event.getByLabel(tracks_, tracksHandle);
  
  // handle to the beam spot
  edm::Handle<reco::BeamSpot> beamSpot;
  event.getByLabel(beamSpot_, beamSpot);
  
  // rare case of no reconstructed primary vertex
  if (pvHandle->empty())
    return;

  // extract the position of the (most probable) reconstructed vertex
  math::XYZPoint pv = (*pvHandle)[0].position();


  if(redoBTag_){

    // Get b tag information (TCHE)
    edm::Handle<reco::JetTagCollection> bTagHandleTCHE;
    event.getByLabel("MyTrackCountingHighEffBJetTags", bTagHandleTCHE);
    const reco::JetTagCollection & bTagsTCHE = *(bTagHandleTCHE.product());
    
    for (unsigned int i = 0; i != bTagsTCHE.size(); ++i) {;
      if (bTagsTCHE[i].first->pt() > jetPtCut_ && std::abs(bTagsTCHE[i].first->eta()) < jetEtaCut_ 
	  // && ZbbUtils::isJetIdOk((bTagsJP[i].first),"loose")
	  // && ZbbUtils::isBJet((bTagsTCHE[i].first),bTagAlgoWP_) 
	  ){
	plots_.discrTCHE_fromTags->Fill(bTagsTCHE[i].second);    
      }
    }
    
    // Get b tag information (TCHP)
    edm::Handle<reco::JetTagCollection> bTagHandleTCHP;
    event.getByLabel("MyTrackCountingHighPurBJetTags", bTagHandleTCHP);
    const reco::JetTagCollection & bTagsTCHP = *(bTagHandleTCHP.product()); 
    
    for (unsigned int i = 0; i != bTagsTCHP.size(); ++i) {;
      if (bTagsTCHP[i].first->pt() > jetPtCut_ && std::abs(bTagsTCHP[i].first->eta()) < jetEtaCut_ 
	  // && ZbbUtils::isJetIdOk((bTagsJP[i].first),"loose")
	  // && ZbbUtils::isBJet((bTagsTCHP[i].first),bTagAlgoWP_) 
	  ){
	plots_.discrTCHP_fromTags->Fill(bTagsTCHP[i].second);    
      }
    }
    
    // Get b tag information (JP)
    edm::Handle<reco::JetTagCollection> bTagHandleJP;
    event.getByLabel("MyJetProbabilityBJetTags", bTagHandleJP);
    const reco::JetTagCollection & bTagsJP = *(bTagHandleJP.product());
    
    // Loop over jets and study b tag info.
    for (unsigned int i = 0; i != bTagsJP.size(); ++i) {
      if (bTagsJP[i].first->pt() > jetPtCut_ && std::abs(bTagsJP[i].first->eta()) < jetEtaCut_ 
	  // && ZbbUtils::isJetIdOk((bTagsJP[i].first),"loose")
	  // && ZbbUtils::isBJet((bTagsJP[i].first),bTagAlgoWP_) 
	  ){
	
	//   std::cout<<" Jet "<< i 
	// 	       <<" has b tag discriminator (jetProbabilityBJetTags) = "<<bTagsJP[i].second
	// 	       << " and jet Pt = "<<bTagsJP[i].first->pt()<<std::endl;
	
	plots_.discrJP_fromTags->Fill(bTagsJP[i].second); 
	plots_.jetPt_fromTags->Fill(bTagsJP[i].first->pt());   
      }
    }
  }
  
  // now go through all jets
  for(pat::JetCollection::const_iterator jet = jetsHandle->begin();jet != jetsHandle->end(); ++jet) {
    
    // only look at jets that pass the pt and eta cut
    if (jet->pt() < jetPtCut_ ||  std::abs(jet->eta()) > jetEtaCut_
	|| !ZbbUtils::isJetIdOk((*jet),"loose") 
	|| !ZbbUtils::isBJet((*jet),bTagAlgoWP_) 
	)
      continue;

    plots_.jetPt->Fill(jet->pt());   

    double discrTCHE  = jet->bDiscriminator("trackCountingHighEffBJetTags");
    double discrTCHP  = jet->bDiscriminator("trackCountingHighPurBJetTags");
    double discrSSV   = jet->bDiscriminator("simpleSecondaryVertexBJetTags");
    double discrSSVHP = jet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
    double discrSSVHE = jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    double discrCSV   = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
    double discrJP    = jet->bDiscriminator("jetProbabilityBJetTags");

    plots_.discrTCHE->Fill(discrTCHE);
    plots_.discrTCHP->Fill(discrTCHP);
    plots_.discrSSV->Fill(discrSSV);
    plots_.discrSSVHP->Fill(discrSSVHP);
    plots_.discrSSVHE->Fill(discrSSVHE);
    plots_.discrCSV->Fill(discrCSV);
    plots_.discrJP->Fill(discrJP);

    // this vector will contain IP value / error pairs
    std::vector<Measurement1D> ipValErr;
    
    // Note: PAT is also able to store associated tracks
    //       within the jet object, so we don't have to do the
    //       matching ourselves
    // (see ->associatedTracks() method)
    // However, using this we can't play with the DeltaR cone
    // withour rerunning the PAT producer
    
    // now loop through all tracks
    for(reco::TrackCollection::const_iterator track = tracksHandle->begin();track != tracksHandle->end(); ++track) {
      
      // check the quality criteria
      if (track->pt() < minPt_ || track->hitPattern().numberOfValidHits() < (int)minTotalHits_ || track->hitPattern().numberOfValidPixelHits() < (int)minPixelHits_)
	continue;
      
      // check the Delta R between jet axis and track
      // (Delta_R^2 = Delta_Eta^2 + Delta_Phi^2)
      double deltaR = ROOT::Math::VectorUtil::DeltaR(jet->momentum(), track->momentum());
      plots_.allDeltaR->Fill(deltaR);
     
      // only look at tracks in jet cone
      if (deltaR > maxDeltaR_)
	continue;
      
      // What follows here is an approximation!
      //
      // The dxy() method of the tracks does a linear
      // extrapolation from the track reference position
      // given as the closest point to the beam spot
      // with respect to the given vertex.
      // Since we are using primary vertices, this
      // approximation works well
      //
      // In order to get better results, the
      // "TransientTrack" and specialised methods have
      // to be used.
      // Or look at the "impactParameterTagInfos",
      // which contains the precomputed information
      // from the official b-tagging algorithms
      //
      // see ->tagInfoTrackIP() method
      
      double ipError = track->dxyError();
      double ipValue = std::abs(track->dxy(pv));
      
      // in order to compute the sign, we check if
      // the point of closest approach to the vertex
      // is in front or behind the vertex.
      // Again, we a linear approximation
      // 
      // dot product between reference point and jet axis

      math::XYZVector closestPoint = track->referencePoint() - beamSpot->position();
      // only interested in transverse component, z -> 0
      closestPoint.SetZ(0.);
      double sign = closestPoint.Dot(jet->momentum());
      
      if (sign < 0) ipValue = -ipValue;
      ipValErr.push_back(Measurement1D(ipValue, ipError));
    }
    
    // now order all tracks by significance (highest first)
    std::sort(ipValErr.begin(), ipValErr.end(), significanceHigher);
    
    plots_.nTracks->Fill(ipValErr.size());    
    // plot all tracks
    
    for(std::vector<Measurement1D>::const_iterator iter = ipValErr.begin();iter != ipValErr.end(); ++iter) {
      
      plots_.allIP->Fill(iter->value());      
      plots_.allIPErr->Fill(iter->error());
      // significance (is really just value / error)
      plots_.allIPSig->Fill(iter->significance());
    }
    
    // check if we have enough tracks to fulfill the
    // n-th track requirement
    if (ipValErr.size() < nThTrack_) continue;
    
    // pick the n-th highest track
    const Measurement1D *trk = &ipValErr[nThTrack_ - 1];
    
    plots_.trackIP->Fill(trk->value());    
    plots_.trackIPErr->Fill(trk->error());
    plots_.trackIPSig->Fill(trk->significance());

    // here we define a "negative tagger", i.e. we take
    // the track with the n-lowest signed IP
    // (i.e. preferrably select tracks that appear to become
    //  from "behind" the primary vertex, supposedly mismeasured
    //  tracks really coming from the primary vertex, and
    //  the contamination of displaced tracks should be small)
    trk = &ipValErr[ipValErr.size() - nThTrack_];
    
    plots_.negativeIP->Fill(trk->value());
    plots_.negativeIPErr->Fill(trk->error());
    plots_.negativeIPSig->Fill(trk->significance());

    //============================================================================
    //
    // Tag Infos
    //
    //============================================================================

    // For debug purposes, print out the offline tracks in all triggered jets.
    const reco::TrackIPTagInfo* tagInfo = jet->tagInfoTrackIP("");
    if (tagInfo != 0) {
      const edm::RefVector<reco::TrackCollection>& tracks = tagInfo->selectedTracks();
      const std::vector<reco::TrackIPTagInfo::TrackIPData>& ip = tagInfo->impactParameterData();
      for (unsigned int itrk = 0; itrk < tracks.size(); itrk++) {
	//std::cout<<"Track "<<itrk<<" pt="<<tracks[itrk]->pt()<<" d0="<<ip[itrk].ip3d.value()<<" pixelHits="<<tracks[itrk]->hitPattern().pixelLayersWithMeasurement();  
	plots_.allIPErrfromTagInfo->Fill(ip[itrk].ip3d.error()); 
	plots_.allIPSigfromTagInfo->Fill(ip[itrk].ip3d.significance());
	plots_.trackptfromTagInfo->Fill(tracks[itrk]->pt());
	plots_.trackpixelhitsFromTagInfo->Fill(tracks[itrk]->hitPattern().pixelLayersWithMeasurement());    
       	plots_.trackIPfromTagInfo->Fill(ip[itrk].ip3d.value());
      } 
    }

    //============================================================================
    //
    // SV Tag Infos
    //
    //============================================================================

    for(pat::JetCollection::const_iterator bjet=jetsHandle->begin(); bjet!=jetsHandle->end(); ++bjet){
      
      // retrieve the "secondary vertex tag infos"
      // this is the output of the b-tagging reconstruction code
      // and contains secondary vertices in the jets
      //const reco::SecondaryVertexTagInfo* svTagInfo = bjet->tagInfoSecondaryVertex("simpleSecondaryVertexHighEffBJetTags");
      const reco::SecondaryVertexTagInfo* svTagInfo = bjet->tagInfoSecondaryVertex("secondaryVertex");
      
      // count the number of secondary vertices
      plots_.SVnVertices->Fill(svTagInfo->nVertices());
      
      // ignore jets without SV from now on
      if (svTagInfo->nVertices() < 1)
	continue;
      
      // pick the first secondary vertex (the "best" one)
      const reco::Vertex &sv = svTagInfo->secondaryVertex(0);
      
      // and plot number of tracks and chi^2
      plots_.SVnTracks->Fill(sv.tracksSize());
      plots_.SVchi2->Fill(sv.chi2());
      
      // the precomputed transverse distance to the primary vertex
      Measurement1D distance = svTagInfo->flightDistance(0, true);
      plots_.SVdist->Fill(distance.value());
      plots_.SVdistErr->Fill(distance.error());       
      plots_.SVdistSig->Fill(distance.significance());
      
      // the precomputed 3D distance to the primary vertex
      Measurement1D distance3D = svTagInfo->flightDistance(0, false);
      plots_.SVdist3D->Fill(distance3D.value());	       
      plots_.SVdistErr3D->Fill(distance3D.error());      
      plots_.SVdistSig3D->Fill(distance3D.significance());

      Double_t signS(0.);
      if(distance3D.value()>0.){
	signS=1;
      } else {
	signS=-1;	
      }

      // plots_.DSSV->Fill(signS*TMath::Log(1+TMath::Abs(distance3D.significance())));        
      plots_.DSSVHE->Fill(btagSSVHE(svTagInfo));
      plots_.DSSVHE->Fill(-btagNegativeSSVHE(svTagInfo));

      plots_.DSSVHP->Fill(btagSSVHP(svTagInfo));
      plots_.DSSVHP->Fill(-btagNegativeSSVHP(svTagInfo));
      
      // the precomputed direction with respect to the primary vertex
      GlobalVector dir = svTagInfo->flightDirection(0);
      
      // unfortunately CMSSW hsa all kinds of vectors,
      // and sometimes we need to convert them *sigh*
      math::XYZVector dir2(dir.x(), dir.y(), dir.z());
      
      // compute a few variables that we are plotting
      double deltaR = ROOT::Math::VectorUtil::DeltaR(bjet->momentum(), dir2);
      plots_.SVJetdeltaR->Fill(deltaR);
      
      // compute the invariant mass from a four-vector sum
      math::XYZTLorentzVector trackFourVectorSum;
      
      // loop over all tracks in the vertex
      for(reco::Vertex::trackRef_iterator track = sv.tracks_begin(); track != sv.tracks_end(); ++track) {
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> > vec;
	vec.SetPx((*track)->px());
	vec.SetPy((*track)->py());
	vec.SetPz((*track)->pz());
	vec.SetM(0.13957);	// pion mass
	trackFourVectorSum += vec;
      }
      
      // get the invariant mass: sqrt(E² - px² - py² - pz²)
      double vertexMass = trackFourVectorSum.M();
      plots_.SVmass->Fill(vertexMass);
      plots_.SVmass2->Fill(sv.p4().M());


    } //closes loop on bjets
  }


  //------------------------------------------------------------------------------------------------------------------------------
  // MC SV MASS TEMPLATES

  if(ismc_){

    //////////////////////////////////////
    // Tagged Jets and Flavor Separator //
    //////////////////////////////////////
    
    int numBottom = 0, numCharm = 0, numLight = 0;
    int numTags = 0, numGoodJets=0;
    double sumVertexMass = 0.;
    // Loop over the jets and find out which are tagged
    const pat::JetCollection::const_iterator kJetEnd = jetsHandle->end();
    for (pat::JetCollection::const_iterator jetIter = jetsHandle->begin(); kJetEnd != jetIter;  ++jetIter){

      if ( jetIter->pt() < jetPtCut_ ||  std::abs(jetIter->eta()) > jetEtaCut_ || !ZbbUtils::isJetIdOk((*jetIter),"loose") ) continue;
      numGoodJets++;
      
      // Is this jet tagged and does it have a good secondary vertex
      if( !ZbbUtils::isBJet((*jetIter),bTagAlgoWP_) ){
	// This jet isn't tagged
	continue;
      }
      
      reco::SecondaryVertexTagInfo const * svTagInfos = jetIter->tagInfoSecondaryVertex("secondaryVertex");
      if ( svTagInfos->nVertices() <= 0 ) {
	// Given that we are using simple secondary vertex
	// tagging, I don't think this should ever happen.
	// Maybe we should put a counter here just to check.
	continue;
      } // if we have no secondary verticies
      
      // count it
      ++numTags;
      
      // What is the flavor of this jet
      int jetFlavor = std::abs( jetIter->partonFlavour() );
      if (5 == jetFlavor){
	++numBottom;
      } // if bottom 
      else if (4 == jetFlavor){
	++numCharm;
      } // if charm
      else if (1 == jetFlavor || 2 == jetFlavor || 3 == jetFlavor || 21 == jetFlavor ){
	++numLight;
      } // if light flavor
      
      ///////////////////////////
      // Calculate SecVtx Mass //
      ///////////////////////////
      ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > sumVec;
      reco::CompositeCandidate vertexCand;
      reco::Vertex::trackRef_iterator kEndTracks = svTagInfos->secondaryVertex(0).tracks_end();
      for (reco::Vertex::trackRef_iterator trackIter = svTagInfos->secondaryVertex(0).tracks_begin(); trackIter != kEndTracks;  ++trackIter ) {
	const double kPionMass = 0.13957018;
	ROOT::Math::LorentzVector< ROOT::Math::PxPyPzM4D<double> > vec;
	vec.SetPx( (*trackIter)->px() );
	vec.SetPy( (*trackIter)->py() );
	vec.SetPz( (*trackIter)->pz() );
	vec.SetM (kPionMass);
	sumVec += vec;
      } // for trackIter
      sumVertexMass += sumVec.M();
      
      if (2 == numTags){
	// We've got enough.  Stop.
	break;
      } // if we have enough tags
    } // for jet
    
    ////////////////////////
    // General Accounting //
    ////////////////////////
    int numJets = jetsHandle->size();
    plots_.MCjetsVstags->Fill(numGoodJets, numJets);

    // If we don't have any tags, don't bother going on
    if ( ! numTags){
      return;
    }

    ///////////////////////////////////////
    // Calculate average SecVtx mass and //
    // fill appropriate histograms.      //
    ///////////////////////////////////////
    sumVertexMass /= numTags;
    if (1 == numTags){
      // single tag
      if      (numBottom)                 plots_.MCsecvtxMass_b->Fill(sumVertexMass); 
      else if (numCharm)               	  plots_.MCsecvtxMass_c->Fill(sumVertexMass); 
      else if (numLight)               	  plots_.MCsecvtxMass_q->Fill(sumVertexMass); 
      else                             	  plots_.MCsecvtxMass_X->Fill(sumVertexMass);
    } else {				  		      
      // double tags			  
      if      (2 == numBottom)         	  plots_.MCsecvtxMass_bb->Fill(sumVertexMass);
      else if (2 == numCharm)          	  plots_.MCsecvtxMass_cc->Fill(sumVertexMass);
      else if (2 == numLight)          	  plots_.MCsecvtxMass_qq->Fill(sumVertexMass);
      else if (numBottom && numCharm)  	  plots_.MCsecvtxMass_bc->Fill(sumVertexMass);
      else if (numBottom && numLight)  	  plots_.MCsecvtxMass_bq->Fill(sumVertexMass);
      else if (numCharm  && numLight)  	  plots_.MCsecvtxMass_cq->Fill(sumVertexMass);
      else                             	  plots_.MCsecvtxMass_XX->Fill(sumVertexMass);
    } // if two tags
    
    plots_.MCsecvtxMass_All->Fill(sumVertexMass);
  } 
}
	
double  PatBJetTagAnalyzer::btagSSVHE(  const reco::SecondaryVertexTagInfo* svTagInfo ) {
  const edm::RefVector<reco::TrackCollection>& tracks = svTagInfo->selectedTracks();
  for (unsigned int itrk = 0; itrk < tracks.size(); itrk++) {
    const reco::TrackIPTagInfo::TrackIPData & ip = svTagInfo->trackIPData(itrk);
    if (tracks.size()<2) continue;
    if (ip.ip3d.error()<0.) continue; 
    if (ip.ip3d.value()<0.) continue;  
    return TMath::Log(1.+ (ip.ip3d.value()/ip.ip3d.error()));
  }
  return -1.;
}

double  PatBJetTagAnalyzer::btagSSVHP(  const reco::SecondaryVertexTagInfo* svTagInfo ) {
  const edm::RefVector<reco::TrackCollection>& tracks = svTagInfo->selectedTracks(); 
  for (unsigned int itrk = 0; itrk < tracks.size(); itrk++) {
    const reco::TrackIPTagInfo::TrackIPData & ip = svTagInfo->trackIPData(itrk);
    if (tracks.size()<3) continue;
    if (ip.ip3d.error()<0.) continue; 
    if (ip.ip3d.value()<0.) continue;  
    return TMath::Log(1.+ (ip.ip3d.value()/ip.ip3d.error()));
  }
  return -1.;
}

double  PatBJetTagAnalyzer::btagNegativeSSVHE(  const reco::SecondaryVertexTagInfo* svTagInfo ) {
  const edm::RefVector<reco::TrackCollection>& tracks = svTagInfo->selectedTracks();
  for (unsigned int itrk = 0; itrk < tracks.size(); itrk++) {
    const reco::TrackIPTagInfo::TrackIPData & ip = svTagInfo->trackIPData(itrk);
    if (tracks.size()<2) continue;
    if (ip.ip3d.error()<=0.) continue; 
    if (ip.ip3d.value()>0.) continue;  
    return -TMath::Log(1.+ (ip.ip3d.value()/ip.ip3d.error()));
  }
  return -1.;
}

double PatBJetTagAnalyzer::btagNegativeSSVHP(  const reco::SecondaryVertexTagInfo* svTagInfo ) {
  const edm::RefVector<reco::TrackCollection>& tracks = svTagInfo->selectedTracks();
  for (unsigned int itrk = 0; itrk < tracks.size(); itrk++) {
    const reco::TrackIPTagInfo::TrackIPData & ip = svTagInfo->trackIPData(itrk);
    if (tracks.size()<2) continue;
    if (ip.ip3d.error()<=0.) continue; 
    if (ip.ip3d.value()>0.) continue;  
    return -TMath::Log(1.+ (ip.ip3d.value()/ip.ip3d.error()));
  }
  return -1.;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PatBJetTagAnalyzer);
