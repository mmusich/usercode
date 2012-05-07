import FWCore.ParameterSet.Config as cms

analyzeBJets = cms.EDAnalyzer("PatBJetTagAnalyzer",
                              # input collections
                              isMC        = cms.bool(False),   
                              bTagAlgoWP  = cms.string('SSVHEM'), 
                              primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                              beamSpot = cms.InputTag("offlineBeamSpot"),
                              tracks = cms.InputTag("generalTracks"),
                              jets = cms.InputTag("patJets"),
                              
                              #jet cuts
                              jetPtCut = cms.double(25),
                              jetEtaCut = cms.double(2.1),
                              
                              #track cuts
                              maxDeltaR = cms.double(0.5),
                              minPt = cms.double(1.0),
                              minPixelHits = cms.uint32(2),
                              minTotalHits = cms.uint32(8),
                              
                              # selecting the second-highest signed IP (i.e. "TrackCountingHighEff")
                              nThTrack = cms.uint32(2)
                              )
