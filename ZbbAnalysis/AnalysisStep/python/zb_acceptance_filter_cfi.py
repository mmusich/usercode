import FWCore.ParameterSet.Config as cms
# specify 'Zee' or 'Zmm'
zb_acceptance_filter = cms.EDFilter("ZbAcceptanceFilter",
                                    genPSrc    = cms.untracked.InputTag("genParticles"), 
                                    flavLept   = cms.untracked.string('DUMMY'),  # available Zmm,Zee,all 
                                    etaMuMax = cms.untracked.double(999),
                                    ptMuMin  = cms.untracked.double(0),
                                    etaEleMax = cms.untracked.double(999),
                                    ptEleMin  = cms.untracked.double(0),   
                                    flavJet = cms.untracked.string('DUMMY'),     # available 'b','c','all' cases
                                    etaGenJetMax = cms.untracked.double(999),
                                    ptGenJetMin  = cms.untracked.double(0),
                                    ngoodGenJet = cms.untracked.int32(0),
                                    isExclusive = cms.untracked.bool(False),     # if True -> == ngoodGenJet, if False -> >= ngoodGenJet
                                    recoJetSrc = cms.untracked.InputTag("patJets"),
                                    genJetSrc  = cms.untracked.InputTag("patJets:genJets"),
                                    electronSrc = cms.untracked.InputTag("offlinePrimaryVertexFromZ:patElectronsFromZ"),
                                    muonSrc = cms.untracked.InputTag("offlinePrimaryVertexFromZ:patMuonsFromZ"),
                                    )
