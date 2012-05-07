import FWCore.ParameterSet.Config as cms

###############################################################
#MC filter: select events with muons from b semileptonic decays 
###############################################################
leptfromdecayfilter = cms.EDFilter("LeptFromDecayFilter",
                                   src = cms.InputTag("genParticles"),
                                   JetCollection = cms.InputTag("cleanPatJetsPF"),
                                   ElectronCollection = cms.InputTag("patElectronsWithTrigger"),
                                   MuonCollection = cms.InputTag("patMuonsWithTrigger"),
                                   ZmmCollection = cms.InputTag('zMMCand'),
                                   DecayChainSelection = cms.untracked.string('b>m')
                                   )
