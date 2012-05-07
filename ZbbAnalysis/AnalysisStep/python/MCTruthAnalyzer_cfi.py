import FWCore.ParameterSet.Config as cms

MCTruthAnalyze = cms.EDAnalyzer("MCTruthAnalyzer",
                                JetCollection = cms.InputTag("cleanPatJets"),
                                genJetSrc = cms.InputTag("patJets:genJets"),
                                src = cms.InputTag("genParticles"),
                                ElectronCollection = cms.InputTag("patElectrons"),
                                MuonCollection = cms.InputTag("patMuons"),
                                ZmmCollection = cms.InputTag('zMMCand'),
                                isMCatNLO = cms.untracked.bool(False)
                                )

