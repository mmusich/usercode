import FWCore.ParameterSet.Config as cms

analyzePAT = cms.EDAnalyzer("PATAnalyzer",
                            #photonSrc = cms.untracked.InputTag("cleanPatPhotons"),
                            #tauSrc    = cms.untracked.InputTag(""),
                            numEventsNames = cms.untracked.vstring('TotalEventCounter','AfterPVFilterCounter', 'AfterNSFilterCounter', 'AfterPATCounter', 'AfterCandidatesCounter', 'AfterJetsCounter'),
                            electronSrc= cms.untracked.InputTag("patElectrons"),
                            muonSrc    = cms.untracked.InputTag("patMuons"),                                             
                            jetSrc     = cms.untracked.InputTag("patJets"),                                      
                            corrLevels = cms.vstring('Uncorrected','L2Relative', 'L3Absolute'), #, 'L2L3Residual')    ## (to veto for mc)                                  
                            metSrc     = cms.untracked.InputTag("patMETs"),
                            ZmmSrc     = cms.untracked.InputTag("zMMCand"),
                            ZeeSrc     = cms.untracked.InputTag("zEECand"),
                            firstRunNumber = cms.untracked.uint32(0),
                            lastRunNumber  = cms.untracked.uint32(999999)
                            )



