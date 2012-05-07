import FWCore.ParameterSet.Config as cms

CalculateAcceptance = cms.EDAnalyzer("ZjetsAcceptanceCalculator",
                                     GenSrc = cms.InputTag("genParticles"),
                                     ElectronCollection = cms.InputTag("patElectrons"),
                                     MuonCollection = cms.InputTag("patMuons"),
                                     JetCollection = cms.InputTag("cleanPatJets"),
                                     genJetSrc = cms.InputTag("patJets:genJets"),
                                     ZmmCollection = cms.InputTag('zMMCand'),
                                     ZeeCollection = cms.InputTag('zEECand'),
                                     unLockDefaultCuts = cms.bool(False),  # bool to unlock default acceptance cuts
                                     jetEtaCut   = cms.double(2.1),
                                     muonEtaCut  = cms.double(2.1),
                                     eleEtaCut   = cms.double(2.5),
                                     jetPtCut    = cms.double(25.),
                                     genJetPtCut = cms.double(0.),
                                     muonPtCut   = cms.double(20.),
                                     elePtCut    = cms.double(25.),
                                     minMassCut  = cms.double(60.),
                                     maxMassCut  = cms.double(120.),
                                     dRLeptonMatch = cms.double(0.3),
                                     dRJetMatch    = cms.double(0.5),
                                     useStatus3forMuons = cms.bool(False), 	    
                                     useStatus3forElectrons = cms.bool(True), 	    
                                     doFSRCorrectionForMuons = cms.bool(False), 
                                     doFSRCorrectionForElectrons = cms.bool(False),                                     
                                     BparticlePtCut  = cms.double(0.),
                                     isMCatNLO     = cms.untracked.bool(False),
                                     applyPUCorr   = cms.untracked.bool(True),
                                     verbose       = cms.untracked.bool(False),
                                     useClopperPearsonErrors = cms.untracked.bool(True),
                                     saveNTuple    = cms.untracked.bool(True),
                                     PartonLevel = cms.untracked.bool(False),
                                     OutfileName = cms.string('CorrectionFactors_Zmm.txt'),
                                     DecayChainSelection = cms.untracked.string('Z_m'),
                                     EffbcMCSSVHE =  cms.VPSet(cms.PSet(ptrange=cms.untracked.vdouble(0,9999),effb_barrel=cms.untracked.double(1),effb_forward=cms.untracked.double(1),effc_barrel=cms.untracked.double(1),effc_forward=cms.untracked.double(1))),   
                                     EffbcMCSSVHP =  cms.VPSet(cms.PSet(ptrange=cms.untracked.vdouble(0,9999),effb_barrel=cms.untracked.double(1),effb_forward=cms.untracked.double(1),effc_barrel=cms.untracked.double(1),effc_forward=cms.untracked.double(1)))   
                                     )
