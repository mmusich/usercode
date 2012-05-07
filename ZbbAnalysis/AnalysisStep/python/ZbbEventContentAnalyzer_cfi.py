import FWCore.ParameterSet.Config as cms

eventcontentanalyze = cms.EDAnalyzer("ZbbEventContentAnalyzer",
                                     OutfileName = cms.string("interestingevents.txt"),
                                     doAllThePlotting = cms.bool(False),                                # true will make the plots for all categories
                                     isMC        = cms.bool(False),                                   # MC switch
                                     applyMCweight = cms.bool(True),                                  # apply mc weight for SHERPA/ aMC@NLO
                                     isMCatNLO     = cms.bool(False),                                 # enables rule for NLO reweighting 
                                     Debug       = cms.bool(False),                                   # enable debug mode
                                     andOr       = cms.bool(False),                                   # if false => single object match (2010) / if true => double object match (2011)
                                     bTagAlgoWP  = cms.string('SSVHEM'),                              # b-tag method applied
                                     unLockDefaultCuts = cms.bool(False),                             # bool to unlock default acceptance cuts         
                                     jetEtaCut   = cms.double(2.1),
                                     muonEtaCut  = cms.double(2.5),
                                     eleEtaCut   = cms.double(2.1),
                                     jetPtCut    = cms.double(25),
                                     muonPtCut   = cms.double(20),
                                     elePtCut    = cms.double(25),
                                     genPSrc     = cms.untracked.InputTag("genParticles"),            # MC truth
                                     genJetSrc   = cms.untracked.InputTag("patJets:genJets"), 
                                     electronSrc = cms.untracked.InputTag("patElectrons"),            # electrons
                                     muonSrc     = cms.untracked.InputTag("patMuons"),                # muons                                   
                                     jetSrc      = cms.untracked.InputTag("patJets"),                 # jets                                             
                                     metSrc      = cms.untracked.InputTag("patMETs"),                 # met
                                     VertexSrc   = cms.untracked.InputTag("offlinePrimaryVertices"),  # vertices
                                     ZVertexSrc  = cms.untracked.InputTag("offlinePrimaryVertexFromZ:ZVertexProducerZVertex"),  # Z vertex
                                     ZmmSrc      = cms.untracked.InputTag("zMMCand"),                 # Z->mm candidates
                                     ZeeSrc      = cms.untracked.InputTag("zEECand"),                 # Z->ee candidates
                                     TriggerEventSrc      = cms.untracked.InputTag("patTriggerEvent"),# Trigger
                                     muHLTRunRangeList    = cms.VPSet(cms.PSet(hltpath=cms.untracked.vstring("*"),runrange=cms.untracked.vuint32(0,9999))), # Trigger to Match
                                     eleHLTRunRangeList   = cms.VPSet(cms.PSet(hltpath=cms.untracked.vstring("*"),runrange=cms.untracked.vuint32(0,9999))), # Trigger to Match
                                     applybEffCalibration = cms.bool(False),   # apply b-eff calibration
                                     bEffCalibrationMethod = cms.string('DUMMY'), bMistagCalibrationMethod = cms.string('DUMMY'), # apply b-eff scale
                                     applyLeptonEfficiency = cms.bool(False),  # apply letpon-eff calibration
                                     bSFkFactor  = cms.int32(0), #SFb multiplicative factor (for uncertainty)
                                     lEffkFactor = cms.int32(0), #light Eff multiplicative factor (for uncertainty)
                                     EffbcMCPtRangeList    = cms.VPSet(cms.PSet(ptrange=cms.untracked.vdouble(0,9999),effb_barrel=cms.untracked.double(1),effb_forward=cms.untracked.double(1),effc_barrel=cms.untracked.double(1),effc_forward=cms.untracked.double(1))),      # MC b/c-eff 
                                     applyPUcorrection    = cms.bool(False), # apply PU-reweighting
                                     minBtags =cms.int32(1)
                                     )
                                     
