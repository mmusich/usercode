import FWCore.ParameterSet.Config as cms

zvertexproducer = cms.EDProducer('ZVertexProducer',
                                 VertexSrc   = cms.untracked.InputTag("offlinePrimaryVertices"),
                                 ZmmSrc      = cms.untracked.InputTag("zMMCand"),
                                 ZeeSrc      = cms.untracked.InputTag("zEECand")
                                 )
