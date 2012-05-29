import FWCore.ParameterSet.Config as cms

LeadingJetsProducer = cms.EDProducer('LeadingJetsProducer',
                                     jets         = cms.InputTag("patJets"),
                                     genJets      = cms.InputTag("patJets:genJets"),
                                     primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                                     verbose      = cms.bool(False),
                                     etaCut       = cms.double(5.0),
                                     betaCut      = cms.double(0.15),  # put here very _LOW_  value to be transparent
                                     betastarCut  = cms.double(0.85),  # put here very _HIGH_ value to be transparent
                                     collectionLength = cms.int32(1)
                                     )
