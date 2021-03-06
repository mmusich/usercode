import FWCore.ParameterSet.Config as cms

totalKinematicsFilter = cms.EDFilter('TotalKinematicsFilter',
                                     src             = cms.InputTag("genParticles"),
                                     tolerance       = cms.double(5),
                                     verbose         = cms.untracked.bool(False)                                   
                                     )


