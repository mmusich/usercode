import FWCore.ParameterSet.Config as cms

ztautaufilter = cms.EDFilter('ZTauTauFilter',
                             src= cms.InputTag("genParticles")
                             )


