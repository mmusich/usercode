import FWCore.ParameterSet.Config as cms
# specify 'Zee' or 'Zmm'
lhezllfilter = cms.EDFilter("LHEZllFilter",                            
                            ZllDecaySelection = cms.untracked.string('DUMMY'),
                            etaMax = cms.untracked.double(999),
                            ptMins = cms.untracked.vdouble(0.,0.)
                            )
