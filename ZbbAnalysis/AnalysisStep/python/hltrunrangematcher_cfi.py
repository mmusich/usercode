import FWCore.ParameterSet.Config as cms

hltrunrange = cms.EDFilter('HLTRunRangeMatcher',
                           HLTRunRangeList=cms.VPSet(cms.PSet(hltpath=cms.untracked.vstring("*","*","*","*"),runrange=cms.untracked.vuint32(0,999999))
                                                     )
                           )
