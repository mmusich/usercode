import FWCore.ParameterSet.Config as cms

## Muon Paths
hltpaths_and_runs = cms.VPSet()
hltpaths_and_runs.extend([
    cms.PSet(hltpath=cms.untracked.string("HLT_Mu3"),runrange=cms.untracked.vuint32(132440,139980)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Mu5"),runrange=cms.untracked.vuint32(140058,140401)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Mu7"),runrange=cms.untracked.vuint32(141956,144114)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Mu9"),runrange=cms.untracked.vuint32(146428,147116)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Mu11"),runrange=cms.untracked.vuint32(147146,148102)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Mu15_v1"),runrange=cms.untracked.vuint32(148783,149442))
    ])


