import FWCore.ParameterSet.Config as cms

## Muon Paths
runs = cms.VPSet()
runs.extend([
    cms.PSet(runrange=cms.untracked.vuint32(132440,139980)),
    cms.PSet(runrange=cms.untracked.vuint32(140058,140401)),
    cms.PSet(runrange=cms.untracked.vuint32(141956,144114)),
    cms.PSet(runrange=cms.untracked.vuint32(146428,147116)),
    cms.PSet(runrange=cms.untracked.vuint32(147146,148102)),
    cms.PSet(runrange=cms.untracked.vuint32(148783,149442))
    ])


