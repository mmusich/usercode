import FWCore.ParameterSet.Config as cms

## Electron and Egamma Paths
runs = cms.VPSet()
runs.extend([
    cms.PSet(runrange=cms.untracked.vuint32(132440,137028)),  
    cms.PSet(runrange=cms.untracked.vuint32(138564,140401)),
    cms.PSet(runrange=cms.untracked.vuint32(141956,144114)),
    cms.PSet(runrange=cms.untracked.vuint32(146428,147116)),
    cms.PSet(runrange=cms.untracked.vuint32(147146,148102)),
    cms.PSet(runrange=cms.untracked.vuint32(148783,149063)),
    cms.PSet(runrange=cms.untracked.vuint32(149181,149442))
    ])






