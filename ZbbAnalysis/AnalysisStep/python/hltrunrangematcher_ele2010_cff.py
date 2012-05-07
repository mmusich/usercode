import FWCore.ParameterSet.Config as cms

## Electron and Egamma Paths
hltpaths_and_runs = cms.VPSet()
hltpaths_and_runs.extend([
    cms.PSet(hltpath=cms.untracked.string("HLT_Photon10_L1R"),runrange=cms.untracked.vuint32(132440,137028)),  # should impose a cut at 15 GeV by hand
    cms.PSet(hltpath=cms.untracked.string("HLT_Photon15_Cleaned_L1R"),runrange=cms.untracked.vuint32(138564,140401)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Ele15_SW_CaloEleId_L1R"),runrange=cms.untracked.vuint32(141956,144114)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Ele17_SW_CaloEleId_L1R"),runrange=cms.untracked.vuint32(146428,147116)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Ele17_SW_TightEleId_L1R"),runrange=cms.untracked.vuint32(147146,148102)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v1"),runrange=cms.untracked.vuint32(148783,149063)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2"),runrange=cms.untracked.vuint32(149181,149442))
    ])





