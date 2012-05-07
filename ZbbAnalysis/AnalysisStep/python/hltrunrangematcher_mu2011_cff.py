import FWCore.ParameterSet.Config as cms

## Muon Paths
hltpaths_and_runs = cms.VPSet()
hltpaths_and_runs.extend([
    cms.PSet(hltpath=cms.untracked.vstring("HLT_DoubleMu6_v1","HLT_DoubleMu6_v*","hltDiMuonL3PreFiltered6","hltDiMuonL3PreFiltered6"),runrange=cms.untracked.vuint32(160410,163268)),
    cms.PSet(hltpath=cms.untracked.vstring("HLT_DoubleMu7_v2","HLT_DoubleMu7_v*","hltDiMuonL3PreFiltered7","hltDiMuonL3PreFiltered7"),runrange=cms.untracked.vuint32(163269,165120)),
    cms.PSet(hltpath=cms.untracked.vstring("HLT_Mu13_Mu8_v2","HLT_Mu13_Mu8_v*","hltSingleMu13L3Filtered13","hltDiMuonL3PreFiltered8"),runrange=cms.untracked.vuint32(165121,170248)),
    cms.PSet(hltpath=cms.untracked.vstring("HLT_Mu13_Mu8_v3","HLT_Mu13_Mu8_v*","hltSingleMu13L3Filtered13","hltDiMuonL3PreFiltered8"),runrange=cms.untracked.vuint32(167039,170248)),
    cms.PSet(hltpath=cms.untracked.vstring("HLT_Mu13_Mu8_v4","HLT_Mu13_Mu8_v*","hltSingleMu13L3Filtered13","hltDiMuonL3PreFiltered8"),runrange=cms.untracked.vuint32(167039,170248)),
    cms.PSet(hltpath=cms.untracked.vstring("HLT_Mu13_Mu8_v6","HLT_Mu13_Mu8_v*","hltSingleMu13L3Filtered13","hltDiMuonL3PreFiltered8"),runrange=cms.untracked.vuint32(170249,173235)),
    cms.PSet(hltpath=cms.untracked.vstring("HLT_Mu13_Mu8_v7","HLT_Mu13_Mu8_v*","hltSingleMu13L3Filtered13","hltDiMuonL3p5PreFiltered8"),runrange=cms.untracked.vuint32(173236,199999))
    ])

