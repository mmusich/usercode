import FWCore.ParameterSet.Config as cms

## Electron Paths
hltpaths_and_runs = cms.VPSet()
hltpaths_and_runs.extend([
    cms.PSet(hltpath=cms.untracked.vstring("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*","hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter","hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter"),runrange=cms.untracked.vuint32(160410,161216)),
    cms.PSet(hltpath=cms.untracked.vstring("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*","hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter","hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter"),runrange=cms.untracked.vuint32(161217,163268)),
    cms.PSet(hltpath=cms.untracked.vstring("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*","hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter","hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter"),runrange=cms.untracked.vuint32(163269,165120)),
    cms.PSet(hltpath=cms.untracked.vstring("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*","hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter","hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter"),runrange=cms.untracked.vuint32(165121,165969)),
    cms.PSet(hltpath=cms.untracked.vstring("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*","hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter","hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter"),runrange=cms.untracked.vuint32(165970,167038)),
    cms.PSet(hltpath=cms.untracked.vstring("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*","hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter","hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter"),runrange=cms.untracked.vuint32(167039,170248)),
    cms.PSet(hltpath=cms.untracked.vstring("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*","hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter","hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter"),runrange=cms.untracked.vuint32(170249,171049)),
    cms.PSet(hltpath=cms.untracked.vstring("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*","hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter","hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter"),runrange=cms.untracked.vuint32(171050,173235)),
    cms.PSet(hltpath=cms.untracked.vstring("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*","hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter","hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter"),runrange=cms.untracked.vuint32(173236,199999)),
    ])

