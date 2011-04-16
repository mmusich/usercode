import FWCore.ParameterSet.Config as cms

process = cms.Process('TEST')
# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/Mu/Mu_All/MergedOutputFile_1_1_9WV.root'))

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.hltrunrange = cms.EDFilter('HLTRunRangeMatcher',
                                   MuonHLTRunRangeList=cms.VPSet(
    ## Muon Paths
    cms.PSet(hltpath=cms.untracked.string("HLT_Mu3"),runrange=cms.untracked.vuint32(132440,139980)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Mu5"),runrange=cms.untracked.vuint32(140058,140401)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Mu7"),runrange=cms.untracked.vuint32(141956,144114)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Mu9"),runrange=cms.untracked.vuint32(146428,147116)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Mu11"),runrange=cms.untracked.vuint32(147146,148102)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Mu15_v1"),runrange=cms.untracked.vuint32(148783,149442))
    ),
                                   ElectronHLTRunRangeList=cms.VPSet(
    ## Electron and Egamma Paths
    cms.PSet(hltpath=cms.untracked.string("HLT_Photon10_L1R"),runrange=cms.untracked.vuint32(132440,137028)),  # should impose a cut at 15 GeV by hand
    cms.PSet(hltpath=cms.untracked.string("HLT_Photon15_Cleaned_L1R"),runrange=cms.untracked.vuint32(138564,140401)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Ele15_SW_CaloEleId_L1R"),runrange=cms.untracked.vuint32(141956,144114)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Ele17_SW_CaloEleId_L1R"),runrange=cms.untracked.vuint32(146428,147116)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Ele17_SW_TightEleId_L1R"),runrange=cms.untracked.vuint32(147146,148102)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v1"),runrange=cms.untracked.vuint32(148783,149063)),
    cms.PSet(hltpath=cms.untracked.string("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2"),runrange=cms.untracked.vuint32(149181,149442))
    ),
                                   isMuChannel=cms.bool(True)
                                   )

process.p = cms.Path(process.hltrunrange)

# Trigger report via MessageLogger service
process.MessageLogger.destinations = ['HLTreport_Mu_All.txt']
from HLTrigger.HLTanalyzers.hlTrigReport_cfi import hlTrigReport
process.hltReport = hlTrigReport.clone()

process.endpath = cms.EndPath(process.hltReport) 
process.MessageLogger.categories.append("HLTrigReport")
