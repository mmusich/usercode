import FWCore.ParameterSet.Config as cms

process = cms.Process('TEST')
# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/Mu/Mu_All/MergedOutputFile_1_1_9WV.root'))

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

from ZbbAnalysis.AnalysisStep.hltrunrangematcher_cfi import hltrunrange
# muons
from ZbbAnalysis.AnalysisStep.hltrunrangematcher_mu2010_cff import hltpaths_and_runs as  hltpaths_and_runs_mu2010
process.mu_hltrunrange = hltrunrange.clone(
    HLTRunRangeList=hltpaths_and_runs_mu2010
)
# electrons
from ZbbAnalysis.AnalysisStep.hltrunrangematcher_ele2010_cff import hltpaths_and_runs as  hltpaths_and_runs_ele2010
process.ele_hltrunrange = hltrunrange.clone(
    HLTRunRangeList=hltpaths_and_runs_ele2010
)



process.p1 = cms.Path(process.mu_hltrunrange)
process.p2 = cms.Path(process.ele_hltrunrange)

process.schedule = cms.Schedule(process.p1,process.p2)

# Trigger report via MessageLogger service
# process.MessageLogger.destinations = ['HLTreport_Mu_All.txt']
# from HLTrigger.HLTanalyzers.hlTrigReport_cfi import hlTrigReport
# process.hltReport = hlTrigReport.clone()

# process.endpath = cms.EndPath(process.hltReport) 
# process.MessageLogger.categories.append("HLTrigReport")
