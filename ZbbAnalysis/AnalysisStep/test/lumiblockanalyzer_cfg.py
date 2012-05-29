import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("ZbbAnalysis.AnalysisStep.muon_ZbbSkimDec22ReReco_PAT397_08Apr11_cff")
process.load("ZbbAnalysis.AnalysisStep.lumiblockanalyzer_cfi")

process.TFileService=cms.Service('TFileService',
                                 fileName=cms.string('MyPlots.root')
                                 )




process.p = cms.Path(process.lumiblockanalysis)
