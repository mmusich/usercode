import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


# run the input file through the end;
# for a limited number of events, replace -1 with the desired number 
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load( "SimGeneral.HepPDTESSource.pythiapdt_cfi" )

process.source = cms.Source( "PoolSource",
                             fileNames = cms.untracked.vstring(
'file:/tmp/emiglior/store/mc/Winter10/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/E7TeV_ProbDist_2010Data_BX156_START39_V8-v1/0078/C0D6C43D-F615-E011-ADBA-E0CB4E19F9A1.root'
#			     'file:/tmp/emiglior/store/mc/Winter10/ZbbToLL_M-40_PtB1-15_TuneZ2_7TeV-madgraph-pythia6/AODSIM/E7TeV_ProbDist_2010Data_BX156_START39_V8-v1/0016/0A4182F2-681D-E011-873D-003048D45FE2.root'
			     )
                           )
	      
# FileService is mandatory, as the following analyzer module 
# will want it, to create output histogram file
# 
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("test.root")
)

# the analyzer itself - empty parameter set 
#
process.TestLHEanalyzer = cms.EDAnalyzer( "LHEAnalyzer" )

process.p1 = cms.Path( process.TestLHEanalyzer )

