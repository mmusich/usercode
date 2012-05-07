#################
#PATMuons analysis
##################

import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'START39_V8::All'

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
# ZbbLL
#     'rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/ZbbToLL/MergedOutputFile_1_1_Q6N.root' 
# ZccLL
      'rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/ZccToLL/MergedOutputFile_1_1_dh4.root'
    )
                            )

# local variable specifying the number of evts to print
nEvtsToDebug = 100

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(nEvtsToDebug)
)

process.TFileService=cms.Service('TFileService',
                                 fileName=cms.string('MCTruth_Zbb_2010_Plots_FromPAT_03Mar11-filter.root')
                                 )


# configure the module for pythia-style particles listing
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
  useMessageLogger = cms.untracked.bool(True),
  maxEventsToPrint = cms.untracked.int32(nEvtsToDebug),
  printVertex = cms.untracked.bool(False),
  printP4 = cms.untracked.bool(False),
  printStatus = cms.untracked.bool(True),
  printIndex = cms.untracked.bool(True),
  printPtEtaPhi = cms.untracked.bool(False),
  src = cms.InputTag("genParticles")
)

process.load("ZbbAnalysis.AnalysisStep.leptfromdecayfilter_cfi")
#process.leptfromdecayfilter.DecayChainSelection = cms.untracked.string('b>m')
#process.leptfromdecayfilter.DecayChainSelection = cms.untracked.string('bc>m')
process.leptfromdecayfilter.DecayChainSelection = cms.untracked.string('c>m')

####################################################
#MC Analyzer: spectra of MC Truth Analyzer
####################################################
#process.MCAnalyzer = cms.EDAnalyzer("MCTruthAnalyzer",
#                                    src = cms.InputTag("genParticles"),
#                                    JetCollection = cms.InputTag("patJets"),
#                                    ElectronCollection = cms.InputTag("patElectronsWithTrigger"),
#                                    MuonCollection = cms.InputTag("patMuonsWithTrigger"),
#                                    ZmmCollection = cms.InputTag('zMMCand')
#                                    )

# configure messages
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring('cout', 'info_messages'),
    destinations = cms.untracked.vstring('cout', 'info_messages'), ## .log automatically appended
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING'),
        noLineBreaks = cms.untracked.bool(True)
    ),
    info_messages = cms.untracked.PSet(
        noLineBreaks = cms.untracked.bool(True),
        INFO = cms.untracked.PSet(limit = cms.untracked.int32(1000000)),
        threshold = cms.untracked.string('INFO')
   ),
   debugModules = cms.untracked.vstring('leptfromdecayfilter')                                     
)


# 
process.p = cms.Path(process.printTree+process.leptfromdecayfilter)

#process.p = cms.Path(process.leptfromdecayfilter)



