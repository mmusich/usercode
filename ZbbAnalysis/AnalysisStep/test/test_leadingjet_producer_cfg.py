import FWCore.ParameterSet.Config as cms
process = cms.Process("LEADINGPT")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

###################################################################
# PoolSource
###################################################################
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_2_1_epd.root'
                                                              )
                            )

###################################################################
# Jet Cleaning
###################################################################
process.cleanPatJets = cms.EDProducer("PATJetCleaner",
                                      src = cms.InputTag("patJets"),
                                      preselection = cms.string('pt > 20.0 && abs(eta) < 2.4'),
                                      checkOverlaps = cms.PSet(ele = cms.PSet(src       = cms.InputTag("offlinePrimaryVertexFromZ:patElectronsFromZ"), #Z daughters
                                                                              #src       = cms.InputTag("goldenElectrons"),
                                                                              algorithm = cms.string("byDeltaR"),
                                                                              preselection        = cms.string(""),
                                                                              deltaR              = cms.double(0.5),
                                                                              checkRecoComponents = cms.bool(False),
                                                                              pairCut             = cms.string(""),
                                                                              requireNoOverlaps = cms.bool(True)
                                                                              ),
                                                               mu = cms.PSet(src       = cms.InputTag("offlinePrimaryVertexFromZ:patMuonsFromZ"), #Z daughters
                                                                             #src       = cms.InputTag("goldenMuons"),
                                                                             algorithm = cms.string("byDeltaR"),
                                                                             preselection        = cms.string(""),
                                                                             deltaR              = cms.double(0.5),
                                                                             checkRecoComponents = cms.bool(False),
                                                                             pairCut             = cms.string(""),
                                                                             requireNoOverlaps = cms.bool(True)
                                                                             ),
                                                               ),
                                      finalCut = cms.string(''),
                                      )

###################################################################
# Leading jets producer
###################################################################
from ZbbAnalysis.Tools.LeadingJetsProducer_cfi import LeadingJetsProducer
process.leadingJetsProducer = LeadingJetsProducer.clone(
    verbose      = cms.bool(False),
    jets         = cms.InputTag("cleanPatJets"),
    genJets      = cms.InputTag("patJets:genJets"),
    primaryVertices= cms.InputTag("offlinePrimaryVertices"),
    etaCut       = cms.double(2.1),
    betaCut      = cms.double(0.15),
    betastarCut  = cms.double(10000),
    collectionLength= cms.int32(1)
    )

###################################################################
# PAT basic distributions
###################################################################
from ZbbAnalysis.AnalysisStep.PATAnalyzer_cfi import analyzePAT
process.analyzePat =  analyzePAT.clone(
    jetSrc = cms.untracked.InputTag("leadingJetsProducer:patLeadingJets")
)

###################################################################
# output file
###################################################################
process.TFileService = cms.Service("TFileService",
                                   fileName= cms.string("myTestPatAnalyze.root")
                                   )

###################################################################
# PoolOutput
###################################################################
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('/tmp/musich/myOutputFile.root'),
                               outputCommands = cms.untracked.vstring('keep *')
                               )

###################################################################
# Paths and endpaths
###################################################################
process.p = cms.Path(process.cleanPatJets*process.leadingJetsProducer*process.analyzePat)
process.e = cms.EndPath(process.out)
