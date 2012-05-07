import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'GR_R_39X_V5::All' # for the DATA
#process.GlobalTag.globaltag = 'START39_V8::All' # for MC


###################################################################
# Source files
###################################################################
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/Mu/Mu_All/MergedOutputFile_1_1_9WV.root" 
                                                              )
                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

###################################################################
# Check duplicates (for DATA)
###################################################################
process.options = cms.untracked.PSet(
  duplicateCheckMode = cms.untracked.string('checkAllFilesOpened')
)

###################################################################
# Trigger Filter 
###################################################################
### Trigger path FIRED in the event
# import HLTrigger.HLTfilters.hltHighLevel_cfi 
# process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone() 
# process.hltFilter.HLTPaths = [ "HLT_DoubleMu3_v*" ]
# process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
# process.hltFilter.throw  = cms.bool(False)

###################################################################
# HLT run range matcher
###################################################################

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


###################################################################
# PAT basic distributions
###################################################################
process.analyzePat = cms.EDAnalyzer("PATAnalyzer",
  #photonSrc = cms.untracked.InputTag("cleanPatPhotons"),
  #tauSrc    = cms.untracked.InputTag(""),                                     
  electronSrc= cms.untracked.InputTag("patElectronsWithTrigger"),
  muonSrc    = cms.untracked.InputTag("patMuonsWithTrigger"),                                             
  jetSrc     = cms.untracked.InputTag("patJetsPF"),                                      
  corrLevels = cms.vstring('Uncorrected','L2Relative', 'L3Absolute'), #, 'L2L3Residual')    ## (to veto for mc)                                  
  metSrc     = cms.untracked.InputTag("patMETs"),
  ZmmSrc     = cms.untracked.InputTag("zMMCand"),
  ZeeSrc     = cms.untracked.InputTag("zEECand"),                                    
  firstRunNumber = cms.untracked.uint32(136030),
  lastRunNumber  = cms.untracked.uint32(140160)                                    
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('analyzePAT_DATA_Mu.root')
                                   )

###################################################################
# Jet Cleaning
###################################################################
process.cleanPatJetsPF = cms.EDProducer("PATJetCleaner",
                                      src = cms.InputTag("patJetsPF"),
                                      preselection = cms.string('pt > 20.0 && abs(eta) < 2.4'),
                                      checkOverlaps = cms.PSet(ele = cms.PSet(src       = cms.InputTag("goldenElectrons"),
                                                                              algorithm = cms.string("byDeltaR"),
                                                                              preselection        = cms.string(""),
                                                                              deltaR              = cms.double(0.5),
                                                                              checkRecoComponents = cms.bool(False),
                                                                              pairCut             = cms.string(""),
                                                                              requireNoOverlaps = cms.bool(True)
                                                                              ),
                                                               mu = cms.PSet(src       = cms.InputTag("goldenMuons"),
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

#process.cleanPatJets = process.cleanPatJetsPF.clone( src = cms.InputTag("patJets") )
#process.cleanPatJetsNoPU = process.cleanPatJetsPF.clone( src = cms.InputTag("patJetsNoPU") )


###################################################################
# Jet number filter... ?
###################################################################
process.jetFilter = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("cleanPatJets"),
                                  minNumber = cms.uint32(2),
                                 )

process.JetCleaningSeq = cms.Sequence( 
    (
    process.cleanPatJetsPF 
#    + process.cleanPatJets + 
#    + process.cleanPatJetsNoPU
     )
    # + process.jetFilter                  # eventually filter on candidate jets 
    )


###################################################################
# PAT basic distributions after cleaning
###################################################################
process.analyzePatAfterCleaning = cms.EDAnalyzer("PATAnalyzer",
  #photonSrc = cms.untracked.InputTag("cleanPatPhotons"),
  #tauSrc    = cms.untracked.InputTag(""),                                     
  electronSrc= cms.untracked.InputTag("patElectronsWithTrigger"),
  muonSrc    = cms.untracked.InputTag("patMuonsWithTrigger"),                                             
  jetSrc     = cms.untracked.InputTag("cleanPatJetsPF"),                                      
  corrLevels = cms.vstring('Uncorrected','L2Relative', 'L3Absolute'), #, 'L2L3Residual')    ## (to veto for mc)                                  
  metSrc     = cms.untracked.InputTag("patMETs"),
  ZmmSrc     = cms.untracked.InputTag("zMMCand"),
  ZeeSrc     = cms.untracked.InputTag("zEECand")                                                
)


###################################################################
# Ntuplizer (not used if EDM ntuple not produced)
###################################################################
process.load("ZbbAnalysis.AnalysisStep.ZJJNtuplizer_cff")
process.edmNtuplesOut.fileName = cms.untracked.string('file:/tmp/musich/ZbbAnalysisCandNtuples_DATA_Merged.root')
#process.edmNtuplesOut.fileName = cms.untracked.string('file:/tmp/musich/ZbbAnalysisCandNtuples_Mu_ByRun.root')
#process.edmNtuplesOut.fileName = cms.untracked.string('file:/tmp/musich/ZbbAnalysisCandNtuples_Ele_ByRun.root')
#process.edmNtuplesOut.fileName = cms.untracked.string('file:/tmp/musich/ZbbAnalysisCandNtuples_Mu.root')
#process.edmNtuplesOut.fileName = cms.untracked.string('file:/tmp/musich/ZbbAnalysisCandNtuples_MC_ZbbToLL.root')
#process.edmNtuplesOut.fileName = cms.untracked.string('file:/tmp/musich/ZbbAnalysisCandNtuples_MC_ZccToLL.root')
#process.edmNtuplesOut.fileName = cms.untracked.string('file:/tmp/musich/ZbbAnalysisCandNtuples_MC_TTJets.root')
#process.edmNtuplesOut.fileName = cms.untracked.string('file:/tmp/musich/ZbbAnalysisCandNtuples_MC_DYJets.root')
#process.edmNtuplesOut.fileName = cms.untracked.string('file:/tmp/musich/ZbbAnalysisCandNtuples_Ele.root')
#process.edmNtuplesOut.fileName = cms.untracked.string('file:/tmp/musich/ZbbAnalysisCandNtuples_Egamma.root')

###################################################################
# Path and Endpath
###################################################################
process.p = cms.Path(process.mu_hltrunrange*~process.ele_hltrunrange*(process.analyzePat+
                                                                      process.JetCleaningSeq+
                                                                      process.analyzePatAfterCleaning)

#process.p = cms.Path(process.mu_hltrunrange*~process.ele_hltrunrange*(process.analyzePat+
#                                                                      process.JetCleaningSeq+
#                                                                      process.ZMMNtuplizer+
#                                                                      process.ZEENtuplizer+
#                                                                      process.ZLLNtuplizer+
#                                                                      process.JetNtuplizer+
#                                                                      process.METNtuplizer)
                     )

# not used if EMD nutple not produced
#process.endPath = cms.EndPath(process.edmNtuplesOut)

#process.schedule = cms.Schedule(
#    process.p,
#    process.endPath
#    )

###################################################################
# Trigger report via MessageLogger service
###################################################################
#process.MessageLogger.destinations = ['HLTreport_Mu_All.txt']
#from HLTrigger.HLTanalyzers.hlTrigReport_cfi import hlTrigReport
#process.hltReport = hlTrigReport.clone()

#process.endpath = cms.EndPath(process.hltReport) 
#process.MessageLogger.categories.append("HLTrigReport")

