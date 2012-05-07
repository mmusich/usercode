# usage w/ command line options
# cmsRun ZbbEventContentAnalyzer_cfg.py  isMuon=True (or isMuon=False)

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("TEST")

# setup 'standard'  options
options = VarParsing.VarParsing()
# setup any defaults you want
options.register('isMuon',True, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool,"True: run on muons; False: run on electrons")
options.parseArguments()


# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#process.GlobalTag.globaltag = 'GR_R_311_V4::All'  # for the DATA in 4_1_X
process.GlobalTag.globaltag = 'GR_R_39X_V5::All'   # for the DATA in 3_9_X
#process.GlobalTag.globaltag = 'START39_V8::All'   # for MC

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
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
#from ZbbAnalysis.AnalysisStep.hltrunrangematcher_mu2011_cff import hltpaths_and_runs as  hltpaths_and_runs_mu2011

process.mu_hltrunrange = hltrunrange.clone(
    HLTRunRangeList=hltpaths_and_runs_mu2010
)

# electrons
from ZbbAnalysis.AnalysisStep.hltrunrangematcher_ele2010_cff import hltpaths_and_runs as  hltpaths_and_runs_ele2010
#from ZbbAnalysis.AnalysisStep.hltrunrangematcher_ele2011_cff import hltpaths_and_runs as  hltpaths_and_runs_ele2011

process.ele_hltrunrange = hltrunrange.clone(
    HLTRunRangeList=hltpaths_and_runs_ele2010
)


###################################################################
# PAT basic distributions
###################################################################
from ZbbAnalysis.AnalysisStep.PATAnalyzer_cfi import analyzePAT
process.analyzePat =  analyzePAT.clone(
    ##for 2010
    firstRunNumber = 132440,
    lastRunNumber  = 139980
    ##for 2011
    #firstRunNumber = 160410,
    #lastRunNumber  = 999999
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
    #+ process.cleanPatJets + 
    #+ process.cleanPatJetsNoPU
     )
    #+ process.jetFilter                  # eventually filter on candidate jets 
    )

###################################################################
# PAT basic distributions after cleaning
###################################################################
process.analyzePatAfterCleaning = process.analyzePat.clone(                                           
    jetSrc     = cms.untracked.InputTag("cleanPatJetsPF")                                     
    )

###################################################################
# Event content analysis
###################################################################
from ZbbAnalysis.AnalysisStep.ZbbEventContentAnalyzer_cfi import eventcontentanalyze
process.finaldistros = eventcontentanalyze.clone(
    # Double object trigger switch false to run on 2010 / true to run on 2011
    andOr       = cms.bool(False),                                
    jetSrc      = cms.untracked.InputTag("cleanPatJetsPF"),
    muHLTRunRangeList  = hltpaths_and_runs_mu2010,
    eleHLTRunRangeList = hltpaths_and_runs_ele2010,
    )

###################################################################
# Source files
###################################################################
muonSource    = ['rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/Mu/Mu_All/MergedOutputFile_1_1_9WV.root']  # 2010

electronSource= ['rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/EG/MergedOutputFile_1_1_kgn.root',         # 2010
                 'rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/Electron/MergedOutputFile_1_1_6ue.root']

readFiles = cms.untracked.vstring()

###################################################################
# if you want to prefilter with the trigger mask
###################################################################
process.TriggerMaskSequence = cms.Sequence()

process.muonTriggerMaskSequence     = cms.Sequence(process.mu_hltrunrange*~process.ele_hltrunrange)
process.electronTriggerMaskSequence = cms.Sequence(~process.mu_hltrunrange*process.ele_hltrunrange)
###################################################################

if options.isMuon:    
    readFiles.extend(muonSource)
    process.TriggerMaskSequence+=process.muonTriggerMaskSequence
    outputRootFileName=cms.string('analyzePAT_DATA2010_muon.root')
    process.finaldistros.OutfileName = cms.string("MuonList.txt")
else:
    readFiles.extend(electronSource)
    process.TriggerMaskSequence+=process.electronTriggerMaskSequence
    outputRootFileName=cms.string('analyzePAT_DATA2010_electron.root')
    process.finaldistros.OutfileName = cms.string("ElectronList.txt")
        
process.source = cms.Source("PoolSource",
                            fileNames = readFiles ,
                            duplicateCheckMode = cms.untracked.string('checkAllFilesOpened')
                            )

process.TFileService = cms.Service("TFileService",
                                   fileName = outputRootFileName
                                   )

###################################################################
# check if events in the list
###################################################################
process.checkList = cms.EDFilter("EventListFilter",
                             Inputfile=cms.untracked.string("list.txt")
                             )

###################################################################
# Path 
###################################################################
process.AnalysisSequence = cms.Sequence(process.analyzePat+
                                        process.JetCleaningSeq+
                                        process.analyzePatAfterCleaning+
                                        process.finaldistros)

process.p = cms.Path( process.TriggerMaskSequence*process.AnalysisSequence )









