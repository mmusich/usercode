import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("TEST")

#setup 'standard'  options
options = VarParsing.VarParsing()
# setup any defaults you want
options.register('Sample',
                 "DY", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string,         # string, int, or float
                 "Sample to be processed (DY,Zbb,Zcc,tt)")
options.parseArguments()

# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 5000


process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#process.GlobalTag.globaltag = 'GR_R_39X_V5::All' # for the DATA
process.GlobalTag.globaltag = 'START41_V0::All' # for MC

###################################################################
# Source files
###################################################################
#DYSource= ['rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/DYJetsToLL_TuneZ2/MergedOutputFile_1_1_hFC.root']
#ZbbSource=['rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/ZbbToLL/MergedOutputFile_1_1_Q6N.root']
#ZccSource=['rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/ZccToLL/MergedOutputFile_1_1_dh4.root']
TTSource=['rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/TTJets_TuneZ2/MergedOutputFile_1_1_85V.root']

# new production
DYSource=[
'rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimMCWinter10_PAT397_22Jun11_V2/DYJetsToLL_TuneZ2_split/MergedOutputFile_1_1_9np.root',
'rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimMCWinter10_PAT397_22Jun11_V2/DYJetsToLL_TuneZ2_split/MergedOutputFile_2_1_xzP.root',
'rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimMCWinter10_PAT397_22Jun11_V2/DYJetsToLL_TuneZ2_split/MergedOutputFile_3_1_GNk.root'
]
ZbbSource=['rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimMCWinter10_PAT397_22Jun11_V2/ZbbToLL/MergedOutputFile_1_1_QYJ.root']
ZccSource=['rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimMCWinter10_PAT397_22Jun11_V2/ZccToLL/MergedOutputFile_1_1_13S.root']

readFiles = cms.untracked.vstring()

process.source = cms.Source("PoolSource",
                            fileNames = readFiles,
                            duplicateCheckMode=cms.untracked.string('checkAllFilesOpened')
                            )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

##############################
# b-tag calibration (currently not used)
##############################

#Data measurements from Fall10
#process.load ("RecoBTag.PerformanceDB.BTagPerformanceDB1011")       # BtagPerformanceESProducer
#process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1011")   # PoolDBESSource
                             
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

# from ZbbAnalysis.AnalysisStep.hltrunrangematcher_cfi import hltrunrange
# # muons
# from ZbbAnalysis.AnalysisStep.hltrunrangematcher_mu2010_cff import hltpaths_and_runs as  hltpaths_and_runs_mu2010
# process.mu_hltrunrange = hltrunrange.clone(
#     HLTRunRangeList=hltpaths_and_runs_mu2010
# )
# # electrons
# from ZbbAnalysis.AnalysisStep.hltrunrangematcher_ele2010_cff import hltpaths_and_runs as  hltpaths_and_runs_ele2010
# process.ele_hltrunrange = hltrunrange.clone(
#     HLTRunRangeList=hltpaths_and_runs_ele2010
# )

###################################################################
# HF filter...  (for DY - status3+ status2 daughter veto)
###################################################################

from ZbbAnalysis.AnalysisStep.HeavyFlavorFilter_cfi import heavyflavorfilter

process.myheavyflavorfilter= heavyflavorfilter.clone(
    src= cms.InputTag("prunedGen"),
    status2 = cms.bool(True),
    status3 = cms.bool(True),
    hDaughterVeto = cms.bool(True),
    zDaughterVeto = cms.bool(False),
    ptcut=cms.double(15)
    )

###################################################################
# PAT basic distributions
###################################################################
from ZbbAnalysis.AnalysisStep.PATAnalyzer_cfi import analyzePAT
process.analyzePat =  analyzePAT.clone() 

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
process.analyzePatAfterCleaning = analyzePAT.clone(                                           
    jetSrc     = cms.untracked.InputTag("cleanPatJetsPF")                                    
    )

###################################################################
# Import b-efficiency scaling factors
###################################################################
from ZbbAnalysis.AnalysisStep.sfbtables_cff import ssvhem_2010
from ZbbAnalysis.AnalysisStep.sfbtables_cff import ssvhpt_2010
from ZbbAnalysis.AnalysisStep.sfbtables_cff import ssvhem_2011
from ZbbAnalysis.AnalysisStep.sfbtables_cff import ssvhpt_2011

###################################################################
# Event content analysis
###################################################################
from ZbbAnalysis.AnalysisStep.ZbbEventContentAnalyzer_cfi import eventcontentanalyze
process.finaldistros = eventcontentanalyze.clone(
    isMC        = cms.bool(True),
    jetSrc      = cms.untracked.InputTag("cleanPatJetsPF"),
    SFbSSVHEMPtRangeList = ssvhem_2010,
    SFbSSVHPTPtRangeList = ssvhpt_2010,
    applybEffCalibration = cms.bool(True),
    bEffCalibrationMethod = cms.string('DUMMY'),
    applyPUcorrection = cms.bool(False),
    PUFileName = cms.string('DUMMY')
    #  muHLTRunRangeList  = hltpaths_and_runs_mu2010,
    #  eleHLTRunRangeList = hltpaths_and_runs_ele2010,
    )

###################################################################
# Definition of analysis sequence
###################################################################

process.AnalysisSequence   = cms.Sequence()

## sequence with hf pt>15 
process.AnalysisWithFilter = cms.Sequence(~process.myheavyflavorfilter*(process.analyzePat+
                                                                        process.JetCleaningSeq+
                                                                        process.analyzePatAfterCleaning+
                                                                        process.finaldistros
                                                                        )
                                          )

## sequence without hf pt>15 
process.AnalysisNoFilter   = cms.Sequence(process.analyzePat+
                                          process.JetCleaningSeq+
                                          process.analyzePatAfterCleaning+
                                          process.finaldistros
                                          )

###################################################################
# switches
###################################################################

if options.Sample=="DY":    
    readFiles.extend(DYSource)
    outputRootFileName=cms.string('analyzePAT_MC_DYJets.root')
    process.AnalysisSequence+=process.AnalysisWithFilter
elif options.Sample=="Zbb":
    readFiles.extend(ZbbSource)
    outputRootFileName=cms.string('analyzePAT_MC_ZbbToLL.root')
    process.AnalysisSequence+=process.AnalysisNoFilter
elif options.Sample=="Zcc":
    readFiles.extend(ZccSource)
    outputRootFileName=cms.string('analyzePAT_MC_ZccToLL.root')
    process.AnalysisSequence+=process.AnalysisNoFilter
elif options.Sample=="tt":
    readFiles.extend(TTSource)
    outputRootFileName=cms.string('analyzePAT_MC_TTJets.root')
    process.AnalysisSequence+=process.AnalysisNoFilter
else :
    print "============================================="
    print "%msg-w: error while exectuing ZbbEventContent"
    print "        Unrecognized sample"
    print "        Please choose among DY,Zbb,Zcc or tt"
    print "============================================="
    outputRootFileName=cms.string('fake.root')

###################################################################
# output file
###################################################################

process.TFileService = cms.Service("TFileService", fileName = outputRootFileName)

###################################################################
# path
###################################################################

process.p = cms.Path(process.AnalysisSequence)



