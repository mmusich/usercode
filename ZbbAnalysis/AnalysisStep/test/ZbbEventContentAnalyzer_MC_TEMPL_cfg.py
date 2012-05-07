import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("TEST")

###################################################################
# Setup 'standard' options
###################################################################
options = VarParsing.VarParsing()
options.register('Sample',
                 "Zl", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string,         # string, int, or float
                 "Sample to be processed (Zl,Zb5f,Zb5fSherpa,Zb5faMCatNLO,Zc,tt,Ztautau,zz,wz)")
options.register('maxEvents',
                 -1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to process (-1 for all)")
options.parseArguments()

###################################################################
# Messages
###################################################################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

###################################################################
# Standard loads
###################################################################
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

###################################################################
# Global tag
###################################################################
#process.GlobalTag.globaltag = 'GR_R_39X_V5::All' # for the DATA
process.GlobalTag.globaltag = 'START42_V13::All'  # for MC

###################################################################
# Additional JP calibration sequence
###################################################################
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),        
             tag = cms.string("TrackProbabilityCalibration_2D_2011_v1_mc"),
             connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU")),
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
             tag = cms.string("TrackProbabilityCalibration_3D_2011_v1_mc"),
             connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU"))
    )

###################################################################
# Source files
###################################################################
readFiles = cms.untracked.vstring()
readFiles.extend(FILESOURCETEMPLATE)
outputRootFileName=cms.string('OUTFILETEMPLATE')

process.source = cms.Source("PoolSource",
                            fileNames = readFiles,
                            duplicateCheckMode=cms.untracked.string('checkAllFilesOpened')
                            )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
                             
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
# HLT run range matcher (for data only)
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
# Z -> tau tau filter
###################################################################
from ZbbAnalysis.AnalysisStep.ZTauTauFilter_cfi import ztautaufilter

process.zttfilter = ztautaufilter.clone(
    src= cms.InputTag("genParticles")
    ) 

###################################################################
# HF filter (for DY - status3+ status2 daughter veto)
###################################################################

from ZbbAnalysis.AnalysisStep.HeavyFlavorFilter_cfi import heavyflavorfilter

process.b_heavyflavorfilter= heavyflavorfilter.clone(
    src= cms.InputTag("genParticles"),
    status2 = cms.bool(True),
    status3 = cms.bool(True),
    bOn = cms.bool(True),
    cOn = cms.bool(False),
    hDaughterVeto = cms.bool(True),
    zDaughterVeto = cms.bool(False),
    ptcut=cms.double(0)
    )

process.notHadr_b_heavyflavorfilter= heavyflavorfilter.clone(
    src= cms.InputTag("genParticles"),
    status2 = cms.bool(True),
    status3 = cms.bool(True),
    bOn = cms.bool(True),
    cOn = cms.bool(False),
    hDaughterVeto = cms.bool(False),
    zDaughterVeto = cms.bool(False),
    ptcut=cms.double(0)
    )

process.c_heavyflavorfilter= heavyflavorfilter.clone(
    src= cms.InputTag("genParticles"),
    status2 = cms.bool(True),
    status3 = cms.bool(True),
    cOn = cms.bool(True),
    bOn = cms.bool(False),
    hDaughterVeto = cms.bool(True),
    zDaughterVeto = cms.bool(False),
    ptcut=cms.double(0)
    )

process.bc_heavyflavorfilter= heavyflavorfilter.clone(
    src= cms.InputTag("genParticles"),
    status2 = cms.bool(True),
    status3 = cms.bool(True),
    cOn = cms.bool(True),
    bOn = cms.bool(True),
    hDaughterVeto = cms.bool(True),
    zDaughterVeto = cms.bool(False),
    ptcut=cms.double(0)
    )

###################################################################
# PAT basic distributions
###################################################################
from ZbbAnalysis.AnalysisStep.PATAnalyzer_cfi import analyzePAT
process.analyzePat =  analyzePAT.clone(
    jetSrc = cms.untracked.InputTag("patJets")
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

#process.cleanPatJets     = process.cleanPatJets.clone( src = cms.InputTag("patJets") )
#process.cleanPatJetsNoPU = process.cleanPatJets.clone( src = cms.InputTag("patJetsNoPU") )

###################################################################
# Jet number filter... ?
###################################################################
process.jetFilter = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("cleanPatJets"),
                                 minNumber = cms.uint32(2),
                                 )

process.JetCleaningSeq = cms.Sequence( 
    (
    process.cleanPatJets 
    #+ process.cleanPatJets + 
    #+ process.cleanPatJetsNoPU
     )
    #+ process.jetFilter                  # eventually filter on candidate jets 
    )

###################################################################
# PAT basic distributions after cleaning
###################################################################
process.analyzePatAfterCleaning = analyzePAT.clone(                                           
    jetSrc     = cms.untracked.InputTag("cleanPatJets")                                    
    )

###################################################################
# Event content analysis
###################################################################

###################################################################
# MC Truth Analyzer
###################################################################
from ZbbAnalysis.AnalysisStep.MCTruthAnalyzer_cfi import MCTruthAnalyze

process.MCTruthAnalyzer = MCTruthAnalyze.clone(
    isMCatNLO = cms.untracked.bool(False)
    )

####################################
# SFb from DB                      #
#                                  #
# RecoBTag/PerformanceDB V00-04-11 #
#   REQUIRED IN CMSSW_4_2_1        #
#                                  #
####################################
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")
process.load("RecoBTag.PerformanceDB.BTagPerformanceDB1107")

###################################################################
# Import b-efficiency scaling factors
###################################################################
from ZbbAnalysis.AnalysisStep.effbctables_cff import mc_effbc_ssvhem_2011
from ZbbAnalysis.AnalysisStep.effbctables_cff import mc_effbc_ssvhpt_2011
from ZbbAnalysis.AnalysisStep.effbctables_cff import mc_effbc_csvm_2011
from ZbbAnalysis.AnalysisStep.effbctables_cff import mc_effbc_csvt_2011

# other corrections available ('MISTAGSSVHPT', 'MISTAGSSVHEM','MISTAGCSVM','MISTAGCSVT'),
#                             ('BTAGSSVHPT'  , 'BTAGSSVHEM'  ,'BTAGCSVM'  ,'BTAGCSVT')

from ZbbAnalysis.AnalysisStep.ZbbEventContentAnalyzer_cfi import eventcontentanalyze

### all before tagging histograms
process.steps_before_btag = eventcontentanalyze.clone(
    isMC        = cms.bool(True),
    doAllThePlotting = cms.bool(True),
    bTagAlgoWP  = cms.string('SSVHPT'),
    jetSrc      = cms.untracked.InputTag("cleanPatJets"),
    applybEffCalibration = cms.bool(BSELTEMPLATE),
    bEffCalibrationMethod = cms.string('BTAGSSVHPT'),
    bMistagCalibrationMethod = cms.string('MISTAGSSVHPT'),    
    EffbcMCPtRangeList = mc_effbc_ssvhpt_2011,
    applyPUcorrection = cms.bool(PUTEMPLATE),
    applyLeptonEfficiency = cms.bool(TNPTEMPLATE),
    minBtags =cms.int32(1),           #minimum b-tags in event
    OutfileName = cms.string("OUTLISTTEMPLATE.txt")
    )

### ssvhem
process.finaldistros_ssvhem = eventcontentanalyze.clone(
    isMC        = cms.bool(True),
    bTagAlgoWP  = cms.string('SSVHEM'),
    jetSrc      = cms.untracked.InputTag("cleanPatJets"),
    applybEffCalibration = cms.bool(BSELTEMPLATE),
    bEffCalibrationMethod = cms.string('BTAGSSVHEM'),
    bMistagCalibrationMethod = cms.string('MISTAGSSVHEM'),    
    EffbcMCPtRangeList = mc_effbc_ssvhem_2011,
    applyPUcorrection = cms.bool(PUTEMPLATE),
    applyLeptonEfficiency = cms.bool(TNPTEMPLATE),
    minBtags =cms.int32(1),           #minimum b-tags in event
    OutfileName = cms.string("OUTLISTTEMPLATE_ssvhem.txt")
    )

### ssvhpt
process.finaldistros_ssvhpt = eventcontentanalyze.clone(
    isMC        = cms.bool(True),
    bTagAlgoWP  = cms.string('SSVHPT'),
    jetSrc      = cms.untracked.InputTag("cleanPatJets"),
    applybEffCalibration = cms.bool(BSELTEMPLATE),
    bEffCalibrationMethod = cms.string('BTAGSSVHPT'),
    bMistagCalibrationMethod = cms.string('MISTAGSSVHPT'),    
    EffbcMCPtRangeList = mc_effbc_ssvhpt_2011,
    applyPUcorrection = cms.bool(PUTEMPLATE),
    applyLeptonEfficiency = cms.bool(TNPTEMPLATE),
    minBtags =cms.int32(1),          #minimum b-tags in event
    OutfileName = cms.string("OUTLISTTEMPLATE_ssvhpt.txt")
    )

### csvm
process.finaldistros_csvm = eventcontentanalyze.clone(
    isMC        = cms.bool(True),
    bTagAlgoWP  = cms.string('CSVM'),
    jetSrc      = cms.untracked.InputTag("cleanPatJets"),
    applybEffCalibration = cms.bool(BSELTEMPLATE),
    bEffCalibrationMethod = cms.string('BTAGCSVM'),
    bMistagCalibrationMethod = cms.string('MISTAGCSVM'),    
    EffbcMCPtRangeList = mc_effbc_csvm_2011,
    applyPUcorrection = cms.bool(PUTEMPLATE),
    applyLeptonEfficiency = cms.bool(TNPTEMPLATE),
    minBtags =cms.int32(1),           #minimum b-tags in event
    OutfileName = cms.string("OUTLISTTEMPLATE_csvm.txt")
    )

### csvt
process.finaldistros_csvt = eventcontentanalyze.clone(
    isMC        = cms.bool(True),
    bTagAlgoWP  = cms.string('CSVT'),
    jetSrc      = cms.untracked.InputTag("cleanPatJets"),
    applybEffCalibration = cms.bool(BSELTEMPLATE),
    bEffCalibrationMethod = cms.string('BTAGCSVT'),
    bMistagCalibrationMethod = cms.string('MISTAGCSVT'),    
    EffbcMCPtRangeList = mc_effbc_csvt_2011,
    applyPUcorrection = cms.bool(PUTEMPLATE),
    applyLeptonEfficiency = cms.bool(TNPTEMPLATE),
    minBtags =cms.int32(1),          #minimum b-tags in event
    OutfileName = cms.string("OUTLISTTEMPLATE_csvt.txt")
    )

##################################################################
# hack to force NLO weights rule
##################################################################
if options.Sample=="Zb5faMCatNLO":
    process.finaldistros_ssvhem.isMCatNLO = cms.bool(True)    
    process.finaldistros_ssvhpt.isMCatNLO = cms.bool(True)    
    process.MCTruthAnalyzer.isMCatNLO = cms.untracked.bool(True)

# final sequence
process.finaldistros = cms.Sequence(
    process.steps_before_btag+
    process.finaldistros_ssvhem+
    process.finaldistros_ssvhpt+
    process.finaldistros_csvm+
    process.finaldistros_csvt
    )

###################################################################
# Filter on the generated acceptance
###################################################################
from ZbbAnalysis.AnalysisStep.zb_acceptance_filter_cfi import zb_acceptance_filter 

#filter for b
process.zplusbAccfilter = zb_acceptance_filter.clone(
    flavLept  = cms.untracked.string('all'),   # available Zmm,Zee,all
    etaMuMax  = cms.untracked.double(2.1),
    ptMuMin   = cms.untracked.double(17.),
    etaEleMax = cms.untracked.double(2.5),
    ptEleMin  = cms.untracked.double(17.),                                     
    flavJet   = cms.untracked.string('b'),     # available b,c cases
    etaGenJetMax = cms.untracked.double(3.5),
    ptGenJetMin  = cms.untracked.double(15),
    ngoodGenJet  = cms.untracked.int32(1),
    isExclusive  = cms.untracked.bool(False),  # if True -> == ngoodGenJet, if False -> >= ngoodGenJet
    jetSrc = cms.untracked.InputTag("patJets")
    )

# filter for c
process.zpluscAccfilter = zb_acceptance_filter.clone(
    flavLept  = cms.untracked.string('all'),   # available Zmm,Zee,all
    etaMuMax  = cms.untracked.double(2.1),
    ptMuMin   = cms.untracked.double(17.),
    etaEleMax = cms.untracked.double(2.5),
    ptEleMin  = cms.untracked.double(17.),                                     
    flavJet   = cms.untracked.string('c'),     # available b,c cases   
    etaGenJetMax = cms.untracked.double(3.5),
    ptGenJetMin  = cms.untracked.double(15),
    ngoodGenJet  = cms.untracked.int32(1),
    isExclusive  = cms.untracked.bool(False),  # if True -> == ngoodGenJet, if False -> >= ngoodGenJet
    jetSrc = cms.untracked.InputTag("patJets")
    )

####################################################################
# Define the b-tag squences for offline reconstruction 
####################################################################
process.btagging = cms.Sequence()

process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
process.MyJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
                                                       process.j2tParametersVX,
                                                       jets = cms.InputTag("cleanPatJets")
                                                       )

process.load("RecoBTag.ImpactParameter.impactParameter_cff")
process.MyImpactParameterTagInfos = process.impactParameterTagInfos.clone(
    jetTracks = "MyJetTracksAssociatorAtVertex"
    )

# re-run track counting b-tagging algorithm he
process.MyTrackCountingHighEffBJetTags = process.trackCountingHighEffBJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("MyImpactParameterTagInfos"))
    )

# re-run track counting b-tagging algorithm hp
process.MyTrackCountingHighPurBJetTags = process.trackCountingHighPurBJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("MyImpactParameterTagInfos"))
    )

# re-run track counting b-tagging algorithm jp
process.MyJetProbabilityBJetTags = process.jetProbabilityBJetTags.clone(
    tagInfos = cms.VInputTag(cms.InputTag("MyImpactParameterTagInfos"))
    )

# re-run btagging sequence
process.btagging += cms.Sequence(process.MyJetTracksAssociatorAtVertex*
                                 (process.MyImpactParameterTagInfos *
                                  (process.MyTrackCountingHighEffBJetTags +
                                   process.MyTrackCountingHighPurBJetTags +
                                   process.MyJetProbabilityBJetTags 
                                   )
                                  )
                                 )

###################################################################
# Definition of analysis sequence
###################################################################
process.AnalysisSequence   = cms.Sequence()

process.keepIfB = cms.Sequence(process.b_heavyflavorfilter*(process.analyzePat+
                                                            process.JetCleaningSeq+
                                                            process.btagging*
                                                            process.MCTruthAnalyzer+
                                                            process.analyzePatAfterCleaning+
                                                            process.finaldistros
                                                            )
                               )

process.keepIfC = cms.Sequence(~process.b_heavyflavorfilter*process.c_heavyflavorfilter*(process.analyzePat+
                                                                                         process.JetCleaningSeq+
                                                                                         process.btagging*
                                                                                         process.MCTruthAnalyzer+
                                                                                         process.analyzePatAfterCleaning+
                                                                                         process.finaldistros
                                                                                         )
                               )

process.dropIfBC = cms.Sequence(~process.bc_heavyflavorfilter*(process.analyzePat+
                                                               process.JetCleaningSeq+
                                                               process.btagging*
                                                               process.MCTruthAnalyzer+
                                                               process.analyzePatAfterCleaning+
                                                               process.finaldistros
                                                               )
                                )

process.keepIfNotHadrB = cms.Sequence(process.notHadr_b_heavyflavorfilter*(process.analyzePat+
                                                                           process.JetCleaningSeq+
                                                                           process.btagging*
                                                                           process.MCTruthAnalyzer+
                                                                           process.analyzePatAfterCleaning+
                                                                           process.finaldistros
                                                                           )
                                      )

## sequence without HeavyFlavourFilter 
process.AnalysisNoFilter = cms.Sequence(process.analyzePat+
                                        process.JetCleaningSeq+
                                        process.btagging*
                                        process.MCTruthAnalyzer+
                                        process.analyzePatAfterCleaning+
                                        process.finaldistros
                                        )

##############################################################################
# Patch for sum(E) violation in PYTHIA/MADGRAPH DY samples
# https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/1489.html
##############################################################################
process.load("ZbbAnalysis.Tools.TotalKinematicsFilter_cfi")

if options.Sample=="Zl" or options.Sample=="Zb5f" or options.Sample=="Zc":
    process.AnalysisSequence+=process.totalKinematicsFilter

###################################################################
# switches
###################################################################
if options.Sample=="Zl":
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.dropIfBC
elif options.Sample=="Zb5f":
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.keepIfB
elif options.Sample=="Zc":
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.keepIfC
elif options.Sample=="tt":
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.AnalysisNoFilter
elif options.Sample=="Zb5fSherpa":
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.keepIfNotHadrB
elif options.Sample=="Zb5faMCatNLO":
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.AnalysisNoFilter
elif options.Sample=="zz":
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.AnalysisNoFilter
elif options.Sample=="wz":
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.AnalysisNoFilter
elif options.Sample=="Ztautau":
    process.AnalysisSequence+=process.zttfilter
    process.AnalysisSequence+=process.AnalysisNoFilter
else :
    print "============================================="
    print "%msg-w: error while exectuing ZbbEventContent"
    print "        Unrecognized sample"
    print "        Please choose among Zl,Zb5f,Zb5fSherpa,Zb5faMCatNLO,Zc or tt"
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



