# usage w/ command line options
# cmsRun ZbbEventContentAnalyzer_cfg.py  isMuon=True (or isMuon=False)

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("TEST")

###################################################################
# Setup 'standard' options
###################################################################
options = VarParsing.VarParsing()
options.register('isMuon',True,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool,"True: run on muons; False: run on electrons")
options.register('maxEvents',-1,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int,"Number of events to process (-1 for all)")
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
process.GlobalTag.globaltag  = 'GR_R_42_V20::All'     # for the DATA in 4_2_X
#process.GlobalTag.globaltag = 'GR_R_311_V4::All'  # for the DATA in 4_1_X
#process.GlobalTag.globaltag = 'GR_R_39X_V5::All'  # for the DATA in 3_9_X
#process.GlobalTag.globaltag = 'START39_V8::All'   # for MC

###################################################################
# Additional JP calibration sequence
###################################################################
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
             tag = cms.string("TrackProbabilityCalibration_2D_2011Data_v1_offline"),
             connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU")),
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
             tag = cms.string("TrackProbabilityCalibration_3D_2011Data_v1_offline"),
             connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU"))
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
#from ZbbAnalysis.AnalysisStep.hltrunrangematcher_mu2010_cff import hltpaths_and_runs as  hltpaths_and_runs_mu2010    #2010
from ZbbAnalysis.AnalysisStep.hltrunrangematcher_mu2011_cff import hltpaths_and_runs as  hltpaths_and_runs_mu2011     #2011
process.mu_hltrunrange = hltrunrange.clone(
    HLTRunRangeList=hltpaths_and_runs_mu2011
)

# electrons
#from ZbbAnalysis.AnalysisStep.hltrunrangematcher_ele2010_cff import hltpaths_and_runs as  hltpaths_and_runs_ele2010  #2010
from ZbbAnalysis.AnalysisStep.hltrunrangematcher_ele2011_cff import hltpaths_and_runs as  hltpaths_and_runs_ele2011   #2011
process.ele_hltrunrange = hltrunrange.clone(
    HLTRunRangeList=hltpaths_and_runs_ele2011
)

###################################################################
# PAT basic distributions
###################################################################
from ZbbAnalysis.AnalysisStep.PATAnalyzer_cfi import analyzePAT
process.analyzePat =  analyzePAT.clone(
    jetSrc     = cms.untracked.InputTag("patJets"),
    ##for 2011
    firstRunNumber = 160400,
    lastRunNumber  = 174000
    ) 

###################################################################
# Jet Cleaning
###################################################################
process.cleanPatJets = cms.EDProducer("PATJetCleaner",
                                      src = cms.InputTag("patJets"),
                                      preselection = cms.string('pt > 20.0 && abs(eta) < 2.4'),
                                      checkOverlaps = cms.PSet(ele = cms.PSet(src = cms.InputTag("offlinePrimaryVertexFromZ:patElectronsFromZ"), #Z daughters
                                                                              #src = cms.InputTag("goldenElectrons"),
                                                                              algorithm = cms.string("byDeltaR"),
                                                                              preselection        = cms.string(""),
                                                                              deltaR              = cms.double(0.5),
                                                                              checkRecoComponents = cms.bool(False),
                                                                              pairCut             = cms.string(""),
                                                                              requireNoOverlaps = cms.bool(True)
                                                                              ),
                                                               mu = cms.PSet(src = cms.InputTag("offlinePrimaryVertexFromZ:patMuonsFromZ"), #Z daughters
                                                                             #src = cms.InputTag("goldenMuons"),
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

#process.cleanPatJets = process.cleanPatJets.clone( src = cms.InputTag("patJets") )
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
# PAT basic distributions after cleaning
###################################################################
process.analyzePatAfterCleaning = analyzePAT.clone(                                           
    jetSrc     = cms.untracked.InputTag("cleanPatJets"),
    ##for 2011
    firstRunNumber = 160400,
    lastRunNumber  = 174000
    )

###################################################################
# Event content analysis
###################################################################
from ZbbAnalysis.AnalysisStep.ZbbEventContentAnalyzer_cfi import eventcontentanalyze

### before tagging histograms
process.steps_before_btag = eventcontentanalyze.clone(
    # Double object trigger switch false to run on 2010 / true to run on 2011
    andOr       = cms.bool(True),
    doAllThePlotting = cms.bool(True),
    bTagAlgoWP  = cms.string('SSVHPT'),
    jetSrc      = cms.untracked.InputTag("cleanPatJets"),
    VertexSrc   = cms.untracked.InputTag("offlinePrimaryVertices"),
    muHLTRunRangeList  = hltpaths_and_runs_mu2011,
    eleHLTRunRangeList = hltpaths_and_runs_ele2011,
    minBtags =cms.int32(1),                 #minimum b-tags in event         
    OutfileName = cms.string("OUTLISTTEMPLATE.txt")
    )

### ssvhem
process.finaldistros_ssvhem = eventcontentanalyze.clone(
    # Double object trigger switch false to run on 2010 / true to run on 2011
    andOr       = cms.bool(True),
    bTagAlgoWP  = cms.string('SSVHEM'),
    jetSrc      = cms.untracked.InputTag("cleanPatJets"),
    VertexSrc   = cms.untracked.InputTag("offlinePrimaryVertices"),
    muHLTRunRangeList  = hltpaths_and_runs_mu2011,
    eleHLTRunRangeList = hltpaths_and_runs_ele2011,
    minBtags =cms.int32(1),                 #minimum b-tags in event         
    OutfileName = cms.string("OUTLISTTEMPLATE_ssvhem.txt")
    )

### ssvhpt
process.finaldistros_ssvhpt = eventcontentanalyze.clone(
    # Double object trigger switch false to run on 2010 / true to run on 2011
    andOr       = cms.bool(True),
    bTagAlgoWP  = cms.string('SSVHPT'),
    jetSrc      = cms.untracked.InputTag("cleanPatJets"),
    VertexSrc   = cms.untracked.InputTag("offlinePrimaryVertices"),
    muHLTRunRangeList  = hltpaths_and_runs_mu2011,
    eleHLTRunRangeList = hltpaths_and_runs_ele2011,
    minBtags =cms.int32(1),                 #minimum b-tags in event         
    OutfileName = cms.string("OUTLISTTEMPLATE_ssvhpt.txt")
    )

### csvm
process.finaldistros_csvm = eventcontentanalyze.clone(
    # Double object trigger switch false to run on 2010 / true to run on 2011
    andOr       = cms.bool(True),
    bTagAlgoWP  = cms.string('CSVM'),
    jetSrc      = cms.untracked.InputTag("cleanPatJets"),
    VertexSrc   = cms.untracked.InputTag("offlinePrimaryVertices"),
    muHLTRunRangeList  = hltpaths_and_runs_mu2011,
    eleHLTRunRangeList = hltpaths_and_runs_ele2011,
    minBtags =cms.int32(1),                 #minimum b-tags in event         
    OutfileName = cms.string("OUTLISTTEMPLATE_csvm.txt")
    )

### csvt
process.finaldistros_csvt = eventcontentanalyze.clone(
    # Double object trigger switch false to run on 2010 / true to run on 2011
    andOr       = cms.bool(True),
    bTagAlgoWP  = cms.string('CSVT'),
    jetSrc      = cms.untracked.InputTag("cleanPatJets"),
    VertexSrc   = cms.untracked.InputTag("offlinePrimaryVertices"),
    muHLTRunRangeList  = hltpaths_and_runs_mu2011,
    eleHLTRunRangeList = hltpaths_and_runs_ele2011,
    minBtags =cms.int32(1),                 #minimum b-tags in event         
    OutfileName = cms.string("OUTLISTTEMPLATE_csvt.txt")
    )

# final sequence
process.finaldistros = cms.Sequence(
    process.steps_before_btag+ 
    process.finaldistros_ssvhem+
    process.finaldistros_ssvhpt+
    process.finaldistros_csvm+
    process.finaldistros_csvt
    )

readFiles = cms.untracked.vstring()

###################################################################
# if you want to prefilter with the trigger mask
###################################################################
process.TriggerMaskSequence = cms.Sequence()

process.muonTriggerMaskSequence     = cms.Sequence(process.mu_hltrunrange*~process.ele_hltrunrange)
process.electronTriggerMaskSequence = cms.Sequence(~process.mu_hltrunrange*process.ele_hltrunrange)

if options.isMuon:    
    process.TriggerMaskSequence+=process.muonTriggerMaskSequence
else:
    process.TriggerMaskSequence+=process.electronTriggerMaskSequence

###################################################################

readFiles.extend(FILESOURCETEMPLATE)
outputRootFileName=cms.string('OUTFILETEMPLATE')
            
process.source = cms.Source("PoolSource",
                            fileNames = readFiles ,
                            duplicateCheckMode = cms.untracked.string('checkAllFilesOpened')
                            )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = outputRootFileName
                                   )
###################################################################
# Path 
###################################################################
process.AnalysisSequence = cms.Sequence(process.analyzePat+
                                        process.JetCleaningSeq+
                                        process.btagging*
                                        process.analyzePatAfterCleaning+
                                        process.finaldistros)

process.p = cms.Path( process.TriggerMaskSequence*process.AnalysisSequence )
#process.p = cms.Path( process.AnalysisSequence )









