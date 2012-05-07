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
process.GlobalTag.globaltag  = 'GR_R_42_V20::All'  # for the DATA in 4_2_X
#process.GlobalTag.globaltag = 'GR_R_311_V4::All'  # for the DATA in 4_1_X
#process.GlobalTag.globaltag = 'GR_R_39X_V5::All'  # for the DATA in 3_9_X
#process.GlobalTag.globaltag = 'START39_V8::All'   # for MC

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

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
    firstRunNumber = 160410,
    lastRunNumber  = 174000,
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
    jetSrc     = cms.untracked.InputTag("cleanPatJets"),
    ##for 2011
    firstRunNumber = 160410,
    lastRunNumber  = 174000,
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
    minBtags =cms.int32(1)             #minimum b-tags in event         
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
    minBtags =cms.int32(1)             #minimum b-tags in event         
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
    minBtags =cms.int32(1)             #minimum b-tags in event         
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
    minBtags =cms.int32(1)             #minimum b-tags in event         
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
    minBtags =cms.int32(1)             #minimum b-tags in event         
    )
#####################################################################

# final sequence
process.finaldistros = cms.Sequence(
    process.steps_before_btag+
    process.finaldistros_ssvhem+
    process.finaldistros_ssvhpt+
    process.finaldistros_csvm+
    process.finaldistros_csvt 
    )

###################################################################
# Source files
###################################################################

electronSource=[
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_1_1_4IW.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_2_1_Uce.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_3_1_Qhk.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_4_1_Nuq.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_5_1_Khn.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_6_1_dg5.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_7_1_bJK.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_8_1_cbt.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_9_2_LrS.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_10_3_pL3.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_11_1_zYu.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_12_2_efg.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_13_1_jwC.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_14_2_UEa.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_15_2_kEJ.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_16_2_vl4.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_17_1_bNs.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_18_2_YvL.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_19_1_Vmj.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_20_1_MkM.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_21_1_pZf.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_22_1_iQD.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_23_1_K2F.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_24_1_w6G.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_25_1_X2o.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_26_1_KuX.root'
    ] 

muonSource=[
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_1_1_Uui.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_2_1_lUC.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_3_1_iNw.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_4_1_Z8m.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_5_1_44d.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_6_1_vyd.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_7_1_717.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_8_1_6bu.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_9_1_2K8.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_10_1_lBg.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_11_1_s3U.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_12_1_SfB.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_13_2_hi7.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_14_2_7Yq.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_15_2_ffA.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_16_3_JAL.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_17_1_Osd.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_18_2_T0M.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_19_2_jlN.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_20_2_zac.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_21_1_WYP.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_22_3_3PF.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_23_1_1xg.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_24_1_1GM.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_25_2_s8x.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_26_2_ZdI.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_27_1_WbB.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_28_2_x4D.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_29_3_ld4.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_30_2_XwM.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_31_2_LRz.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_32_4_je7.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_33_3_LKF.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_34_1_KNH.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_35_2_QAN.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_36_2_lzV.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_37_2_z2G.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_38_1_neq.root'
    ]

readFiles = cms.untracked.vstring()

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
# if you want to prefilter with the trigger mask
###################################################################
process.TriggerMaskSequence = cms.Sequence()

process.muonTriggerMaskSequence     = cms.Sequence(process.mu_hltrunrange*~process.ele_hltrunrange)
process.electronTriggerMaskSequence = cms.Sequence(~process.mu_hltrunrange*process.ele_hltrunrange)
###################################################################

if options.isMuon:    
    readFiles.extend(muonSource)
    process.TriggerMaskSequence+=process.muonTriggerMaskSequence
    outputRootFileName=cms.string('analyzePAT_DATA2011_muon.root')
    process.finaldistros_ssvhem.OutfileName = cms.string("MuonList2011_ssvhem.txt")
    process.finaldistros_ssvhpt.OutfileName = cms.string("MuonList2011_ssvhpt.txt")
    process.finaldistros_csvm.OutfileName = cms.string("MuonList2011_cvsm.txt")
    process.finaldistros_csvt.OutfileName = cms.string("MuonList2011_csvt.txt")
else:
    readFiles.extend(electronSource)
    process.TriggerMaskSequence+=process.electronTriggerMaskSequence
    outputRootFileName=cms.string('analyzePAT_DATA2011_electron.root')
    process.finaldistros_ssvhem.OutfileName = cms.string("ElectronList2011_ssvhem.txt")
    process.finaldistros_ssvhpt.OutfileName = cms.string("ElectronList2011_ssvhpt.txt")
    process.finaldistros_csvm.OutfileName = cms.string("ElectronList2011_csvm.txt")
    process.finaldistros_csvt.OutfileName = cms.string("ElectronList2011_csvt.txt")
        
process.source = cms.Source("PoolSource",
                            fileNames = readFiles ,
                            duplicateCheckMode = cms.untracked.string('checkAllFilesOpened')
                            )

process.TFileService = cms.Service("TFileService",
                                   fileName = outputRootFileName
                                   )

###################################################################
# check if 
###################################################################
process.checkList = cms.EDFilter("EventListFilter",
                                 Inputfile=cms.untracked.string("eventlist.txt")
                                 )

###################################################################
# Path 
###################################################################
process.AnalysisSequence = cms.Sequence(process.analyzePat+
                                        process.JetCleaningSeq+
                                        process.btagging*
                                        process.analyzePatAfterCleaning+
                                        process.finaldistros)

process.p = cms.Path(process.TriggerMaskSequence*process.AnalysisSequence )
#process.p = cms.Path( process.AnalysisSequence )









