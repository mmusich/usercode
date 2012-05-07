import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.patTemplate_cfg import *

#process.patDefaultSequence = cms.Sequence() ##39X

# add trigger information (trigTools)
from PhysicsTools.PatAlgos.tools.trigTools import *

#######switchOnTrigger(process,sequence='candidateSequence',hltProcess = '*') ##39X
switchOnTrigger(process,sequence='patDefaultSequence',hltProcess = '*')

#####################
#                   #
#  RUN PARAMETERS   #  
#                   #
#####################

# load the coreTools of PAT
from PhysicsTools.PatAlgos.tools.coreTools import *


### For MC only ################
isMC = True
################################

process.load('Configuration.StandardSequences.Reconstruction_cff') # to define ak5PFJets
process.load('Configuration.StandardSequences.Services_cff')
process.load("RecoLuminosity.LumiProducer.lumiProducer_cff") 

#process.GlobalTag.globaltag = cms.string('GR_R_39X_V5::All')
#process.GlobalTag.globaltag = cms.string('GR_R_311_V4::All')
process.GlobalTag.globaltag = cms.string('START42_V13::All')


process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.suppressWarning += ['patTrigger'] # 1 warning per event on old runs!

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                      SkipEvent = cms.untracked.vstring("ProductNotFound"))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
    #'rfio:/castor/cern.ch/user/m/musich/ZBBTest/Skim_Ele_ZplusJet.root'
    #'rfio:/castor/cern.ch/user/m/musich/ZBBTest/pickeventsEle2011v1_1_3_WzW.root',
    #'rfio:/castor/cern.ch/user/m/musich/ZBBTest/Skim_Ele_ZplusJet.root'
    'file:/tmp/emiglior/FEF00155-78CD-E011-B8F6-0026189437F8.root'
    #/tmp/castello/F85EB302-CC7B-E011-BC89-003048678FF8.root'
    )
)

#####################
#                   #
#      OUTPUT       #  
#                   #
#####################

process.out.fileName = cms.untracked.string('zbbPATSkim_42X.root')
process.out.outputCommands =  cms.untracked.vstring(
        'drop *',

        # lumi block info ------------
        'keep *_*_*_PAT',
        'keep edmMergeableCounter_*_*_*',
        'keep LumiSummary_lumiProducer_*_*',

        # pat candidate -------------
        'keep *_patMuons_*_*',
        'keep *_patElectrons_*_*',
        'keep *_patMuonsWithTrigger_*_*',
        'keep *_patElectronsWithTrigger_*_*',
        'keep *_selectedPatElectrons_*_*',
        'keep *_selectedPatMuons_*_*',
        'keep *_goldenElectrons_*_*',
        'keep *_goldenMuons_*_*',
        'keep *_cleanPatElectrons_*_*',
        
        # diLepton candidate -------------
        'keep *_zEECand_*_*',
        'keep *_zMMCand_*_*',
        'keep *_zllCand_*_*',
        
        # Jets ---------------------------
        
        'keep *_kt6PFJetsForIso_*_*',
        'keep *_kt6PFJets_*_*',

        'keep patJets_patJets_*_*',
        'keep patJets_patJetsAK5PFOffset_*_*',
        'keep recoPFCandidates_patJets_*_*',
        'keep recoTracks_generalTracks_*_*',

        'keep recoJetedmRefToBaseProdrecoTracksrecoTrackrecoTracksTorecoTrackedmrefhelperFindUsingAdvanceedmRefVectorsAssociationVector_*_*_*',
        'keep recoBaseTagInfosOwned_selectedPatJets_tagInfos_PAT',
        'keep recoBaseTagInfosOwned_patJets_tagInfos_PAT',
        'keep recoSecondaryVertexTagInfos_secondaryVertexTagInfosAK5PFOffset__PAT',
        'keep recoSecondaryVertexTagInfos_secondaryVertexTagInfosAOD__PAT',
        'keep recoBaseTagInfosOwned_selectedPatJetsAK5PFOffset_tagInfos_PAT',
        'keep recoBaseTagInfosOwned_patJetsAK5PFOffset_tagInfos_PAT',
  
        #'keep patJets_patJetsPF_*_*',
        #'keep patJets_patJetsNoPU_*_*',
        #'keep patJets_patJetsJPT_*_*',
 
        'keep patJets_cleanPatJets_*_*',
        'keep patJets_cleanPatJetsAK5PFOffset_*_*',
        #'keep patJets_cleanPatJetsPF_*_*',
        #'keep patJets_cleanPatJetsNoPU_*_*',
        #'keep patJets_cleanPatJetsJPT_*_*',
        
        # Tracking ---------------------
        'keep *_offlinePrimaryVertexFromZ_*_*',
        'keep *_offlinePrimaryVertices_*_*',
        'keep *_offlinePrimaryVerticesDA100um_*_*',
        'keep *_offlinePrimaryVerticesWithBS_*_*',
        'keep *_offlineBeamSpot_*_*',
        
        # MET --------------------------
        'keep *_tcMet_*_*',
        'keep *_patMETs*_*_*',
        'keep *_pfMet_*_*',
        
        # MC --------------------------
        'keep GenEventInfoProduct_generator_*_*',
        'keep *_genMetTrue_*_*',
        'keep recoGenJets_ak5GenJets_*_*',
        'keep *_addPileupInfo_*_*',
        'keep LHEEventProduct_*_*_*',
        'keep *_genParticles_*_*',
        
        # Trigger ---------------------
        'keep *_patTrigger*_*_*',
        'keep *_TriggerResults_*_*',
        'keep *_vertexMapProd_*_*',
        #'keep *_hltTriggerSummaryAOD_*_*',

        #Taus ------------------------
        'keep *_selectedPatTaus_*_*',
    )
process.out.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('zbbAllPath'))

#############################################################################
#                                                                           #
#      DETERMINISTIC ANNEALING VERTICES (100um distance; 4.2.X config)      #
#                                                                           #
#############################################################################

## process.offlinePrimaryVerticesDA100um = cms.EDProducer("PrimaryVertexProducer",
##   verbose = cms.untracked.bool(False),
##   algorithm = cms.string('AdaptiveVertexFitter'),
##   TrackLabel = cms.InputTag("generalTracks"),
##   useBeamConstraint = cms.bool(False),
##   beamSpotLabel = cms.InputTag("offlineBeamSpot"),
##   minNdof  = cms.double(0.0),
##   PVSelParameters = cms.PSet(
##       maxDistanceToBeam = cms.double(1.0)
##   ),
##   TkFilterParameters = cms.PSet(
##       algorithm=cms.string('filter'),
##       maxNormalizedChi2 = cms.double(20.0),
##       minPixelLayersWithHits=cms.int32(2),
##       minSiliconLayersWithHits = cms.int32(5),
##       maxD0Significance = cms.double(5.0),
##       minPt = cms.double(0.0),
##       trackQuality = cms.string("any")
##   ),

##   TkClusParameters = cms.PSet(
##       algorithm   = cms.string("DA"),
##       TkDAClusParameters = cms.PSet(
##           coolingFactor = cms.double(0.6),  #  moderate annealing speed
##           Tmin = cms.double(4.),            #  end of annealing
##           vertexSize = cms.double(0.01),    #  ~ resolution / sqrt(Tmin)
##           d0CutOff = cms.double(3.),        # downweight high IP tracks
##           dzCutOff = cms.double(4.)         # outlier rejection after freeze-out (T<Tmin)
##       )
##   )
## )

### To change the PV with respect to which vertexing variables are evaluated
#process.patElectrons.pvSrc = cms.InputTag("offlinePrimaryVerticesDA100um") 

### To change the PV with respect to which vertexing variables are evaluated
#process.patMuons.pvSrc = cms.InputTag("offlinePrimaryVerticesDA100um") 

###  based on https://twiki.cern.ch/twiki/bin/viewauth/CMS/Collisions2010Recipes#Generic_Recommendations

process.primaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
    filter = cms.bool(True)
    )


process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )

process.goodvertexSkim = cms.Sequence(process.primaryVertexFilter + process.noscraping)



###  based on https://twiki.cern.ch/twiki/bin/viewauth/CMS/Collisions2010Recipes#Generic_Recommendations
## process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
##                                                       vertexCollection = cms.InputTag('offlinePrimaryVerticesDA100um'),
##                                                       ##vertexCollection = cms.InputTag('offlinePrimaryVertices'),
##                                                       minimumNDOF = cms.uint32(4) ,
##  						      maxAbsZ = cms.double(24),	
##  						      maxd0 = cms.double(2)	
##                                                       )

## process.noscraping = cms.EDFilter("FilterOutScraping",
##                                   applyfilter = cms.untracked.bool(True),
##                                   debugOn = cms.untracked.bool(False),
##                                   numtrack = cms.untracked.uint32(10),
##                                   thresh = cms.untracked.double(0.25)
##                                   )

## process.goodvertexSkim = cms.Sequence(process.primaryVertexFilter + process.noscraping)

#####################
#                   #
# PU ISO CORRECTION #
#                   #
#####################

process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJetsForIso = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIso.Rho_EtaMax = cms.double(2.5)

#Trigger Matching
#process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
#process.patTrigger.onlyStandAlone = True
#process.patTrigger.processName  = '*' 
#process.patTriggerEvent.processName = '*'


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
    hDaughterVeto = cms.bool(False),
    zDaughterVeto = cms.bool(False),
    ptcut=cms.double(0)
    )

##################
#                #
# ELECTRON PATH  #
#                #
##################


ELE_PRESEL_CUT=( 
    "pt > 5.0 && abs(eta) < 2.5 &&"  +
    #"(isEE || isEB) && !isEBEEGap &&" +
    "electronID('eidVBTFRel85') == 7" 
    )

ELE_TIGHT_CUT=(
    "pt > 25.0 && abs(eta) < 2.5 &&"                               +
    #"(isEE || isEB) && !isEBEEGap &&"                              +
    "electronID('eidVBTFRel85') == 7" 
    )

process.load("PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi")
#process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")


##### ELE ID

###### Code for the liklihood eID
process.load("RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi")
process.liklihoodID = process.eidLikelihoodExt.clone() 


###### Code for VTBF electronID
import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi as vbtfid
process.eidVBTFRel95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95relIso' )
process.eidVBTFRel80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80relIso' )
process.eidVBTFRel85 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '85relIso' )## added for final analyzer
process.eidVBTFCom95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95cIso'   )
process.eidVBTFCom80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80cIso'   )


process.patElectrons.electronIDSources = cms.PSet(
    eidVBTFRel95 = cms.InputTag("eidVBTFRel95"),
    eidVBTFRel80 = cms.InputTag("eidVBTFRel80"),
    eidVBTFRel85 = cms.InputTag("eidVBTFRel85"),
    eidVBTFCom95 = cms.InputTag("eidVBTFCom95"),
    eidVBTFCom80 = cms.InputTag("eidVBTFCom80")
)


###### Electron sequences
process.vbtfIDSequence = cms.Sequence(
    process.eidVBTFRel95 +
    process.eidVBTFRel80 +
    process.eidVBTFRel85 +
    process.eidVBTFCom95 +
    process.eidVBTFCom80 
)

process.convValueMapProd = cms.EDProducer('ConvValueMapProd',
    gsfLabel = cms.untracked.InputTag("gsfElectrons"),
    tkLabel = cms.untracked.InputTag("generalTracks")
)

###### EleTriggerMatchHLT ########################################

process.eleTriggerMatchHLT = cms.EDProducer( "PATTriggerMatcherDRDPtLessByR",
    src     = cms.InputTag( "patElectrons" ),
    matched = cms.InputTag( "patTrigger" ),
    andOr          = cms.bool( False ),
    filterIdsEnum  = cms.vstring( '*' ),
    filterIds      = cms.vint32( 0 ),
    filterLabels   = cms.vstring( '*' ),
    pathNames      = cms.vstring( '*' ),
    #matchedCuts = cms.string('coll("hltL1IsoRecoEcalCandidate")||coll("hltL1NonIsoRecoEcalCandidate")'), ##39X
    matchedCuts = cms.string('path("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*")'),                                       
    maxDPtRel = cms.double( 0.5 ),
    maxDeltaR = cms.double( 0.5 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( True )
)

### patElectronsWithTrigger ###########################################
process.patElectronsWithTrigger = cms.EDProducer( "PATTriggerMatchElectronEmbedder",
    src     = cms.InputTag("patElectrons"),
    matches = cms.VInputTag(cms.InputTag('eleTriggerMatchHLT'))
)

### The electron selector (replace the default one in the patDefaultSequence) ######
process.selectedPatElectrons.src = ("patElectronsWithTrigger")
process.selectedPatElectrons.cut = (ELE_PRESEL_CUT)

### goldenElectrons (for cleaning) ######################################################
process.goldenElectrons = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string(ELE_TIGHT_CUT),
)

process.patElectrons.embedPFCandidate = False
process.patElectrons.embedSuperCluster = True
process.patElectrons.embedTrack = True
process.patElectrons.addElectronID = True
process.electronMatch.matched = "genParticles"
process.patElectrons.userData.userFloats.src = cms.VInputTag(
    cms.InputTag("convValueMapProd","dist"),
    cms.InputTag("convValueMapProd","dcot"),
    cms.InputTag("liklihoodID")
)

if isMC: 
    process.preElectronSequence = cms.Sequence(
        process.convValueMapProd +
        process.vbtfIDSequence +  
        process.liklihoodID +
        process.patTrigger + 
        process.electronMatch 
    )
else:

    #removeMCMatching(process, ['Electrons'])  # already removed at the beginning 
    process.preElectronSequence = cms.Sequence(
        process.convValueMapProd +
        process.vbtfIDSequence +  
        process.liklihoodID +
        process.patTrigger 
    )

##from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cff import *
process.patDefaultSequence.replace(
   process.patElectrons,
   process.patElectrons*(process.eleTriggerMatchHLT)*process.patElectronsWithTrigger*process.goldenElectrons
   )

##############
#            #
# MUON PATH  #
#            #
##############

MUON_PRESEL_CUT=( 
    "pt > 5. && isGlobalMuon && isTrackerMuon && globalTrack().normalizedChi2 < 15 &&" +
    "innerTrack().hitPattern().numberOfValidTrackerHits > 10 && "                      +
    "innerTrack().hitPattern().numberOfValidPixelHits > 0 && "                         +
    "globalTrack().hitPattern().numberOfValidMuonHits > 0 && "                         +
    "abs(dB) < 0.2 && "                                                                     +
    "trackIso + caloIso < 0.15 * pt && "                                               +
    "numberOfMatches > 1 && abs(eta) < 2.1" 
    )

MUON_TIGHT_CUT=(
    "pt > 20. && isGlobalMuon && isTrackerMuon && globalTrack().normalizedChi2 < 10 &&" +
    "innerTrack().hitPattern().numberOfValidTrackerHits > 10 && "                      +
    "innerTrack().hitPattern().numberOfValidPixelHits > 0 && "                         +
    "globalTrack().hitPattern().numberOfValidMuonHits > 0 && "                         +
    "abs(dB) < 0.2 && "                                                                     +
    "trackIso + caloIso < 0.15 * pt && "                                               +
    "numberOfMatches > 1 && abs(eta) < 2.1" 
    )


### muonTriggerMatchHLT ################################################################

process.muonTriggerMatchHLT = cms.EDProducer( 'PATTriggerMatcherDRDPtLessByR',
    src     = cms.InputTag( 'patMuons' ),
    matched = cms.InputTag( 'patTrigger' ),
    andOr          = cms.bool( False ),
    filterIdsEnum  = cms.vstring( '*' ),
    filterIds      = cms.vint32( 0 ),
    filterLabels   = cms.vstring( '*' ),
    pathNames      = cms.vstring( '*' ),
    matchedCuts = cms.string('path("HLT_DoubleMu7")'),
    maxDPtRel = cms.double( 0.5 ),
    maxDeltaR = cms.double( 0.1 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( False )
)

### patMuonsWithTrigger ############################################################

process.patMuonsWithTrigger = cms.EDProducer( 'PATTriggerMatchMuonEmbedder',
    src     = cms.InputTag(  'patMuons' ),
    matches = cms.VInputTag('muonTriggerMatchHLT')
)

process.patMuons.embedPFCandidate = False
process.patMuons.embedTrack = True
process.muonMatch.matched = "genParticles"

if isMC: 
##     if False: ## Turn this on to get extra info on muon MC origin, on GEN-SIM-RECO
##         process.load("MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi")
##         from MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi import addUserData as addClassByHits
##         addClassByHits(process.patMuons, labels=['classByHitsGlb'], extraInfo = True)
##         process.muonClassificationByHits = cms.Sequence(process.mix * process.trackingParticlesNoSimHits * process.classByHitsGlb )
        
##         process.preMuonSequence = cms.Sequence(
##             process.prunedGen *
##             process.muonMatch +
##             process.patTrigger +
##             process.muonClassificationByHits
##             )
##         process.MessageLogger.suppressWarning += ['classByHitsGlb'] # kill stupid RPC hit associator warning
##     else:
        
     process.preMuonSequence = cms.Sequence(
            process.muonMatch +
            process.patTrigger
            )
else:
    process.preMuonSequence = cms.Sequence(
        process.patTrigger
        )
    
### golden muons #################################
process.goldenMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(MUON_TIGHT_CUT)
    )

#from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cff import *
process.patDefaultSequence.replace(
    process.patMuons,
    process.patMuons *process.muonTriggerMatchHLT *process.patMuonsWithTrigger*process.goldenMuons
    )


### The muon selector (replace the default one in the patDefaultSequence)
process.selectedPatMuons.src = ("patMuonsWithTrigger")
process.selectedPatMuons.cut = (MUON_PRESEL_CUT)

######################
#                    #
# ELECTRON CLEANING  #
#                    # 
######################

#### Cleaner for muons reconstructed as electron (electron is discared) #######
process.cleanPatElectrons = cms.EDProducer("PATElectronCleaner",
    src = cms.InputTag("selectedPatElectrons"), 
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(
        ele = cms.PSet(
           src       = cms.InputTag("selectedPatMuons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.01),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(True),
          )
        ),
    finalCut = cms.string(''),
)

######################
#                    #
#  PRE-PAT SEQUENCE  #
#                    # 
######################

process.prePATSequence= cms.Sequence(
    process.preMuonSequence +
    process.preElectronSequence
    )

########################
#                      #
# 2l CANDIDATES        # 
#                      #
########################

## process.allLeps = cms.EDProducer("CandViewMerger",
##     src = cms.VInputTag(cms.InputTag("preselMuons"), cms.InputTag("preselElectrons"))
## )

## process.zEECand = cms.EDProducer("CandViewShallowCloneCombiner",
##     decay = cms.string('preselElectronsWithTrigger@+ preselElectronsWithTrigger@-'),
##     cut = cms.string('mass > 8 && daughter(0).pt > 0 && daughter(1).pt > 0'),
##     checkCharge = cms.bool(True)
## )

process.zEECand = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('cleanPatElectrons@+ cleanPatElectrons@-'),
    cut = cms.string('mass > 40 && daughter(0).pt > 0 && daughter(1).pt > 0'),
    checkCharge = cms.bool(True)
)

process.zMMCand = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('selectedPatMuons@+ selectedPatMuons@-'),
    cut = cms.string('mass > 40 && daughter(0).pt > 0 && daughter(1).pt > 0'),
    checkCharge = cms.bool(True)
)

process.zllCand = cms.EDProducer("CandViewMerger",
     src = cms.VInputTag(cms.InputTag("zMMCand"), cms.InputTag("zEECand"))
      )

## process.diLeps = cms.EDProducer("CandViewShallowCloneCombiner",
##     decay = cms.string('allLeps allLeps'),
##     cut = cms.string('mass > 20 && daughter(0).pt > 0 && daughter(1).pt > 0'),
##     checkCharge = cms.bool(False)
## )


process.DiLepFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("zllCand"), minNumber = cms.uint32(0))

process.allCandSequence = cms.Sequence(

        process.zEECand  +
        process.zMMCand  +
        process.zllCand  +
        process.DiLepFilter
)

################
#              #
#  Z Vertex    # 
#              #
################

from ZbbAnalysis.Tools.zvertexproducer_cfi import zvertexproducer 
process.offlinePrimaryVertexFromZ =  zvertexproducer.clone(
       VertexSrc   = cms.untracked.InputTag("offlinePrimaryVertices"),
       ZmmSrc      = cms.untracked.InputTag("zMMCand"),
       ZeeSrc      = cms.untracked.InputTag("zEECand")
       )

################
#              #
#  JET PATH    # 
#              #
################

#------------------------------------------------------------------------------------------
# Jet energy corrections to use (3_9_X)
#inputJetCorrLabel = ('AK5PF', ['L1Offset', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
# -----------------------------------------------------------------------------------------
# corrections for 4_1_X (pileup)
#inputJetCorrLabel = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']) 
# -----------------------------------------------------------------------------------------
# no residual correction in 4_2_X ( https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCor2011 )
inputJetCorrLabel = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute'])


# Add PF jets
from PhysicsTools.PatAlgos.tools.jetTools import *
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')

process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax = cms.double(5.0)
process.kt6PFJets.Ghost_EtaMax = cms.double(5.0)

##process.ak5PFJets.doRhoFastjet = True ## added
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Rho_EtaMax = cms.double(2.5)

process.patJetCorrFactors.rho = cms.InputTag('kt6PFJets','rho')

## see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#Jet_Tools

addJetCollection(process,cms.InputTag('ak5PFJets'),
               'AK5','PFOffset',
               doJTA        = True,
               doBTagging   = True,
               jetCorrLabel = ('AK5PF', cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute'])),
               doType1MET   = False,
               doL1Cleaning   = False,
               doL1Counters   = True,
               genJetCollection=cms.InputTag("ak5GenJets"),
               doJetID      = True,
               jetIdLabel   = "ak5"
               )

### the default collection become this one, i.e. patJets is replaced by this values

switchJetCollection(process,cms.InputTag('ak5PFJets'),
               doJTA        = True,
               doBTagging   = True,
               jetCorrLabel = inputJetCorrLabel,
               doType1MET   = True,
               genJetCollection=cms.InputTag("ak5GenJets"),
               doJetID      = True
               )

process.patJets.addGenJetMatch = True
process.patJets.addGenPartonMatch = True
process.patJets.getJetMCFlavour = True
process.patJets.embedPFCandidates = True
process.patJets.addAssociatedTracks = True
process.patJets.addTagInfos = True
process.patJets.tagInfoSources  = cms.VInputTag(
  cms.InputTag("secondaryVertexTagInfosAOD"),
  )

############# original part (BEGIN) #############################

## from PhysicsTools.PatAlgos.tools.jetTools import *

## process.load("PhysicsTools.PFCandProducer.pfNoPileUp_cff")  
## process.pfPileUp.Enable = True    
## process.ak5PFJetsNoPU = process.ak5PFJets.clone(    
##     src =   cms.InputTag("pfNoPileUp")  
## )
## process.jetSequence = cms.Sequence(   
##     process.pfNoPileUpSequence  +   
##     process.ak5PFJetsNoPU   
## )  

## if not isMC:
##     removeMCMatching(process, ['Jets'])

## addJetCollection(
##     process,
##     cms.InputTag("ak5PFJetsNoPU"),
##     algoLabel    = "NoPU",
##     typeLabel    = "",
##     doJTA        = True,
##     doBTagging   = True,
##     jetCorrLabel = ('AK5PF',cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute','L2L3Residual'])), ### remove L2L3Residuals for MC
##     doL1Cleaning = False,
##     doL1Counters = True,                 
##     doType1MET   = True,
##     genJetCollection=cms.InputTag("ak5GenJets"),
##     doJetID      = True,
##     jetIdLabel   = "ak5"
## )

## addJetCollection(
##     process,
##     cms.InputTag("ak5PFJets"),
##     algoLabel    = "PF",
##     typeLabel    = "",
##     doJTA        = True,
##     doBTagging   = True,
##     jetCorrLabel = ('AK5PF',cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute','L2L3Residual'])),  ### remove L2L3Residuals for MC
##     doL1Cleaning = False,
##     doL1Counters = True,                 
##     doType1MET   = True,
##     genJetCollection=cms.InputTag("ak5GenJets"),
##     doJetID      = True,
##     jetIdLabel   = "ak5"
## )

## ## switchJetCollection(
## ##     process,
## ##     cms.InputTag('ak5PFJets'),
## ##     doJTA        = True,
## ##     doBTagging   = True,
## ##     jetCorrLabel = ('AK5', 'PF'),
## ##     doType1MET   = True,
## ##     genJetCollection=cms.InputTag("ak5GenJets"),
## ##     doJetID      = True
## ## )

## process.patJets.addGenJetMatch = False
## process.patJetsPF.addGenJetMatch = False
## process.patJetsNoPU.addGenJetMatch = False

## if not isMC:
##     process.patJetsPF.addGenJetMatch = False
##     process.patJetsPF.addGenPartonMatch = False
##     process.patJetsPF.getJetMCFlavour = False
##     process.patJetsNoPU.addGenJetMatch = False
##     process.patJetsNoPU.addGenPartonMatch = False
##     process.patJetsNoPU.getJetMCFlavour = False


## process.patJets.addTagInfos = False
## process.patJetsPF.addTagInfos = False
## process.patJetsNoPU.addTagInfos = False

## process.patJets.embedPFCandidates = False
## process.patJetsPF.embedPFCandidates = False
## process.patJetsNoPU.embedPFCandidates = False

## process.patJets.addAssociatedTracks = False
## process.patJetsPF.addAssociatedTracks = False
## process.patJetsNoPU.addAssociatedTracks = False

## process.patJetCorrFactors.levels = cms.vstring(['L2Relative', 'L3Absolute','L2L3Residual']) ## remove L2L3Residual for MC
## process.patJetCorrFactorsPF.levels = cms.vstring(['L2Relative', 'L3Absolute','L2L3Residual']) 
## process.patJetCorrFactorsNoPU.levels = cms.vstring(['L2Relative', 'L3Absolute','L2L3Residual'])

############# original part (END) #############################

process.cleanPatJets = cms.EDProducer("PATJetCleaner",
    src = cms.InputTag("patJets"),
    preselection = cms.string('pt > 20.0 && abs(eta) < 2.4'),
    checkOverlaps = cms.PSet(
        ele = cms.PSet(
           src       = cms.InputTag("goldenElectrons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.5),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(False)
        ),
        mu = cms.PSet(
           src       = cms.InputTag("goldenMuons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.5),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(False)
        ),
    ),
    finalCut = cms.string(''),
)

process.cleanPatJetsAK5PFOffset = cms.EDProducer("PATJetCleaner",
                         src = cms.InputTag("patJetsAK5PFOffset"),
                         preselection = cms.string('pt > 20.0 && abs(eta) < 2.4'),
                         checkOverlaps = cms.PSet(
                            muons = cms.PSet(
                                        src = cms.InputTag("selectedPatMuons"),
                                        algorithm = cms.string("byDeltaR"),
                                        preselection = cms.string(""),
                                        deltaR = cms.double(0.5),
                                        checkRecoComponents = cms.bool(False),
                                        pairCut = cms.string(""),
                                        requireNoOverlaps = cms.bool(False),
                            ),
                            electrons = cms.PSet(
                                            src = cms.InputTag("selectedPatElectrons"),
                                            algorithm = cms.string("byDeltaR"),
                                            preselection = cms.string(""),
                                            deltaR = cms.double(0.5),
                                            checkRecoComponents = cms.bool(False),
                                            pairCut = cms.string(""),
                                            requireNoOverlaps = cms.bool(False),
                           )
                         ),
                         finalCut = cms.string('')
)

#process.cleanPatJetsPF = process.cleanPatJets.clone( src = cms.InputTag("patJetsPF") )
#process.cleanPatJetsNoPU = process.cleanPatJets.clone( src = cms.InputTag("patJetsNoPU") )

## process.jetFilter = cms.EDFilter("CandViewCountFilter",
##     src = cms.InputTag("cleanPatJetsPF"),
##     minNumber = cms.uint32(1),
## )

process.patMETs.metSource = cms.InputTag("pfMet")
process.patMETs.addGenMET = cms.bool(False)

# add MET (for PF objects)
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')

#process.jetSequence *= process.makePatJets
#process.jetSequence *=(process.cleanPatJets + process.cleanPatJetsPF + process.cleanPatJetsNoPU)
#process.jetSequence += process.patMETs

#############
#           #
# TAU PATH  # 
#           #
#############

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(
   process,
   pfTauLabelOld = cms.InputTag('shrinkingConePFTauProducer'),
   pfTauLabelNew = cms.InputTag('hpsPFTauProducer'),
   postfix=""
)

process.selectedPatTaus.cut = (
   "pt > 15 && " +
   #"tauID('leadingTrackFinding') > 0.5 &&
   "tauID('byLooseIsolation') > 0.5"
   )

#process.tauSequence = cms.Sequence(process.makePatTaus * process.selectedPatTaus)

############################################################################

process.TotalEventCounter = cms.EDProducer("EventCountProducer")
process.AfterPVFilterCounter = cms.EDProducer("EventCountProducer")
process.AfterNSFilterCounter = cms.EDProducer("EventCountProducer")
process.AfterPATCounter = cms.EDProducer("EventCountProducer")
process.AfterCandidatesCounter = cms.EDProducer("EventCountProducer")
process.AfterJetsCounter = cms.EDProducer("EventCountProducer")

#### Analyzer ################################################################
from ZbbAnalysis.AnalysisStep.PATAnalyzer_cfi import analyzePAT
process.analyzeBasicPat =  analyzePAT.clone(
    jetSrc     = cms.untracked.InputTag("patJets"),
    ##for 2010+2011
    firstRunNumber = 140000,
    lastRunNumber  = 170000
    ) 

#### B-jet analyzer (NOT YET COMMITTED) ###########################################################

## from ZbbAnalysis.AnalysisStep.PatBJetTagAnalyzer_cfi import analyzeBJets
## process.bjetanalysis =  analyzeBJets.clone(
##     primaryVertices = cms.InputTag("offlinePrimaryVertexFromZ:ZVertexProducerZVertex:PAT"),
##     tracks = cms.InputTag("generalTracks"),
##     jets = cms.InputTag("patJets"),
##     )

process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string('analyzePatBasics.root')
                                   )

#### Re-btag ###################################################################

process.load("RecoBTag.Configuration.RecoBTag_cff")
process.impactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertexFromZ:ZVertexProducerZVertex:PAT")
process.ak5JetTracksAssociatorAtVertex.jets = cms.InputTag("patJets") 

## process.candidateSequence+=( process.elSequence + process.muSequence + process.tauSequence )
## process.candidateSequence*=process.cleanPatElectrons
## process.candidateSequence+=process.AfterPATCounter
## process.candidateSequence+=process.jetSequence 

## process.readAK5PF = cms.EDAnalyzer('JetCorrectorDBReader', 
##          payloadName    = cms.untracked.string('AK5PF'),
##          printScreen    = cms.untracked.bool(False),
##          createTextFile = cms.untracked.bool(True),
##         globalTag      = cms.untracked.string('GR_R_39X_V5')
## ##         globalTag      = cms.untracked.string('GR_R_38X_V1')
##    )

##################
#                #
#   FINAL PATH   # 
#                #
##################

process.zbbAllPath = cms.Path(
    ##  process.lumiProducer +
    process.TotalEventCounter +
    ## process.b_heavyflavorfilter *
    ##process.offlinePrimaryVertices*
    process.primaryVertexFilter *
    process.AfterPVFilterCounter +
    process.noscraping *
    process.AfterNSFilterCounter +
    process.kt6PFJetsForIso +
    process.recoPFJets *
    process.prePATSequence *
    process.patDefaultSequence*
    process.AfterPATCounter*
    process.cleanPatJets *
    process.cleanPatJetsAK5PFOffset *
    ### Selection on the basis of candidates/jets
    process.allCandSequence +
    process.offlinePrimaryVertexFromZ +
    process.AfterCandidatesCounter +
    ## process.jetFilter +
    process.impactParameterTagInfos *
    process.ak5JetTracksAssociatorAtVertex *
    process.btagging *
    process.AfterJetsCounter +
    process.analyzeBasicPat
    ##*process.bjetanalysis ## NOT yet committed
    ## +process.readAK5PF
    )

# uncomment next three lines to dump the full python configuration 
#outFile = open("/tmp/musich/tmpConfig.py","w")
#outFile.write(process.dumpPython())
#outFile.close()


