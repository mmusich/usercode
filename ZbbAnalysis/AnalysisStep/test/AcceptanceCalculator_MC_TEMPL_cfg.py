import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("TEST")

###################################################################
# Setup 'standard' options
###################################################################
options = VarParsing.VarParsing()
options.register('Sample',
                 "MadGraph", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string,         # string, int, or float
                 "Sample to be processed (ZIncl,MadGraph,Sherpa,aMCatNLO)")
options.register('maxEvents',
                 -1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to process (-1 for all)")
options.register('channel',
                 "Z_e",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Z decay channel? (Z_e,Z_m)")
options.parseArguments()

###################################################################
# Messages
###################################################################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))
 
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'START42_V13::All'

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


###############################################################
# Producer: produce needed collections of genJets
###############################################################
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.genParticlesForPartonJets = process.genParticlesForJets.clone()
process.genParticlesForPartonJets.partonicFinalState = True
process.genParticlesForPartonJets.excludeFromResonancePids = cms.vuint32(11, 12, 13, 14, 15, 16)

process.genParticlesForJets.ignoreParticleIDs = cms.vuint32(
    1000022, 2000012, 2000014,
    2000016, 1000039, 5000039,
    4000012, 9900012, 9900014,
    9900016, 39, 12, 14, 16
    )

process.load("RecoJets.JetProducers.ak5GenJets_cfi")
process.ak5GenJetsNoNuBSM  =  process.ak5GenJets

process.ak7GenJetsNoNuBSM  =  process.ak5GenJets.clone()
process.ak7GenJetsNoNuBSM.rParam = cms.double(0.7)

process.ak5PartonJets  =  process.ak5GenJets.clone()
process.ak5PartonJets.src = cms.InputTag("genParticlesForPartonJets")

process.ak7PartonJets  =  process.ak5GenJets.clone()
process.ak7PartonJets.rParam = cms.double(0.7)
process.ak7PartonJets.src = cms.InputTag("genParticlesForPartonJets")

process.pgenjets = cms.Sequence(
    process.genParticlesForJets
    *process.ak5GenJetsNoNuBSM
    #*process.ak7GenJetsNoNuBSM
    *process.genParticlesForPartonJets
    *process.ak5PartonJets   #exclude for the moment
    #*process.ak7PartonJets
    )

###############################################################
# Producer: create new collections of leptons from HF semileptonic decays
###############################################################
# process.LeptFromDecayProducer = cms.EDProducer("LeptFromDecayProducer",
#                                             MuonCollection = cms.InputTag("patMuonsWithTrigger")
#     )

###############################################################
#MC filter: select events with muons from b semileptonic decays 
###############################################################
# process.MCFilter = cms.EDFilter("LeptFromDecayFilter",
#                                 src = cms.InputTag("prunedGen"),
#                                 JetCollection = cms.InputTag("cleanPatJetsPF"),
#                                 ElectronCollection = cms.InputTag("patElectronsWithTrigger"),
#                                 MuonCollection = cms.InputTag("patMuonsWithTrigger"),
#                                 ZmmCollection = cms.InputTag('zMMCand'),
#                                 DecayChainSelection = cms.untracked.string('c>m')
#                                 )

####################################################
#MC Analyzer: spectra of MC Truth Analyzer
####################################################
# process.LeptFromDecayAnalyzer = cms.EDAnalyzer("LeptFromDecayAnalyzer",
#                                                src = cms.InputTag("LeptFromDecayProducer"),
#                                                JetCollection = cms.InputTag("patJets"),
#                                                ElectronCollection = cms.InputTag("patElectronsWithTrigger"),
#                                                MuonCollection = cms.InputTag("patMuonsWithTrigger"),
#                                                genParticleCollection = cms.InputTag("prunedGen"),
#                                                MuonFromBCollection = cms.InputTag("LeptFromDecayProducer:patMuonsFromB:TEST"),
#                                                MuonFromCCollection = cms.InputTag("LeptFromDecayProducer:patMuonsFromC:TEST"),
#                                                MuonFromBCCollection = cms.InputTag("LeptFromDecayProducer:patMuonsFromBC:TEST"),
#                                                ZmmCollection = cms.InputTag('zMMCand'),
#                                                HFSelection = cms.untracked.string("b")
#                                                )

# process.out = cms.OutputModule("PoolOutputModule",
#                                fileName = cms.untracked.string('/tmp/musich/ZbbToLL_PAT_NewCollections.root'),
#                                outputCommands = cms.untracked.vstring('keep *')
#                                )

###################################################################
# Jet Cleaning
###################################################################
process.cleanPatJets = cms.EDProducer("PATJetCleaner",
                                      src = cms.InputTag("patJets"),
                                      #preselection  = cms.string(''),
                                      preselection = cms.string('pt > 0.0 && abs(eta) < 5.'),
                                      checkOverlaps = cms.PSet(ele = cms.PSet(src       = cms.InputTag("offlinePrimaryVertexFromZ:patElectronsFromZ"), #Z daughters
                                                                              #src       = cms.InputTag("goldenElectrons"),
                                                                              algorithm = cms.string("byDeltaR"),
                                                                              preselection        = cms.string(""),
                                                                              deltaR              = cms.double(0.),
                                                                              checkRecoComponents = cms.bool(False),
                                                                              pairCut             = cms.string(""),
                                                                              requireNoOverlaps = cms.bool(True)
                                                                              ),
                                                               mu = cms.PSet(src       = cms.InputTag("offlinePrimaryVertexFromZ:patMuonsFromZ"), #Z daughters
                                                                             #src       = cms.InputTag("goldenMuons"),
                                                                             algorithm = cms.string("byDeltaR"),
                                                                             preselection        = cms.string(""),
                                                                             deltaR              = cms.double(0.),
                                                                             checkRecoComponents = cms.bool(False),
                                                                             pairCut             = cms.string(""),
                                                                             requireNoOverlaps = cms.bool(True)
                                                                             ),
                                                               ),
                                      finalCut = cms.string(''),
                                      )

process.JetCleaningSeq = cms.Sequence( process.cleanPatJets )

###################################################################
# HF Filter
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

###################################################################
# MC Truth Analyzer
###################################################################
# process.MCTruthAnalyzer = cms.EDAnalyzer("MCTruthAnalyzer",
#                                          JetCollection = cms.InputTag("cleanPatJets"),
#                                          genJetSrc = cms.InputTag("patJets:genJets"),
#                                          src = cms.InputTag("genParticles"),
#                                          ElectronCollection = cms.InputTag("patElectrons"),
#                                          MuonCollection = cms.InputTag("patMuons"),
#                                          ZmmCollection = cms.InputTag('zMMCand'),
#                                          isMCatNLO = cms.untracked.bool(False)
#                                          )

# process.ZLONLOHistogrammer = cms.EDAnalyzer("ZLONLOHistogrammer",
#                                             genParticles = cms.InputTag("genParticles"),
#                                             nbinsMass = cms.untracked.uint32(100),
#                                             nbinsPt = cms.untracked.uint32(100),
#                                             nbinsAng = cms.untracked.uint32(100),
#                                             massMax = cms.untracked.double(300.),
#                                             ptMax = cms.untracked.double(200.),
#                                             angMax = cms.untracked.double(3.142),
#                                             accPtMin = cms.untracked.double(20.),
#                                             accMassMin = cms.untracked.double(40.),
#                                             accMassMax = cms.untracked.double(300.),
#                                             accEtaMin = cms.untracked.double(-2.1),
#                                             accEtaMax = cms.untracked.double(2.1),
#                                             isMCatNLO = cms.untracked.bool(False)
#                                             )

#################################################################
# SFb from DB                      
##################################################################
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")
process.load("RecoBTag.PerformanceDB.BTagPerformanceDB1107")

###################################################################
# Import b-efficiency scaling factors
###################################################################
from ZbbAnalysis.AnalysisStep.effbctables_cff import mc_effbc_ssvhem_2011
from ZbbAnalysis.AnalysisStep.effbctables_cff import mc_effbc_ssvhpt_2011

###################################################################
# Acceptance Calculator
###################################################################
process.calculatePartonAcceptance = cms.EDAnalyzer("ZbbAcceptanceCalculator",
                                                   GenSrc = cms.InputTag("genParticles"),
                                                   ElectronCollection = cms.InputTag("patElectrons"),
                                                   MuonCollection = cms.InputTag("patMuons"),
                                                   JetCollection = cms.InputTag("cleanPatJets"),
                                                   # genJetSrc = cms.InputTag("patJets:genJets"),
                                                   genJetSrc  = cms.InputTag("ak5GenJetsNoNuBSM"),
                                                   ZmmCollection = cms.InputTag('zMMCand'),
                                                   ZeeCollection = cms.InputTag('zEECand'),
                                                   unLockDefaultCuts = cms.bool(False),  # bool to unlock default acceptance cuts
                                                   useStatus3forMuons = cms.bool(True), 	    
                                                   useStatus3forElectrons = cms.bool(True), 	    
                                                   doFSRCorrectionForMuons = cms.bool(False), 
                                                   doFSRCorrectionForElectrons = cms.bool(False),
                                                   jetEtaCut   = cms.double(2.1),
                                                   muonEtaCut  = cms.double(2.1),
                                                   eleEtaCut   = cms.double(2.5),
                                                   jetPtCut    = cms.double(25.),
                                                   genJetPtCut = cms.double(0.),
                                                   muonPtCut   = cms.double(20.),
                                                   elePtCut    = cms.double(25.),
                                                   minMassCut  = cms.double(60.),
                                                   maxMassCut  = cms.double(120.),
                                                   dRLeptonMatch = cms.double(0.3),
                                                   dRJetMatch    = cms.double(0.5),
                                                   BparticlePtCut= cms.double(0.),
                                                   isMCatNLO     = cms.untracked.bool(False),
                                                   applyPUCorr   = cms.untracked.bool(True),
                                                   verbose       = cms.untracked.bool(False),
                                                   useClopperPearsonErrors = cms.untracked.bool(True),
                                                   saveNTuple    = cms.untracked.bool(True),
                                                   PartonLevel   = cms.untracked.bool(False),
                                                   DecayChainSelection = cms.untracked.string("CHANNELTEMPLATE"),
                                                   OutfileName = cms.string('OUTLISTTEMPLATE.txt'),
                                                   EffbcMCSSVHE = mc_effbc_ssvhem_2011,
                                                   EffbcMCSSVHP = mc_effbc_ssvhpt_2011
                                                   )

process.calculateHadronAcceptance = process.calculatePartonAcceptance.clone()
process.calculateHadronAcceptance.genJetSrc = cms.InputTag("ak5GenJetsNoNuBSM")
process.calculateHadronAcceptance.PartonLevel = cms.untracked.bool(False)

##################################################################
# hack to force NLO weights rule
##################################################################
#if options.Sample=="aMCatNLO":        
# process.ZLONLOHistogrammer =cms.untracked.bool(True)     
# process.MCTruthAnalyzer.isMCatNLO = cms.untracked.bool(True)
# process.calculateAcceptance.isMCatNLO = cms.untracked.bool(True)   

###################################################################
# Definition of analysis sequence
###################################################################
#process.AnalysisSequence   = cms.Sequence(process.b_heavyflavorfilter*process.JetCleaningSeq)

process.AnalysisSequence   = cms.Sequence(process.JetCleaningSeq)

process.calcAcc= cms.Sequence(process.pgenjets*
                              #process.ZLONLOHistogrammer* 
                              #process.MCTruthAnalyzer*
                              process.calculatePartonAcceptance*
                              process.calculateHadronAcceptance
                              )

##############################################################################
# Patch for sum(E) violation in PYTHIA/MADGRAPH DY samples
# https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/1489.html
##############################################################################
process.load("ZbbAnalysis.Tools.TotalKinematicsFilter_cfi")

# apply kin filter to the MadGraph sample
if options.Sample=="MadGraph":
    process.AnalysisSequence+=process.totalKinematicsFilter

process.AnalysisSequence+=process.calcAcc
###################################################################
# switches
###################################################################
# if options.Sample=="ZIncl":    
#     readFiles.extend(ZjetsInclSrc)
#     process.AnalysisSequence+=process.calcAcc
#     if options.channel=="Z_m":
#         outputRootFileName=cms.string('Acceptance_ZInclToMM_MADGRAPH.root')
#     elif options.channel=="Z_e":
#         outputRootFileName=cms.string('Acceptance_ZInclToEE_MADGRAPH.root')
# elif options.Sample=="MadGraph":    
#     readFiles.extend(MadGraphSrc)
#     process.AnalysisSequence+=process.calcAcc
#     if options.channel=="Z_m":
#         outputRootFileName=cms.string('Acceptance_ZbbToMM_MADGRAPH.root')
#     elif options.channel=="Z_e":
#         outputRootFileName=cms.string('Acceptance_ZbbToEE_MADGRAPH.root')
# elif options.Sample=="Sherpa":
#     process.AnalysisSequence+=process.calcAcc
#     readFiles.extend(SherpaSrc)
#     if options.channel=="Z_m":
#         outputRootFileName=cms.string('Acceptance_ZbbToMM_SHERPA.root')
#     elif options.channel=="Z_e":
#         outputRootFileName=cms.string('Acceptance_ZbbToEE_SHERPA.root')
# elif options.Sample=="aMCatNLO":
#     process.AnalysisSequence+=process.calcAcc
#     readFiles.extend(aMCatNLOSrc)
#     if options.channel=="Z_m":
#         outputRootFileName=cms.string('Acceptance_ZbbToMM_aMCatNLO.root')
#     elif options.channel=="Z_e":
#         outputRootFileName=cms.string('Acceptance_ZbbToEE_aMCatNLO.root')
 
###################################################################
# output file
###################################################################
process.TFileService = cms.Service("TFileService", fileName = outputRootFileName)

###################################################################
# check if events in the list
###################################################################
process.checkList = cms.EDFilter("EventListFilter",
                                 Inputfile=cms.untracked.string("list.txt")
                                 )

process.p = cms.Path(process.AnalysisSequence)
#process.p = cms.Path(process.checkList*process.printTree)
#process.p = cms.Path(process.b_heavyflavorfilter)
#process.e = cms.EndPath(process.out)




