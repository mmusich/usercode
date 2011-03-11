import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("ZbbAnalysis.ZJJNtuplizer_cff")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'GR_R_39X_V5::All'


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    
    ## Ernesto's file merge from Bari
    "rfio:/castor/cern.ch/cms/store/user/emiglior/test/zzSkim/test/PAT_03Mar11/Mu_ZbbSkimDec22ReReco_PAT397_03Mar11.root"
    #"rfio:/castor/cern.ch/cms/store/user/emiglior/test/zzSkim/test/PAT_03Mar11/Electron_ZbbSkimDec22ReReco_PAT397_03Mar11.root"
    #"rfio:/castor/cern.ch/cms/store/user/emiglior/test/zzSkim/test/PAT_03Mar11/EG_ZbbSkimDec22ReReco_PAT397_03Mar11.root"
    )
)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.edmNtuplesOut.fileName = cms.untracked.string('file:/tmp/musich/ZJJCandNtuples_Mu_CleanedJets.root')
#process.edmNtuplesOut.fileName = cms.untracked.string('file:/tmp/musich/ZJJCandNtuples_Ele_CleanedJets.root')
#process.edmNtuplesOut.fileName = cms.untracked.string('file:/tmp/musich/ZJJCandNtuples_Egamma_CleanedJets.root.root')

#    _      _      _____ _                            
#   (_)    | |    / ____| |                           
#    _  ___| |_  | |    | | ___  __ _ _ __   ___ _ __ 
#   | |/ _ \ __| | |    | |/ _ \/ _` | '_ \ / _ \ '__|
#   | |  __/ |_  | |____| |  __/ (_| | | | |  __/ |   
#   | |\___|\__|  \_____|_|\___|\__,_|_| |_|\___|_|   
#  _/ |                                               
# |__/                                                

process.cleanPatJets = cms.EDProducer("PATJetCleaner",
    src = cms.InputTag("cleanPatJets"),
    preselection = cms.string('pt > 2.0 && abs(eta) < 2.1'),
    checkOverlaps = cms.PSet(
        ele = cms.PSet(
           src       = cms.InputTag("preselElectronsWithTrigger"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.5),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(True),
        ),
        mu = cms.PSet(
           src       = cms.InputTag("preselMuonsWithTrigger"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.5),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(True),
        ),
    ),
    finalCut = cms.string(''),
)

process.cleanPatJetsPF = process.cleanPatJets.clone( src = cms.InputTag("cleanPatJetsPF") )
process.cleanPatJetsNoPU = process.cleanPatJets.clone( src = cms.InputTag("cleanPatJetsNoPU") )
process.cleanPatJetsJPT = process.cleanPatJets.clone( src = cms.InputTag("cleanPatJetsJPT") )

#  ___ _ _ _                         ___              _ _    _      _          
# | __(_) | |_ ___ _ _   ___ _ _    / __|__ _ _ _  __| (_)__| |__ _| |_ ___ ___
# | _|| | |  _/ -_) '_| / _ \ ' \  | (__/ _` | ' \/ _` | / _` / _` |  _/ -_|_-<
# |_| |_|_|\__\___|_|   \___/_||_|  \___\__,_|_||_\__,_|_\__,_\__,_|\__\___/__/

# process.jetFilter = cms.EDFilter("CandViewCountFilter",
#                                  src = cms.InputTag("cleanPatJets"),
#                                  minNumber = cms.uint32(2),
#                                  )

process.autreSeq = cms.Sequence( 
    (process.cleanPatJets + 
    process.cleanPatJetsPF + 
    process.cleanPatJetsNoPU + 
    process.cleanPatJetsJPT)
    # + process.jetFilter                  # eventually filter on candidate jets 
    )

process.p = cms.Path(
    process.autreSeq+
    # process.muonAnalysis+                # Add here your analyzer if any
    process.ZMMNtuplizer+
    process.ZEENtuplizer+    
    process.JetNtuplizer+
    # process.PrimaryVerticesNtuplizer+    # To be fixed! 
    # process.BeamSpotNtuplizer+           # To be fixed! 
    process.METNtuplizer
    )

process.endPath = cms.EndPath(process.edmNtuplesOut)

process.schedule = cms.Schedule(
    process.p,
    process.endPath
    )



