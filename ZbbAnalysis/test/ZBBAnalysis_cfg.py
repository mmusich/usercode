import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("ZZAnalysis.Examples.ZJJNtuplizer_cff")

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

process.p = cms.Path(
    process.autreSeq+
   # process.muonAnalysis+                 # Add here your analyzer if any
    process.ZMMNtuplizer+
    process.ZEENtuplizer+    
    process.JetNtuplizer+
   # process.PrimaryVerticesNtuplizer+     # To be fixed! 
   # process.BeamSpotNtuplizer+            # To be fixed!
    process.METNtuplizer
    )

process.endPath = cms.EndPath(process.edmNtuplesOut)

process.schedule = cms.Schedule(
    process.p,
    process.endPath
    )



