import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'GR_R_39X_V5::All'

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    
    ## Ernesto's from Bari
    "rfio:/castor/cern.ch/cms/store/user/emiglior/test/zzSkim/test/PAT_03Mar11/Mu_ZbbSkimDec22ReReco_PAT397_03Mar11.root",
    "rfio:/castor/cern.ch/cms/store/user/emiglior/test/zzSkim/test/PAT_03Mar11/Electron_ZbbSkimDec22ReReco_PAT397_03Mar11.root",
    "rfio:/castor/cern.ch/cms/store/user/emiglior/test/zzSkim/test/PAT_03Mar11/EG_ZbbSkimDec22ReReco_PAT397_03Mar11.root"
    )
)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


process.analyzeBasicPat = cms.EDAnalyzer("PatBasicAnalyzer",
  #photonSrc   = cms.untracked.InputTag("cleanPatPhotons"),
  electronSrc = cms.untracked.InputTag("patElectronsWithTrigger"),
  muonSrc     = cms.untracked.InputTag("patMuonsWithTrigger"),                                             
  #tauSrc      = cms.untracked.InputTag(""),
  jetSrc      = cms.untracked.InputTag("cleanPatJets"),
  #corrLevel   = cms.string('Uncorrected'),
  corrLevels = cms.vstring('Uncorrected','L2Relative', 'L3Absolute', 'L2L3Residual')                                       
  ##metSrc      = cms.untracked.InputTag("pfMet")
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('analyzePatBasics.root')
)

process.p = cms.Path(process.analyzeBasicPat)




