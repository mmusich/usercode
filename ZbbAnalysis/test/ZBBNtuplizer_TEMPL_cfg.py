import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'GLOBALTAGTEMPLATE'  

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    "INPUTFILETEMPLATE"
    )
                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.TFileService=cms.Service('TFileService',
                                 fileName=cms.string('PLOTOUTFILETEMPLATE')
                                 )


##################
# TRIGGER FILTER #
##################
### Trigger path FIRED in the event
# import HLTrigger.HLTfilters.hltHighLevel_cfi 
# process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone() 
# process.hltFilter.HLTPaths = [ "HLT_DoubleMu3_v*" ]
# process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
# process.hltFilter.throw  = cms.bool(False)
 
###################################################################
#First Analyzer: general spectra for all PATMuons + filter studies
###################################################################
process.MuAllSpectra = cms.EDAnalyzer("MuAnalyzer",
                                      MuonCollection = cms.untracked.string('patMuonsWithTrigger'),
                                      JetCollectionPF = cms.untracked.string('cleanPatJetsPF'),
                                      ZllCollection = cms.untracked.string('zMMCand'),
                                      numEventsNames = cms.untracked.vstring('TotalEventCounter','AfterPVFilterCounter', 'AfterNSFilterCounter', 'AfterPATCounter', 'AfterCandidatesCounter', 'AfterJetsCounter')
                                      )


###################################################################
#Second Analyzer: spectra about preselection analysis
###################################################################
process.MuPreselSpectra = cms.EDAnalyzer("PreselAnalyzer",
                                      MuonCollection = cms.untracked.string('patMuonsWithTrigger'),
                                      ZllCollection = cms.untracked.string('zMMCand'),
                                                                           )

###################################################################
#Third Analyzer: Jet collections studies
###################################################################
process.JetLepton = cms.EDAnalyzer("JetLeptonCleanAnalyzer",
                                   MuonCollection = cms.untracked.string('patMuonsWithTrigger'),
                                   ZllCollection = cms.untracked.string('zMMCand'),
                                   JetCollection = cms.untracked.string('cleanPatJets'),
                                   JetCollectionPU = cms.untracked.string('cleanPatJetsNoPU'),
                                   JetCollectionJPT = cms.untracked.string('cleanPatJetsJPT'),
                                   JetCollectionPF = cms.untracked.string('cleanPatJetsPF'),
                                   )

###################################################################
#Ntuplizer
###################################################################
process.load("ZZAnalysis.Examples.ZJJNtuplizer_cff")
process.edmNtuplesOut.fileName = cms.untracked.string('OUTPUTFILETEMPLATE')

process.p = cms.Path(
    process.MuAllSpectra+
    process.MuPreselSpectra+
    process.JetLepton+
    process.ZMMNtuplizer+
    process.ZEENtuplizer+    
    process.JetNtuplizer+
    process.METNtuplizer
    )

process.endPath = cms.EndPath(process.edmNtuplesOut)

process.schedule = cms.Schedule(
    process.p,
    process.endPath
    )



