import FWCore.ParameterSet.Config as cms



lumiblockanalysis = cms.EDAnalyzer("LumiBlockAnalyzer",
                                   runRangeList   = cms.VPSet(cms.PSet(runrange=cms.untracked.vuint32(0,999999))),
                                   numEventsNames = cms.untracked.vstring('TotalEventCounter','AfterPVFilterCounter', 'AfterNSFilterCounter', 'AfterPATCounter', 'AfterCandidatesCounter', 'AfterJetsCounter'),
                                   LumiSummaryTag = cms.untracked.string('lumiProducer::RECO')
                                   )



