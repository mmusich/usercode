import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# run the input file through the end;
# for a limited number of events, replace -1 with the desired number

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/c/castello/Zbb/PATtuple2011/DYJetsToLL/MergedOutputFile_10_1_9BS.root',
                                                                           'rfio:/castor/cern.ch/user/c/castello/Zbb/PATtuple2011/DYJetsToLL/MergedOutputFile_11_1_Cxe.root',
                                                                           'rfio:/castor/cern.ch/user/c/castello/Zbb/PATtuple2011/DYJetsToLL/MergedOutputFile_1_2_gyj.root',
                                                                           'rfio:/castor/cern.ch/user/c/castello/Zbb/PATtuple2011/DYJetsToLL/MergedOutputFile_2_1_MuE.root',
                                                                           'rfio:/castor/cern.ch/user/c/castello/Zbb/PATtuple2011/DYJetsToLL/MergedOutputFile_3_1_0rq.root',
                                                                           'rfio:/castor/cern.ch/user/c/castello/Zbb/PATtuple2011/DYJetsToLL/MergedOutputFile_4_1_bmc.root',
                                                                           'rfio:/castor/cern.ch/user/c/castello/Zbb/PATtuple2011/DYJetsToLL/MergedOutputFile_5_1_En5.root',
                                                                           'rfio:/castor/cern.ch/user/c/castello/Zbb/PATtuple2011/DYJetsToLL/MergedOutputFile_6_1_S10.root',
                                                                           'rfio:/castor/cern.ch/user/c/castello/Zbb/PATtuple2011/DYJetsToLL/MergedOutputFile_7_1_ShC.root',
                                                                           'rfio:/castor/cern.ch/user/c/castello/Zbb/PATtuple2011/DYJetsToLL/MergedOutputFile_8_1_wsO.root',
                                                                           'rfio:/castor/cern.ch/user/c/castello/Zbb/PATtuple2011/DYJetsToLL/MergedOutputFile_9_1_C1G.root'
                                                                           )
                            )

# FileService is mandatory, as the following analyzer module 
# will want it, to create output histogram file
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root")
                                   )

# Load and configure the filter
from ZbbAnalysis.AnalysisStep.zb_acceptance_filter_cfi import zb_acceptance_filter 

process.AccZplusbmmfilter = zb_acceptance_filter.clone(
    flavLept   = cms.untracked.string('Zmm'),
    etaMuMax = cms.untracked.double(2.1),
    ptMuMin  = cms.untracked.double(17.),
    etaEleMax = cms.untracked.double(2.5),
    ptEleMin  = cms.untracked.double(17.), 
    flavJet = cms.untracked.string('b'),
    etaGenJetMax = cms.untracked.double(3.5),
    ptGenJetMin  = cms.untracked.double(15),
    ngoodGenJet = cms.untracked.int32(1),
    isExclusive = cms.untracked.bool(False),  # if True -> == ngoodGenJet, if False -> >= ngoodGenJet
    jetSrc = cms.untracked.InputTag("patJets")
    )

process.AccZplusbeefilter = zb_acceptance_filter.clone(
    flavLept   = cms.untracked.string('Zee'),
    etaMuMax = cms.untracked.double(2.1),
    ptMuMin  = cms.untracked.double(17.),
    etaEleMax = cms.untracked.double(2.5),
    ptEleMin  = cms.untracked.double(17.),                                    
    flavJet = cms.untracked.string('b'),
    etaGenJetMax = cms.untracked.double(3.5),
    ptGenJetMin  = cms.untracked.double(15),
    ngoodGenJet = cms.untracked.int32(1),
    isExclusive = cms.untracked.bool(False),  # if True -> == ngoodGenJet, if False -> >= ngoodGenJet
    jetSrc = cms.untracked.InputTag("patJets")
    )

process.AccZplusbfilter = zb_acceptance_filter.clone(
    flavLept  = cms.untracked.string('all'),
    etaMuMax  = cms.untracked.double(2.1),
    ptMuMin   = cms.untracked.double(17.),
    etaEleMax = cms.untracked.double(2.5),
    ptEleMin  = cms.untracked.double(17.),                                     
    flavJet   = cms.untracked.string('b'),
    etaGenJetMax = cms.untracked.double(3.5),
    ptGenJetMin  = cms.untracked.double(15),
    ngoodGenJet  = cms.untracked.int32(1),
    isExclusive  = cms.untracked.bool(False),  # if True -> == ngoodGenJet, if False -> >= ngoodGenJet
    jetSrc = cms.untracked.InputTag("patJets")
    )

process.p1 = cms.Path(process.AccZplusbmmfilter)
process.p2 = cms.Path(process.AccZplusbeefilter)
process.p3 = cms.Path(process.AccZplusbfilter)

