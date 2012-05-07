import FWCore.ParameterSet.Config as cms

# Give the process a name
process = cms.Process("MergeEDMFiles")

# Tell the process which files to use as the source
process.source = cms.Source ("PoolSource",
          fileNames = cms.untracked.vstring ("/store/relval/CMSSW_3_9_5/RelValTTbar/GEN-SIM-RECO/START39_V6-v1/0008/0AEEDFA4-88FA-DF11-B6FF-001A92811718.root")
)

# tell the process to only run over 100 events (-1 would mean run over
#  everything
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32 (-1)
)

# Tell the process what filename to use to save the output
process.Out = cms.OutputModule("PoolOutputModule",
         fileName = cms.untracked.string ("MergedOutputFile.root")
)

# make sure everything is hooked up
process.end = cms.EndPath(process.Out)

