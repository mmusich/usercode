
# Tell the process which files to use as the source
#process.source = cms.Source ("PoolSource",
#          fileNames = cms.untracked.vstring ("/store/relval/CMSSW_3_9_5/RelValTTbar/GEN-SIM-RECO/START39_V6-v1/0008/0AEEDFA4-88FA-DF11-B6FF-001A92811718.root")
#)

#readFiles = cms.untracked.vstring()

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring("ProductNotFound")    # make this exception fatal
    ,fileMode  =  cms.untracked.string('FULLMERGE') # any file order: caches all lumi/run products (memory!)
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",fileNames = myfiles.readFiles)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Tell the process what filename to use to save the output

process.Out = cms.OutputModule("PoolOutputModule")
process.Out.fileName = cms.untracked.string(myfiles.zbbPAT)

# make sure everything is hooked up
process.end = cms.EndPath(process.Out)

