import FWCore.ParameterSet.Config as cms

def customiseCrossingFrame(process):
    process.mix.mixObjects.mixVertices.makeCrossingFrame = cms.untracked.bool(True)
    process.mix.mixObjects.mixTracks.makeCrossingFrame = cms.untracked.bool(True)

    process.theMixObjects.mixVertices.makeCrossingFrame = cms.untracked.bool(True)
    process.theMixObjects.mixTracks.makeCrossingFrame = cms.untracked.bool(True)

    process.mix.mixObjects.mixSH.crossingFrames = cms.untracked.vstring('BSCHits', 
                                                                        'FP420SI', 
                                                                        'MuonCSCHits', 
                                                                        'MuonDTHits', 
                                                                        'MuonRPCHits', 
                                                                        'TotemHitsRP', 
                                                                        'TotemHitsT1', 
                                                                        'TotemHitsT2Gem', 
                                                                        'TrackerHitsPixelBarrelHighTof', 
                                                                        'TrackerHitsPixelBarrelLowTof', 
                                                                        'TrackerHitsPixelEndcapHighTof', 
                                                                        'TrackerHitsPixelEndcapLowTof', 
                                                                        'TrackerHitsTECHighTof', 
                                                                        'TrackerHitsTECLowTof', 
                                                                        'TrackerHitsTIBHighTof', 
                                                                        'TrackerHitsTIBLowTof', 
                                                                        'TrackerHitsTIDHighTof', 
                                                                        'TrackerHitsTIDLowTof', 
                                                                        'TrackerHitsTOBHighTof', 
                                                                        'TrackerHitsTOBLowTof')
    
