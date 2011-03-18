import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("Alignment.OfflineValidation.DATASETTEMPLATE");
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

 ##
 ## Get the Magnetic Field
 ##
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

 ##
 ## Get the Geometry
 ##
from Geometry.CommonDetUnit.globalTrackingGeometry_cfi import *
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")

####### Begin snippet for choosing the GlobalPositionRcd
## load the Global Position Rcd
from CondCore.DBCommon.CondDBSetup_cfi import *
process.globalPosition = cms.ESSource("PoolDBESSource",CondDBSetup,
                      toGet = cms.VPSet(cms.PSet(
                      record = cms.string('GlobalPositionRcd'),
                      tag= cms.string('IdealGeometry')
                       )),
                      connect =cms.string('frontier://FrontierProd/CMS_COND_31X_FROM21X')
                      )
process.es_prefer_GPRcd = cms.ESPrefer("PoolDBESSource","globalPosition")
####### End snippet

 ##
 ## Get the BeamSpot
 ##
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")

 ##
 ## Get the GlogalTag
 ##
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GLOBALTAGTEMPLATE"  # take your favourite

 ##
 ## Get Alignment constants
 ##
from CondCore.DBCommon.CondDBSetup_cfi import *
process.trackerAlignment = cms.ESSource("PoolDBESSource",CondDBSetup,
                                        connect = cms.string('ALIGNOBJTEMPLATE'),
                                        timetype = cms.string("runnumber"),
                                        toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentRcd'),
                                                                   tag = cms.string('GEOMTAGTEMPLATE')
                                                                   ))
                                        )
process.es_prefer_trackerAlignment = cms.ESPrefer("PoolDBESSource", "trackerAlignment")

 ##
 ## Get APE
 ##
from CondCore.DBCommon.CondDBSetup_cfi import *
process.setAPE = cms.ESSource("PoolDBESSource",CondDBSetup,
                                        connect = cms.string('APEOBJTEMPLATE'),
                                        timetype = cms.string("runnumber"),
                                        toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentErrorRcd'),
                                                                   tag = cms.string('ERRORTAGTEMPLATE')
                                                                   ))
                                        )
process.es_prefer_setAPE = cms.ESPrefer("PoolDBESSource", "setAPE")

 ##
 ## Minimal Event Selection
 ##
process.oneGoodVertexFilter = cms.EDFilter("VertexSelector",
                                           # src = cms.InputTag("hiSelectedVertex"),      # for HI
                                           src = cms.InputTag("offlinePrimaryVertices"),  # for pp
                                           cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"), # tracksSize() > 3 for the older cut
                                           filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
                                           ) 

process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

 ##
 ## Load and Configure TrackRefitter1
 ##
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
import RecoTracker.TrackProducer.TrackRefitters_cff
process.TrackRefitter1 = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone()
process.TrackRefitter1.src = "TRACKTYPETEMPLATE"
process.TrackRefitter1.TrajectoryInEvent = True
process.TrackRefitter1.TTRHBuilder = "WithTrackAngle"

process.PVValidation = cms.EDAnalyzer("VALIDATIONMODULETEMPLATE",
                                      TrackCollectionTag = cms.InputTag("TrackRefitter1"),
                                      OutputFileName = cms.string("OUTFILETEMPLATE"),
                                      Debug = cms.bool(False),
                                      TkFilterParameters = cms.PSet(maxNormalizedChi2 = cms.double(5.0),       ## chi2ndof < 5   
                                                                    minSiliconLayersWithHits = cms.int32(7),   ## hits > 7 
                                                                    maxD0Significance = cms.double(1000000.0), ## fake cut (requiring 1 PXB hit)  
                                                                    minPt = cms.double(1.0),                   ## better for softish events 
                                                                    minPixelLayersWithHits = cms.int32(2),     ## hits > 2
                                                                    trackQuality = cms.string("any")
                                                                    ),
                                      TkClusParameters=cms.PSet(algorithm=cms.string('DA'),
                                                                TkDAClusParameters = cms.PSet(coolingFactor = cms.double(0.8),  #  slow annealing
                                                                                              Tmin = cms.double(4.0),           #  freezeout temperature
                                                                                              vertexSize = cms.double(0.05)     #  ~ resolution / sqrt(Tmin)
                                                                                              )
                                                                )
                                      )

process.p = cms.Path(process.oneGoodVertexFilter*process.hltPhysicsDeclared*process.offlineBeamSpot*process.TrackRefitter1*process.PVValidation)
