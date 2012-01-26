import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

options = VarParsing.VarParsing()
options.register('isDA',True,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool,"True: run DA clustering; False: run on GAP clustering")
options.register('isLight',True,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool,"True: run light; False: run on full")
options.register('isSafeAPE',True,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool,"True: run with safe APE; False: run with zero APE")
options.register('maxEvents',-1,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int,"Number of events to process (-1 for all)")
options.parseArguments()

process.load("Alignment.OfflineValidation.ALCARECO_MinBias_QCD_42X_cff");
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100
                                                      
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

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

from CondCore.DBCommon.CondDBSetup_cfi import *

 ##
 ## Load Global position record
 ##
# process.globalPosition = cms.ESSource("PoolDBESSource",CondDBSetup,
#                       toGet = cms.VPSet(cms.PSet(
#                       record = cms.string('GlobalPositionRcd'),
#                       tag= cms.string('IdealGeometry')
#                        )),
#                       connect =cms.string('frontier://FrontierProd/CMS_COND_31X_FROM21X')
#                       )
# process.es_prefer_GPRcd = cms.ESPrefer("PoolDBESSource","globalPosition")


 ##
 ## Load beamspot
 ##
process.beamspot = cms.ESSource("PoolDBESSource",CondDBSetup,
                                toGet = cms.VPSet(cms.PSet( record = cms.string('BeamSpotObjectsRcd'),
                                                            tag= cms.string('Realistic7TeVCollisions2011_START311_V2_v2_mc')
                                                            )),
                                connect =cms.string('frontier://FrontierProd/CMS_COND_31X_BEAMSPOT')
                                )

process.es_prefer_beamspot = cms.ESPrefer("PoolDBESSource","beamspot")

 ##
 ## Get the BeamSpot
 ##
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")

 ##
 ## Get the GlogalTag
 ##
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START42_V12::All"  # take your favourite

 ##
 ## Get Alignment constants
 ##
process.trackerAlignment = cms.ESSource("PoolDBESSource",CondDBSetup,
                                        connect = cms.string('frontier://FrontierProd/CMS_COND_31X_ALIGNMENT'),
                                        timetype = cms.string("runnumber"),
                                        toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentRcd'),
                                                                   tag = cms.string('TrackerAlignment_2010Realistic_mc')
                                                                   ))
                                        )
process.es_prefer_trackerAlignment = cms.ESPrefer("PoolDBESSource", "trackerAlignment")


root_out_file='PVValidation_Startup_MinBias_QCD_42X_'

if options.isSafeAPE:
     
     ##
     ## Get safe APE
     ##
     process.setAPE = cms.ESSource("PoolDBESSource",CondDBSetup,
                                   connect = cms.string('sqlite_file:/afs/cern.ch/user/m/musich/public/moduleDepenentAPE/SafeAPEforFirstCollisions.db'),
                                   timetype = cms.string("runnumber"),
                                   toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentErrorRcd'),
                                                              tag = cms.string('AlignmentErrors')
                                                              ))
                                   )
     root_out_file+='_safeAPE'
     
else:
     
     ##
     ## Get zero APE
     ##
     process.setAPE = cms.ESSource("PoolDBESSource",CondDBSetup,
                                   connect = cms.string('frontier://FrontierProd/CMS_COND_31X_FROM21X'),
                                   timetype = cms.string("runnumber"),
                                   toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentErrorRcd'),
                                                              tag = cms.string('TrackerIdealGeometryErrors210_mc')
                                                              ))
                                   )
     root_out_file+='_zeroAPE'

process.es_prefer_setAPE = cms.ESPrefer("PoolDBESSource", "setAPE")


 ##
 ## Minimal Event Selection
 ##

###  based on https://twiki.cern.ch/twiki/bin/viewauth/CMS/Collisions2010Recipes#Generic_Recommendations

# process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
#                                                       ##vertexCollection = cms.InputTag('offlinePrimaryVerticesDA100um'),
#                                                       vertexCollection = cms.InputTag('offlinePrimaryVertices'),
#                                                       minimumNDOF = cms.uint32(4) ,
#  						        maxAbsZ = cms.double(24),	
#  						        maxd0 = cms.double(2)	
#                                                       )

 ##
 ## Load and Configure event selection
 ##
process.primaryVertexFilter = cms.EDFilter("VertexSelector",
                                           src = cms.InputTag("offlinePrimaryVertices"),
                                           cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
                                           filter = cms.bool(True)
                                           )

process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  src =  cms.untracked.InputTag("ALCARECOTkAlMinBias"),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )

process.goodvertexSkim = cms.Sequence(process.primaryVertexFilter + process.noscraping)

 ##
 ## Load and Configure TrackRefitter1
 ##
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
import RecoTracker.TrackProducer.TrackRefitters_cff
process.TrackRefitter1 = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone()
process.TrackRefitter1.src = "ALCARECOTkAlMinBias"
process.TrackRefitter1.TrajectoryInEvent = True
process.TrackRefitter1.TTRHBuilder = "WithTrackAngle"

if options.isDA:
     root_out_file+='_DAClustering.root'
     process.PVValidation = cms.EDAnalyzer("MultiPVValidation",
                                           TrackCollectionTag = cms.InputTag("TrackRefitter1"),
                                           OutputFileName = cms.string(root_out_file),
                                           Debug = cms.bool(False),
                                           isLightNtuple = cms.bool(options.isLight),

                                           TkFilterParameters = cms.PSet(algorithm=cms.string('filter'),
                                                                         maxNormalizedChi2 = cms.double(20.0),
                                                                         minPixelLayersWithHits=cms.int32(2),
                                                                         minSiliconLayersWithHits = cms.int32(5),
                                                                         maxD0Significance = cms.double(10000.0), 
                                                                         minPt = cms.double(0.0),
                                                                         trackQuality = cms.string("any")
                                                                         ),
                                           
                                           TkClusParameters = cms.PSet(algorithm   = cms.string("DA"),
                                                                       TkDAClusParameters = cms.PSet(coolingFactor = cms.double(0.6),  #  moderate annealing speed
                                                                                                     Tmin = cms.double(4.),            #  end of annealing
                                                                                                     vertexSize = cms.double(0.01),    #  ~ resolution / sqrt(Tmin)
                                                                                                     d0CutOff = cms.double(300.),        # downweight high IP tracks 
                                                                                                     dzCutOff = cms.double(400.)         # outlier rejection after freeze-out (T<Tmin)
                                                                                                     )
                                                                       )
                                           )
     
else:
     root_out_file+='_GAPClustering.root'
     process.PVValidation = cms.EDAnalyzer("MultiPVValidation",
                                           TrackCollectionTag = cms.InputTag("TrackRefitter1"),
                                           OutputFileName = cms.string(root_out_file),
                                           Debug = cms.bool(False),
                                           isLightNtuple = cms.bool(options.isLight),
                                           
                                           TkFilterParameters = cms.PSet(maxNormalizedChi2 = cms.double(5.0),       ## chi2ndof < 5   
                                                                         minSiliconLayersWithHits = cms.int32(7),   ## hits > 7 
                                                                         maxD0Significance = cms.double(10000.),    ## fake cut (requiring 1 PXB hit)  
                                                                         minPt = cms.double(1.0),                   ## better for softish events 
                                                                         minPixelLayersWithHits = cms.int32(1),     ## hits > 2
                                                                         trackQuality = cms.string("any")
                                                                         ),
                                           
                                           TkClusParameters = cms.PSet(algorithm   = cms.string('gap'),
                                                                       TkGapClusParameters = cms.PSet(zSeparation = cms.double(0.2)  ## 0.2 cm max separation betw. clusters
                                                                                                      ) 
                                                                       )
                                           )

process.p = cms.Path(process.goodvertexSkim*process.offlineBeamSpot*process.TrackRefitter1*process.PVValidation)

