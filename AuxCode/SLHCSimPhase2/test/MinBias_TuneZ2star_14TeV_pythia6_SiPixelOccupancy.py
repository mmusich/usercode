# Based on Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: FourMuPt_1_200_cfi --conditions auto:upgradePLS3 -n 10 --eventcontent FEVTDEBUG,DQM --relval 10000,100 -s GEN,SIM,DIGI:pdigi_valid,L1,DIGI2RAW,RAW2DIGI,L1Reco,RECO,VALIDATION,DQM --datatier GEN-SIM-DIGI-RAW-RECO,DQM --beamspot Gauss --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023Muon --geometry Extended2023Muon,Extended2023MuonReco --magField 38T_PostLS1 --fileout file:step_all.root --no_exec
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing()

options.register('maxEvents',
                 -1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to process (-1 for all)")

options.register('BPixThr',
                 2000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "BPix Clusterizer Threshold (2000 e- is default)")

options.register('PixElePerADC',
                 135,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Pixel Digitizer e per ADC (135 e- is default)")

options.register('PixMaxADC',
                 255,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Pixel Digitizer ADC max (255 / 8 bit is default)")

options.register('ChannelThreshold',
                 1000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Pixel Clusterizer channel threshold (1000 e- is default)")

options.register('SeedThreshold',
                 1000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Pixel Clusterizer seed threshold (1000 e- is default)")

options.register('ClusterThreshold',
                 4000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Pixel Clusterizer cluster threshold (4000 e- is default)")


options.register('AgeingScenario',
                 "NoAgeing", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string,         # string, int, or float
                 "Ageing scenario (NoAgeing is default)")

options.register('MySeed',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "my seed for the job (1 is default)")

options.parseArguments()

print "my seed is: ", options.MySeed

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
# for gaussian smeared vertex
#process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('Configuration.StandardSequences.VtxSmearedNoSmear_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 50
#process.MessageLogger.detailedInfo = cms.untracked.PSet(threshold  = cms.untracked.string('DEBUG'))
#process.MessageLogger.debugModules = cms.untracked.vstring('make_ntuple')

# configure seed
process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(options.MySeed)
process.RandomNumberGeneratorService.generator.engineName = cms.untracked.string('TRandom3')

process.RandomNumberGeneratorService.VtxSmeared.initialSeed = cms.untracked.uint32(options.MySeed)
process.RandomNumberGeneratorService.VtxSmeared.engineName = cms.untracked.string('TRandom3') 

# set maximum events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

myFirstEvent = (options.maxEvents * options.MySeed)+1
print "first Event:",myFirstEvent

# Input source
process.source = cms.Source("EmptySource",
                            firstEvent = cms.untracked.uint32(myFirstEvent),
                            firstLuminosityBlock = cms.untracked.uint32(options.MySeed+1)
                            )

process.options = cms.untracked.PSet(
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('FourMuPt_1_200_cfi nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition
outrootfile='file:MinBias_cff_py_GEN_SIM_RECO_age_'+str(options.AgeingScenario)+'_BPixThr_'+str(options.BPixThr)+'_'+str(options.maxEvents)+'_evts_seed_'+str(options.MySeed)+'.root'
outntuplefile='OccupancyPlotsTest_MinBias_ntuple_age_'+str(options.AgeingScenario)+'_BPixThr_'+str(options.BPixThr)+'_'+str(options.maxEvents)+'_evts_seed_'+str(options.MySeed)+'.root'
print 'output file name:', outrootfile, outntuplefile

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string(outrootfile),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-RECO')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

###### begin of the ntuplizer fragment 
# based on http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Brownson/SLHCUpgradeSimulations/test/resolutionPlotter/
process.ReadLocalMeasurement = cms.EDAnalyzer("StdPixelHitNtuplizer",
   src = cms.InputTag("siPixelRecHits"),
   ### if using simple (non-iterative) or old (as in 1_8_4) tracking
   trackProducer = cms.InputTag("generalTracks"),
   OutputFile = cms.string(outntuplefile),
   ### for using track hit association
   associatePixel = cms.bool(True),
   associateStrip = cms.bool(False),
   associateRecoTracks = cms.bool(False),
   ROUList = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof',
                         'g4SimHitsTrackerHitsPixelBarrelHighTof',
                         'g4SimHitsTrackerHitsPixelEndcapLowTof',
                         'g4SimHitsTrackerHitsPixelEndcapHighTof')
)
# since cmsDriver generates a "Schedule" the ntuplizer needs to be defined as "Path"
process.make_ntuple = cms.Path(process.ReadLocalMeasurement)
###### end of ntuplizer fragment


###### begin occupancy plots fragment
# based on http://cmslxr.fnal.gov/lxr/source/DPGAnalysis/SiStripTools/test/OccupancyPlotsTest_phase2_cfg.py?v=CMSSW_6_2_0_SLHC17
from DPGAnalysis.SiStripTools.occupancyplotsselections_phase2_cff import *

process.spclusmultprod = cms.EDProducer("SiPixelClusterMultiplicityProducer",
                                        clusterdigiCollection = cms.InputTag("siPixelClusters"),
                                        wantedSubDets = cms.VPSet()
                                        )
process.spclusmultprod.wantedSubDets.extend(OccupancyPlotsPixelWantedSubDets)

process.spclusoccuprod = cms.EDProducer("SiPixelClusterMultiplicityProducer",
                                        clusterdigiCollection = cms.InputTag("siPixelClusters"),
                                        withClusterSize = cms.untracked.bool(True),
                                        wantedSubDets = cms.VPSet()
                                        )
process.spclusoccuprod.wantedSubDets.extend(OccupancyPlotsPixelWantedSubDets)

process.seqMultProd = cms.Sequence( process.spclusmultprod + process.spclusoccuprod )

process.load("DPGAnalysis.SiStripTools.occupancyplots_cfi")
process.occupancyplots.file = cms.untracked.FileInPath('SLHCUpgradeSimulations/Geometry/data/PhaseII/Pixel10D/PixelSkimmedGeometry.txt')

process.pixeloccupancyplots = process.occupancyplots.clone()
process.pixeloccupancyplots.wantedSubDets = cms.VPSet()
process.pixeloccupancyplots.wantedSubDets.extend(OccupancyPlotsPixelWantedSubDets)
process.pixeloccupancyplots.multiplicityMaps = cms.VInputTag(cms.InputTag("spclusmultprod"))
process.pixeloccupancyplots.occupancyMaps = cms.VInputTag(cms.InputTag("spclusoccuprod"))

# process.siStripQualityESProducer.ListOfRecordToMerge=cms.VPSet(
# 	#    cms.PSet( record = cms.string("SiStripDetVOffRcd"),    tag    = cms.string("") ),
# 	cms.PSet( record = cms.string("SiStripDetCablingRcd"), tag    = cms.string("") ),
# 	cms.PSet( record = cms.string("RunInfoRcd"),           tag    = cms.string("") ),
# 	cms.PSet( record = cms.string("SiStripBadChannelRcd"), tag    = cms.string("") ),
# 	cms.PSet( record = cms.string("SiStripBadFiberRcd"),   tag    = cms.string("") ),
# 	cms.PSet( record = cms.string("SiStripBadModuleRcd"),  tag    = cms.string("") )
# 	)
 
process.SiStripDetInfoFileReader = cms.Service("SiStripDetInfoFileReader")
process.seqAnalyzers = cms.Path(process.pixeloccupancyplots) 
process.seqProducers = cms.Path(process.seqMultProd)
###### end of occupancy plots

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string(outntuplefile)
                                   )

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(14000.0),
    crossSection = cms.untracked.double(79150000000),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTU(21)=1     ! Check on possible errors during program execution', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)', 
            'MSTP(52)=2     ! work with LHAPDF', 
            'PARP(82)=1.921 ! pt cutoff for multiparton interactions', 
            'PARP(89)=1800. ! sqrts for which PARP82 is set', 
            'PARP(90)=0.227 ! Multiple interactions: rescaling power', 
            'MSTP(95)=6     ! CR (color reconnection parameters)', 
            'PARP(77)=1.016 ! CR', 
            'PARP(78)=0.538 ! CR', 
            'PARP(80)=0.1   ! Prob. colored parton from BBR', 
            'PARP(83)=0.356 ! Multiple interactions: matter distribution parameter', 
            'PARP(84)=0.651 ! Multiple interactions: matter distribution parameter', 
            'PARP(62)=1.025 ! ISR cutoff', 
            'MSTP(91)=1     ! Gaussian primordial kT', 
            'PARP(93)=10.0  ! primordial kT-max', 
            'MSTP(81)=21    ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=4     ! Defines the multi-parton model'),
        processParameters = cms.vstring('MSEL=0         ! User defined processes', 
            'MSUB(11)=1     ! Min bias process', 
            'MSUB(12)=1     ! Min bias process', 
            'MSUB(13)=1     ! Min bias process', 
            'MSUB(28)=1     ! Min bias process', 
            'MSUB(53)=1     ! Min bias process', 
            'MSUB(68)=1     ! Min bias process', 
            'MSUB(92)=1     ! Min bias process, single diffractive', 
            'MSUB(93)=1     ! Min bias process, single diffractive', 
            'MSUB(94)=1     ! Min bias process, double diffractive', 
            'MSUB(95)=1     ! Min bias process'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi_valid)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

######################################################################################

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.seqProducers,process.seqAnalyzers,process.make_ntuple,process.endjob_step,process.FEVTDEBUGoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023Muon 

#call to customisation function cust_2023Muon imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023Muon(process)

###################################################################################################
# N.B. This lines were add from the step3 (DIGI-RAW)
# produced by the cmsDriver command from runTheMatrix.py --what upgrade -l 10000
# taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Recipes620SLHC#Recipes_without_pileup
# If not used the job will crash! (NO MORE IN USE IN 620_SLHC17 ??)
###################################################################################################

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

# 620_SLHC17_patch1 Extended2023 MuonGEM added to digitizer
process.mix.mixObjects.mixSH.crossingFrames.extend(['MuonGEMHits', 'MuonME0Hits'])

if options.AgeingScenario!="NoAgeing":
    # Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.aging
    from SLHCUpgradeSimulations.Configuration.aging import * 

    #call to customisation function customise_aging_500 imported from SLHCUpgradeSimulations.Configuration.aging

    if options.AgeingScenario=="100":
        process = customise_aging_100(process)
    elif options.AgeingScenario=="200":    
        process = customise_aging_200(process)
    elif options.AgeingScenario=="300":    
        process = customise_aging_300(process)
    elif options.AgeingScenario=="400":    
        process = customise_aging_400(process)
    elif options.AgeingScenario=="500":    
        process = customise_aging_500(process)
    elif options.AgeingScenario=="600":    
        process = customise_aging_600(process)
    elif options.AgeingScenario=="700":    
        process = customise_aging_700(process)
    elif options.AgeingScenario=="1000":    
        process = customise_aging_1000(process)
    elif options.AgeingScenario=="3000":    
        process = customise_aging_3000(process)
    else:
        print "Unrecognized Ageing scenario, using default (=NoAgeing)"

# customize the CPE errors
from AuxCode.SLHCSimPhase2.PixelCPE_tables_cff import *
process.PixelCPEGenericESProducer.PixelCPEList = PixelCPE_dict['pixelCPE_100x150_upgrade']

## Superseeded by process = setCrossingFrameOn(process)
# customize to make crossingFrames available (needed for turning on g4SimHits in the ROUList)
#from AuxCode.SLHCSimPhase2.crossingFrameCustoms import *
#customiseCrossingFrame(process)

# Uncomment next two lines to change pixel DIGI threshold
process.mix.digitizers.pixel.ThresholdInElectrons_BPix = cms.double(options.BPixThr)
process.mix.digitizers.pixel.ThresholdInElectrons_BPix_L1 = cms.double(options.BPixThr)

# Added line from A. Tricomi to switch off the pixel inefficiencies from python
process.mix.digitizers.pixel.AddPixelInefficiencyFromPython = cms.bool(False)

# Customise noise
process.mix.digitizers.pixel.AddNoise = cms.bool(False)	
process.mix.digitizers.pixel.AddThresholdSmearing = cms.bool(False)
process.mix.digitizers.pixel.AddNoisyPixels = cms.bool(False)


#process.mix.digitizers.pixel.FluctuateCharge = cms.untracked.bool(False)


# from   SimTracker/SiPixelDigitizer/plugins/SiPixelDigitizerAlgorithm.cc
# // ADC calibration 1adc count(135e.
# // Corresponds to 2adc/kev, 270[e/kev]/135[e/adc](2[adc/kev]
# // Be carefull, this parameter is also used in SiPixelDet.cc to
# // calculate the noise in adc counts from noise in electrons.
# // Both defaults should be the same.
# theElectronPerADC(conf.getParameter<double>("ElectronPerAdc")),
# // ADC saturation value, 255(8bit adc.
# //theAdcFullScale(conf.getUntrackedParameter<int>("AdcFullScale",255)),
# theAdcFullScale(conf.getParameter<int>("AdcFullScale")),
# theAdcFullScaleStack(conf.exists("AdcFullScaleStack")?conf.getParameter<int>("AdcFullScaleStack"):255), <- this is for phase2 PS modules 

# CAREFUL: gain is encoded in
#http://cmslxr.fnal.gov/lxr/source/RecoLocalTracker/SiPixelClusterizer/src/PixelThresholdClusterizer.cc?v=%EF%BB%BFCMSSW_6_2_0_SLHC15#260
process.mix.digitizers.pixel.AdcFullScale = cms.int32(options.PixMaxADC)        # phase1 default 255 (8 bit)
process.mix.digitizers.pixel.ElectronPerAdc = cms.double(options.PixElePerADC)  # phase1 default 135 (255*135=34425)

# Customise clusterizer
process.siPixelClusters.ChannelThreshold = cms.int32(options.ChannelThreshold)  # phase1 default 1000e
process.siPixelClusters.SeedThreshold = cms.int32(options.SeedThreshold)        # phase1 default 1000e
process.siPixelClusters.ClusterThreshold = cms.double(options.ClusterThreshold) # phase1 default 4000e
process.siPixelClusters.ElectronPerAdc = cms.double(options.PixElePerADC)
# End of customisation functions

#print process.dumpPython()
