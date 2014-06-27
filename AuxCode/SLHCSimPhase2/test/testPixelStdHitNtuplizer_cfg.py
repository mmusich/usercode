# Based on Auto generated configuration file
# using: 
# Revision: 1.14 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/FourteenTeV/TenMuE_0_200_cfi.py --no_exec -s GEN,SIM,DIGI:pdigi_valid,L1,DIGI2RAW,RAW2DIGI,L1Reco,RECO --conditions auto:upgrade2017 --magField 38T_PostLS1 --eventcontent FEVTDEBUG --beamspot NoSmear --geometry Extended2017 --relval 10000,100 --datatier GEN-SIM-RECO -n 500  --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2017 --fileout file:TenMuE_0_200_cff_py_GEN_SIM_RECO.root

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
process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2017_cff')
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
print "firs Event:",myFirstEvent

# Input source
process.source = cms.Source("EmptySource",
                            firstEvent = cms.untracked.uint32(myFirstEvent),
                            firstLuminosityBlock = cms.untracked.uint32(options.MySeed+1)
                            )

process.options = cms.untracked.PSet(
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.14 $'),
    annotation = cms.untracked.string('Configuration/GenProduction/python/FourteenTeV/TenMuE_0_200_cfi.py nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition

outrootfile='file:TenMuE_0_200_cff_py_GEN_SIM_RECO_age_'+str(options.AgeingScenario)+'_BPixThr_'+str(options.BPixThr)+'_'+str(options.maxEvents)+'_evts_seed_'+str(options.MySeed)+'.root'
outntuplefile='stdgrechitfullph1g_ntuple_age_'+str(options.AgeingScenario)+'_BPixThr_'+str(options.BPixThr)+'_'+str(options.maxEvents)+'_evts_seed_'+str(options.MySeed)+'.root'
print 'output file name:', outrootfile, outntuplefile

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string(outrootfile),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statement
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2017', '')

process.generator = cms.EDProducer("FlatRandomEGunProducer",
    PGunParameters = cms.PSet(
        PartID = cms.vint32(-13, -13, -13, -13, -13),
	#PartID = cms.vint32(-13),
        MaxEta = cms.double(2.7),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-2.7),
        MinE = cms.double(0.0),
        MinPhi = cms.double(-3.14159265359),
        MaxE = cms.double(200.0),
	MinPt = cms.double(0.9)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('Ten mu e 0 to 200'),
    AddAntiParticle = cms.bool(True),				   
    firstRun = cms.untracked.uint32(1)
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
### This fragment is meant to produce ntuples for the calibration of the pixel CPE ###
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Brownson/SLHCUpgradeSimulations/test/resolutionPlotter/
process.ReadLocalMeasurement = cms.EDAnalyzer("PixelStdHitNtuplizer",
   src = cms.InputTag("siPixelRecHits"),
   stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
   rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
   matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
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
######################################################################################

# Schedule definition
# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.make_ntuple,process.endjob_step,process.FEVTDEBUGoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2017 

#call to customisation function cust_2017 imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2017(process)

###################################################################################################
# N.B. This lines were add from the step3 (DIGI-RAW)
# produced by the cmsDriver command from runTheMatrix.py --what upgrade -l 10000
# taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Recipes620SLHC#Recipes_without_pileup
# If not used the job will crash!
###################################################################################################

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

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


#print process.dumpPython()
# End of customisation functions

