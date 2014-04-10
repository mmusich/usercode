# Auto generated configuration file
# using: 
# Revision: 1.14 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step_digitodqm --no_exec -s DIGI,L1,DIGI2RAW,RAW2DIGI,L1Reco,RECO,VALIDATION,DQM --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,SLHCUpgradeSimulations/Configuration/phase1TkCustoms.customise,SLHCUpgradeSimulations/Configuration/customise_mixing.customise_NoCrossing --conditions auto:upgrade2017 --datatier GEN-SIM-RECO,DQM -n 100 --geometry Extended2017 --eventcontent FEVTDEBUGHLT,DQM --filein file:/lustre/cms/store/user/traverso/UpgradeSamples/step1_TTtoAnything_1k_evts.root --pileup AVE_140_BX_25ns --pileup_input file:MinBias_TuneZ2star_14TeV_pythia6_cff_py_GEN_SIM.root --fileout file:step_digitodqm.root
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing()

options.register('InputFileName',
                 "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/612_slhc8/Extended2017/TTbar/step1_TTtoAnything_14TeV_pythia6_15k_evts.root", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string,         # string, int, or float
                 "name of the input file ")

options.register('firstEvent',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "First event to process")

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

options.register('PixelCPE',
                 "pixel_CPE_100x150_upgrade", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string,         # string, int, or float
                 "Pixel CPE (pixel_CPE_100x150_upgrade is default)")

options.register('OutFileName',
                 "step_digitodqm.root", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string,         # string, int, or float
                 "name of the output file (step_digitodqm.root is default)")

options.register('PUScenario',
                 "NoPU", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string,         # string, int, or float
                 "PU scenario (NoPileUp is default)")

options.register('AgeingScenario',
                 "NoAgeing", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string,         # string, int, or float
                 "Ageing scenario (NoAgeing is default)")

options.parseArguments()

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')

# import of configurations for PU
if options.PUScenario!="NoPU":

    process.load('SimGeneral.MixingModule.mix_E8TeV_AVE_16_BX_25ns_cfi')
    process.mix.input.fileNames = cms.untracked.vstring(['root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/612_slhc8/Extended2017/MinBias/step1_MinBias_TuneZ2star_14TeV_pythia6_15k_evts.root'])
    process.mix.bunchspace = cms.int32(25)
    process.mix.minBunch = cms.int32(-12)
    process.mix.maxBunch = cms.int32(3)

    if "140" in options.PUScenario:        
        process.mix.input.nbPileupEvents.averageNumber = cms.double(140.000000)
        print "PU = 140"
    elif "10" in options.PUScenario:
        process.mix.input.nbPileupEvents.averageNumber = cms.double(10.000000)
    elif "25" in options.PUScenario:
        process.mix.input.nbPileupEvents.averageNumber = cms.double(25.000000)
    elif "35" in options.PUScenario:
        process.mix.input.nbPileupEvents.averageNumber = cms.double(35.000000)
    elif "50" in options.PUScenario:
        process.mix.input.nbPileupEvents.averageNumber = cms.double(50.000000)
    elif "70" in options.PUScenario:
        process.mix.input.nbPileupEvents.averageNumber = cms.double(70.000000)
    elif "75" in options.PUScenario:
        process.mix.input.nbPileupEvents.averageNumber = cms.double(75.000000)
    elif "100" in options.PUScenario:
        process.mix.input.nbPileupEvents.averageNumber = cms.double(100.000000)
    elif "125" in options.PUScenario:
        process.mix.input.nbPileupEvents.averageNumber = cms.double(125.000000)
    elif "150" in options.PUScenario:
        process.mix.input.nbPileupEvents.averageNumber = cms.double(150.000000)
    elif "175" in options.PUScenario:
        process.mix.input.nbPileupEvents.averageNumber = cms.double(175.000000)
    elif "200" in options.PUScenario:
        process.mix.input.nbPileupEvents.averageNumber = cms.double(200.000000)
    else:
        print "Unrecognized PU scenario, using default (=NoPU)"
        process.load('SimGeneral.MixingModule.mixNoPU_cfi')
       
else:
    process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
                            secondaryFileNames = cms.untracked.vstring(),
                            fileNames = cms.untracked.vstring(options.InputFileName),
                            firstEvent = cms.untracked.uint32(options.firstEvent)
                            )

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.14 $'),
    annotation = cms.untracked.string('step_digitodqm nevts:100'),
    name = cms.untracked.string('Applications')
)

# Output definition

outputRootFileName=cms.untracked.string(options.OutFileName)

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = outputRootFileName,
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)

DqmFileName=str(options.OutFileName)
NewDqmFileName=DqmFileName.replace("digitodqm","digitodqm_inDQM")
dqmRootFileName=cms.untracked.string(NewDqmFileName)

process.DQMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.DQMEventContent.outputCommands,
    fileName = dqmRootFileName,
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('DQM')
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2017', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.prevalidation_step = cms.Path(process.prevalidation)
process.dqmoffline_step = cms.Path(process.DQMOffline)
process.validation_step = cms.EndPath(process.validation)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,
                                process.L1simulation_step,
                                process.digi2raw_step,
                                process.raw2digi_step,
                                process.L1Reco_step,
                                process.reconstruction_step,
                                process.prevalidation_step,
                                process.validation_step,
                                process.dqmoffline_step,
                                process.endjob_step,
                                process.FEVTDEBUGHLToutput_step,
                                process.DQMoutput_step)


# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.phase1TkCustoms
from SLHCUpgradeSimulations.Configuration.phase1TkCustoms import customise 

#call to customisation function customise imported from SLHCUpgradeSimulations.Configuration.phase1TkCustoms
process = customise(process)

# call to customisation function customise imported from SLHCUpgradesSimulations.Configuration.aging
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

# for light sequence (keeping only tracker validation sequences)
from AuxCode.SLHCSimPhase2.TkOnlyValidationCustoms import customise_tkonly
process = customise_tkonly(process)

# read-in the parametrization of pixel CPE generic from external cff file
# TODO 1) the choice of the array shotld be steerable from somewhere
#      2) the next lines could be fed into the process via the customize option of the cmsDriver.py 

print options.PixelCPE
from AuxCode.SLHCSimPhase2.PixelCPE_tables_cff import *
process.PixelCPEGenericESProducer.PixelCPEList = PixelCPE_dict[options.PixelCPE]

# Uncomment next two lines to change pixel DIGI threshold
process.mix.digitizers.pixel.ThresholdInElectrons_BPix = cms.double(options.BPixThr)
process.mix.digitizers.pixel.ThresholdInElectrons_BPix_L1 = cms.double(options.BPixThr)

# End of customisation functions
