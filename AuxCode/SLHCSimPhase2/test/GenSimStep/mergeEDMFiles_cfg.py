import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# Give the process a name
process = cms.Process("MergeEDMFiles")

###################################################################
# Setup 'standard' options
###################################################################
options = VarParsing.VarParsing()

options.register('input',
                 "TTbar", 
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string, # string, int, or float
                 "name of the sample to be merged")

options.parseArguments()

fnames = None


fileNames_minbias = cms.untracked.vstring("root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_0.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_1.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_2.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_3.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_4.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_5.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_6.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_7.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_8.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_9.root"
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_10.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_11.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_12.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_13.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_14.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_15.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_16.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_17.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_18.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_19.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_20.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_21.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_22.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_23.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_24.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_25.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_26.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_27.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_28.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_29.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_30.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_31.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_32.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_33.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_34.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_35.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_36.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_37.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_38.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_39.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_40.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_41.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_42.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_43.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_44.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_45.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_46.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_47.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_48.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_49.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_50.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_51.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_52.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_53.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_54.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_55.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_56.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_57.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_58.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_59.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_60.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_61.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_62.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_63.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_64.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_65.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_66.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_67.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_68.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_69.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_70.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_71.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_72.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_73.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_74.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_75.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_76.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_77.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_78.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_79.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_80.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_81.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_82.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_83.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_84.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_85.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_86.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_87.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_88.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_89.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_90.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_91.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_92.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_93.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_94.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_95.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_96.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_97.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_98.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_99.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_100.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_101.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_102.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_103.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_104.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_105.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_106.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_107.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_108.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_109.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_110.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_111.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_112.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_113.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_114.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_115.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_116.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_117.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_118.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_119.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_120.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_121.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_122.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_123.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_124.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_125.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_126.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_127.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_128.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_129.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_130.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_131.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_132.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_133.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_134.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_135.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_136.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_137.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_138.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_139.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_140.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_141.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_142.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_143.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_144.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_145.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_146.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_147.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_148.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_149.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_150.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_151.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_152.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_153.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_154.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_155.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_156.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_157.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_158.root",
                                          # "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc17_patch1/Extended2023Muon/Thick_L0L1_0.285/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_159.root"
                                          )

outrootfile='file:/tmp/emiglior/step1_'+str(options.input)+'_merged.root'
print 'output file name:', outrootfile

if options.input == 'TTbar' :
    fnames = fileNames_ttbar
    print "selecting ttbar"
elif options.input == 'muons' :
    fnames = fileNames_4mu
    print "selecting muons"
elif options.input == 'MinBias' :
    fnames = fileNames_minbias
    print "selecting minbias"
else :
    print "unrecognized input sample"
    
# Tell the process which files to use as the source
process.source = cms.Source ("PoolSource",
                             fileNames = fnames 
                             )

# tell the process to only run over 100 events (-1 would mean run over
#  everything
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32 (-1)
)

# Tell the process what filename to use to save the output
process.Out = cms.OutputModule("PoolOutputModule",
         fileName = cms.untracked.string (outrootfile)
)

# make sure everything is hooked up
process.end = cms.EndPath(process.Out)
