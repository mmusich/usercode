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

# fileNames_4mu = cms.untracked.vstring("root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_10_evts_seed_0.root",
#                                       "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_10_evts_seed_1.root")

fileNames_4mu = cms.untracked.vstring("root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_0.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_1.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_10.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_11.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_12.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_13.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_14.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_15.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_16.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_17.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_18.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_19.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_2.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_20.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_21.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_22.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_23.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_24.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_25.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_26.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_27.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_28.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_29.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_3.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_30.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_31.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_32.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_33.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_34.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_35.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_36.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_37.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_38.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_39.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_4.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_5.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_6.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_7.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_8.root",
                                      "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_FourMuPartGun_2500_evts_seed_9.root"
                                      )

# fileNames_minbias = cms.untracked.vstring("root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_10_evts_seed_0.root",
#                                           "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_10_evts_seed_1.root"
#                                           )

fileNames_minbias = cms.untracked.vstring("root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_0.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_1.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_10.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_11.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_12.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_13.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_14.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_15.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_16.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_17.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_18.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_19.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_2.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_20.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_21.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_22.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_23.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_24.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_25.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_26.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_27.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_28.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_29.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_3.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_30.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_31.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_32.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_33.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_34.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_35.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_36.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_37.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_38.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_39.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_4.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_5.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_6.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_7.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_8.root",
                                          "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_MinBias_TuneZ2star_14TeV_pythia6_300_evts_seed_9.root"
                                          )

# fileNames_ttbar = cms.untracked.vstring("root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_10_evts_seed_0.root",
#                                         "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_10_evts_seed_1.root"
#                                         )

fileNames_ttbar = cms.untracked.vstring("root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_0.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_1.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_10.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_11.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_12.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_13.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_14.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_15.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_16.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_17.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_18.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_19.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_2.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_20.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_21.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_22.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_23.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_24.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_25.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_26.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_27.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_28.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_29.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_3.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_30.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_31.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_32.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_33.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_34.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_35.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_36.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_37.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_38.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_39.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_4.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_5.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_6.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_7.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_8.root",
                                        "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/620_slhc11/ExtendedPhaseIPixel/step1_TTtoAnything_300_evts_seed_9.root"
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
