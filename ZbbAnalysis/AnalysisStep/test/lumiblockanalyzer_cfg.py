import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

###################################################################
# Setup 'standard' options
###################################################################
options = VarParsing.VarParsing()
options.register('isMuon',True,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool,"True: run on muons; False: run on electrons")
options.register('maxEvents',-1,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int,"Number of events to process (-1 for all)")
options.parseArguments()

###################################################################
# Messages
###################################################################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

#process.load("ZbbAnalysis.AnalysisStep.muon_ZbbSkimDec22ReReco_PAT397_08Apr11_cff")

###################################################################
# Source files
###################################################################

electronSource=[
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_1_1_4IW.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_2_1_Uce.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_3_1_Qhk.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_4_1_Nuq.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_5_1_Khn.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_6_1_dg5.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_7_1_bJK.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_8_1_cbt.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_9_2_LrS.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_10_3_pL3.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_11_1_zYu.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_12_2_efg.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_13_1_jwC.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_14_2_UEa.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_15_2_kEJ.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_16_2_vl4.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_17_1_bNs.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_18_2_YvL.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_19_1_Vmj.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_20_1_MkM.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_21_1_pZf.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_22_1_iQD.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_23_1_K2F.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_24_1_w6G.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_25_1_X2o.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_26_1_KuX.root'
    ] 

muonSource=[
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_1_1_4IW.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_2_1_Uce.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_3_1_Qhk.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_4_1_Nuq.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_5_1_Khn.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_6_1_dg5.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_7_1_bJK.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_8_1_cbt.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_9_2_LrS.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_10_3_pL3.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_11_1_zYu.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_12_2_efg.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_13_1_jwC.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_14_2_UEa.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_15_2_kEJ.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_16_2_vl4.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_17_1_bNs.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_18_2_YvL.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_19_1_Vmj.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_20_1_MkM.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_21_1_pZf.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_22_1_iQD.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_23_1_K2F.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_24_1_w6G.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_25_1_X2o.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_26_1_KuX.root',
    #muons start here:
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_1_1_Uui.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_2_1_lUC.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_3_1_iNw.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_4_1_Z8m.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_5_1_44d.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_6_1_vyd.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_7_1_717.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_8_1_6bu.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_9_1_2K8.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_10_1_lBg.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_11_1_s3U.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_12_1_SfB.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_13_2_hi7.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_14_2_7Yq.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_15_2_ffA.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_16_3_JAL.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_17_1_Osd.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_18_2_T0M.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_19_2_jlN.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_20_2_zac.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_21_1_WYP.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_22_3_3PF.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_23_1_1xg.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_24_1_1GM.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_25_2_s8x.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_26_2_ZdI.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_27_1_WbB.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_28_2_x4D.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_29_3_ld4.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_30_2_XwM.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_31_2_LRz.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_32_4_je7.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_33_3_LKF.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_34_1_KNH.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_35_2_QAN.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_36_2_lzV.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_37_2_z2G.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_38_1_neq.root'
    ]

readFiles = cms.untracked.vstring()

if options.isMuon:
    readFiles.extend(muonSource)

else:
    readFiles.extend(electronSource)    

process.source = cms.Source("PoolSource",
                            fileNames = readFiles ,
                            duplicateCheckMode = cms.untracked.string('checkAllFilesOpened')
                            )

process.load("ZbbAnalysis.AnalysisStep.lumiblockanalyzer_cfi")

process.TFileService=cms.Service('TFileService',
                                 fileName=cms.string('MyPlots.root')
                                 )

process.p = cms.Path(process.lumiblockanalysis)
