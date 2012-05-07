import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("TEST")

options = VarParsing.VarParsing()
options.register('Sample',
                 "MadGraph", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string,         # string, int, or float
                 "Sample to be processed (ZIncl,MadGraph,Sherpa,aMCatNLO)")
options.register('maxEvents',
                 -1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to process (-1 for all)")
options.register('channel',
                 "Z_m",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Z decay channel? (Z_e,Z_m)")
options.parseArguments()

# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))
 
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'START42_V13::All'

ZjetsInclSrc = [
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_100_1_mjJ.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_104_1_CSZ.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_109_1_lcq.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_111_1_4rL.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_113_5_e0d.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_117_5_syP.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_118_5_Ifc.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_11_1_P45.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_121_1_wWm.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_122_1_0Ln.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_128_1_ODY.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_130_4_Mib.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_146_4_2pN.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_151_2_e6d.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_158_4_Ewk.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_161_6_ira.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_162_4_LzI.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_1_1_ug7.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_28_4_5aJ.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_2_1_Bw5.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_39_4_hQG.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_48_5_VjS.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_49_5_5ii.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_50_5_ysU.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_51_5_ycw.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_53_5_W9Z.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_60_2_mdT.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_72_5_Fbk.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_73_5_TeW.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_76_5_GjH.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_85_1_yQv.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_8_4_ksk.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_93_1_ZJU.root'
    ]

SherpaSrc = [
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_10_2_179.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_11_3_gDi.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_12_2_rdZ.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_13_3_BmO.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_14_3_OB6.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_15_2_mON.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_16_2_GKJ.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_17_3_nUD.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_18_2_D8S.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_19_3_s4T.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_1_3_KkT.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_20_2_sSs.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_21_2_A4l.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_22_3_pDS.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_23_2_Dbl.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_24_3_H9s.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_25_1_vvi.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_26_1_xGf.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_27_2_PRJ.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_28_2_xSE.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_29_1_P0y.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_2_3_JT1.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_30_2_GAI.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_31_2_F5k.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_32_2_ObC.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_33_1_Esr.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_34_2_q7S.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_35_2_hRQ.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_36_2_rti.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_37_2_iht.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_38_1_Vk4.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_3_3_cnb.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_4_2_lKb.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_5_3_0PV.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_6_1_3qS.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_7_1_Tkk.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_8_3_qeG.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_9_3_rqB.root'
    ]

MadGraphSrc = [
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_10_2_IOo.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_11_2_zfe.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_12_2_0i3.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_13_2_oFF.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_14_2_WxD.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_15_2_V9P.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_16_2_ySu.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_17_1_Et3.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_18_2_jtw.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_19_1_jVD.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_1_2_PX9.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_20_2_IR4.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_21_2_Tu2.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_22_2_dGJ.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_23_1_pmz.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_24_2_HqO.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_25_1_hxw.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_26_1_bgq.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_27_1_GPa.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_28_1_CR9.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_29_1_cCs.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_2_2_1yW.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_30_1_k7E.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_31_1_tqe.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_32_2_WpE.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_33_1_nMC.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_34_1_uBP.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_35_1_UZE.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_36_1_rbC.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_37_1_dqJ.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_38_1_MNR.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_39_1_0o2.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_3_1_w2Z.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_4_1_qIM.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_5_1_oz4.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_6_2_Fv1.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_7_2_oVO.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_8_1_1c1.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_9_2_S0I.root'
    ]

aMCatNLOSrc = [
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_10_1_GM3.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_11_1_luU.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_12_1_4dS.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_13_1_YZB.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_14_1_5S6.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_15_1_Lri.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_16_1_HJ8.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_17_3_4bl.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_18_1_IiM.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_19_1_vkw.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_1_1_rKV.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_20_1_Pa2.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_21_1_3if.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_22_3_VME.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_23_1_uE6.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_24_3_C5Y.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_25_1_dPO.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_26_1_DRt.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_27_1_rYT.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_28_3_moU.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_29_3_oTS.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_2_1_mWo.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_30_3_l6d.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_31_3_Wvo.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_32_1_xeR.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_33_1_bDC.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_34_1_kwV.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_35_3_6a5.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_36_1_091.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_37_1_1g9.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_38_3_Zfb.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_39_3_Kj4.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_3_1_8Wp.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_40_1_R9N.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_41_1_hQS.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_42_1_o4V.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_43_3_i9j.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_44_1_STU.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_45_1_j0U.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_46_1_omL.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_47_1_TO0.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_48_1_TXs.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_49_1_MG6.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_4_1_k6w.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_50_1_OU8.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_51_1_0YC.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_52_1_Vey.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_53_1_s5q.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_54_1_Ctt.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_55_1_kwD.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_56_1_Sjo.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_57_1_arq.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_58_1_2oT.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_59_1_VzW.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_5_1_XfD.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_60_1_0Cy.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_61_1_h4c.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_62_1_L6h.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_6_1_cbz.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_7_1_Kk4.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_8_1_rnF.root',
    'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_9_1_Ld9.root'
    ]

readFiles = cms.untracked.vstring()

process.source = cms.Source("PoolSource",
                            fileNames = readFiles
                            #firstEvent = cms.untracked.uint32(6391783),
                            #firstRun = cms.untracked.uint32(1)
                            )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

###############################################################
# Producer: produce needed collections of genJets
###############################################################
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.genParticlesForPartonJets = process.genParticlesForJets.clone()
process.genParticlesForPartonJets.partonicFinalState = True
process.genParticlesForPartonJets.excludeFromResonancePids = cms.vuint32(11, 12, 13, 14, 15, 16)

process.genParticlesForJets.ignoreParticleIDs = cms.vuint32(
    1000022, 2000012, 2000014,
    2000016, 1000039, 5000039,
    4000012, 9900012, 9900014,
    9900016, 39, 12, 14, 16
    )

process.load("RecoJets.JetProducers.ak5GenJets_cfi")
process.ak5GenJetsNoNuBSM  =  process.ak5GenJets

process.ak7GenJetsNoNuBSM  =  process.ak5GenJets.clone()
process.ak7GenJetsNoNuBSM.rParam = cms.double(0.7)

process.ak5PartonJets  =  process.ak5GenJets.clone()
process.ak5PartonJets.src = cms.InputTag("genParticlesForPartonJets")

process.ak7PartonJets  =  process.ak5GenJets.clone()
process.ak7PartonJets.rParam = cms.double(0.7)
process.ak7PartonJets.src = cms.InputTag("genParticlesForPartonJets")

process.pgenjets = cms.Sequence(
    process.genParticlesForJets
    *process.ak5GenJetsNoNuBSM
    *process.ak7GenJetsNoNuBSM
    *process.genParticlesForPartonJets
    #   *process.ak5PartonJets   #exclude for the moment
    #   *process.ak7PartonJets
    )

###############################################################
# Producer: create new collections of leptons from HF semileptonic decays
###############################################################
# process.LeptFromDecayProducer = cms.EDProducer("LeptFromDecayProducer",
#                                             MuonCollection = cms.InputTag("patMuonsWithTrigger")
#     )

###############################################################
#MC filter: select events with muons from b semileptonic decays 
###############################################################
# process.MCFilter = cms.EDFilter("LeptFromDecayFilter",
#                                 src = cms.InputTag("prunedGen"),
#                                 JetCollection = cms.InputTag("cleanPatJetsPF"),
#                                 ElectronCollection = cms.InputTag("patElectronsWithTrigger"),
#                                 MuonCollection = cms.InputTag("patMuonsWithTrigger"),
#                                 ZmmCollection = cms.InputTag('zMMCand'),
#                                 DecayChainSelection = cms.untracked.string('c>m')
#                                 )

####################################################
#MC Analyzer: spectra of MC Truth Analyzer
####################################################
# process.LeptFromDecayAnalyzer = cms.EDAnalyzer("LeptFromDecayAnalyzer",
#                                                src = cms.InputTag("LeptFromDecayProducer"),
#                                                JetCollection = cms.InputTag("patJets"),
#                                                ElectronCollection = cms.InputTag("patElectronsWithTrigger"),
#                                                MuonCollection = cms.InputTag("patMuonsWithTrigger"),
#                                                genParticleCollection = cms.InputTag("prunedGen"),
#                                                MuonFromBCollection = cms.InputTag("LeptFromDecayProducer:patMuonsFromB:TEST"),
#                                                MuonFromCCollection = cms.InputTag("LeptFromDecayProducer:patMuonsFromC:TEST"),
#                                                MuonFromBCCollection = cms.InputTag("LeptFromDecayProducer:patMuonsFromBC:TEST"),
#                                                ZmmCollection = cms.InputTag('zMMCand'),
#                                                HFSelection = cms.untracked.string("b")
#                                                )

# process.out = cms.OutputModule("PoolOutputModule",
#                                fileName = cms.untracked.string('/tmp/musich/ZbbToLL_PAT_NewCollections.root'),
#                                outputCommands = cms.untracked.vstring('keep *')
#                                )

###################################################################
# Jet Cleaning
###################################################################
process.cleanPatJets = cms.EDProducer("PATJetCleaner",
                                      src = cms.InputTag("patJets"),
                                      #preselection  = cms.string(''),
                                      preselection = cms.string('pt > 20.0 && abs(eta) < 2.4'),
                                      checkOverlaps = cms.PSet(ele = cms.PSet(src       = cms.InputTag("offlinePrimaryVertexFromZ:patElectronsFromZ"), #Z daughters
                                                                              #src       = cms.InputTag("goldenElectrons"),
                                                                              algorithm = cms.string("byDeltaR"),
                                                                              preselection        = cms.string(""),
                                                                              deltaR              = cms.double(0.),
                                                                              checkRecoComponents = cms.bool(False),
                                                                              pairCut             = cms.string(""),
                                                                              requireNoOverlaps = cms.bool(True)
                                                                              ),
                                                               mu = cms.PSet(src       = cms.InputTag("offlinePrimaryVertexFromZ:patMuonsFromZ"), #Z daughters
                                                                             #src       = cms.InputTag("goldenMuons"),
                                                                             algorithm = cms.string("byDeltaR"),
                                                                             preselection        = cms.string(""),
                                                                             deltaR              = cms.double(0.),
                                                                             checkRecoComponents = cms.bool(False),
                                                                             pairCut             = cms.string(""),
                                                                             requireNoOverlaps = cms.bool(True)
                                                                             ),
                                                               ),
                                      finalCut = cms.string(''),
                                      )

process.JetCleaningSeq = cms.Sequence( process.cleanPatJets )

###################################################################
# HF Filter
###################################################################
from ZbbAnalysis.AnalysisStep.HeavyFlavorFilter_cfi import heavyflavorfilter

process.b_heavyflavorfilter= heavyflavorfilter.clone(
    src= cms.InputTag("genParticles"),
    status2 = cms.bool(True),
    status3 = cms.bool(True),
    bOn = cms.bool(True),
    cOn = cms.bool(False),
    hDaughterVeto = cms.bool(False),
    zDaughterVeto = cms.bool(False),
    ptcut=cms.double(0)
    )

###################################################################
# MC Truth Analyzer
###################################################################
# process.MCTruthAnalyzer = cms.EDAnalyzer("MCTruthAnalyzer",
#                                          JetCollection = cms.InputTag("cleanPatJets"),
#                                          genJetSrc = cms.InputTag("patJets:genJets"),
#                                          src = cms.InputTag("genParticles"),
#                                          ElectronCollection = cms.InputTag("patElectrons"),
#                                          MuonCollection = cms.InputTag("patMuons"),
#                                          ZmmCollection = cms.InputTag('zMMCand'),
#                                          isMCatNLO = cms.untracked.bool(False)
#                                          )

# process.ZLONLOHistogrammer = cms.EDAnalyzer("ZLONLOHistogrammer",
#                                             genParticles = cms.InputTag("genParticles"),
#                                             nbinsMass = cms.untracked.uint32(100),
#                                             nbinsPt = cms.untracked.uint32(100),
#                                             nbinsAng = cms.untracked.uint32(100),
#                                             massMax = cms.untracked.double(300.),
#                                             ptMax = cms.untracked.double(200.),
#                                             angMax = cms.untracked.double(3.142),
#                                             accPtMin = cms.untracked.double(20.),
#                                             accMassMin = cms.untracked.double(40.),
#                                             accMassMax = cms.untracked.double(300.),
#                                             accEtaMin = cms.untracked.double(-2.1),
#                                             accEtaMax = cms.untracked.double(2.1),
#                                             isMCatNLO = cms.untracked.bool(False)
#                                             )

#################################################################
# SFb from DB                      
##################################################################
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")
process.load("RecoBTag.PerformanceDB.BTagPerformanceDB1107")

###################################################################
# Import b-efficiency scaling factors
###################################################################
from ZbbAnalysis.AnalysisStep.effbctables_cff import mc_effbc_ssvhem_2011
from ZbbAnalysis.AnalysisStep.effbctables_cff import mc_effbc_ssvhpt_2011

###################################################################
# Acceptance Calculator
###################################################################
process.calculateAcceptance = cms.EDAnalyzer("ZbbAcceptanceCalculator",
                                             GenSrc = cms.InputTag("genParticles"),
                                             ElectronCollection = cms.InputTag("patElectrons"),
                                             MuonCollection = cms.InputTag("patMuons"),
                                             JetCollection = cms.InputTag("cleanPatJets"),
                                             # genJetSrc = cms.InputTag("patJets:genJets"),
                                             genJetSrc  = cms.InputTag("ak5GenJetsNoNuBSM"),
                                             ZmmCollection = cms.InputTag('zMMCand'),
                                             ZeeCollection = cms.InputTag('zEECand'),
                                             unLockDefaultCuts = cms.bool(False),  # bool to unlock default acceptance cuts
                                             useStatus3forMuons = cms.bool(True), 	    
                                             useStatus3forElectrons = cms.bool(True), 	    
                                             doFSRCorrectionForMuons = cms.bool(False), 
                                             doFSRCorrectionForElectrons = cms.bool(False),
                                             jetEtaCut   = cms.double(2.1),
                                             muonEtaCut  = cms.double(2.1),
                                             eleEtaCut   = cms.double(2.5),
                                             jetPtCut    = cms.double(25.),
                                             genJetPtCut = cms.double(0.),
                                             muonPtCut   = cms.double(20.),
                                             elePtCut    = cms.double(25.),
                                             minMassCut  = cms.double(60.),
                                             maxMassCut  = cms.double(120.),
                                             dRLeptonMatch = cms.double(0.3),
                                             dRJetMatch    = cms.double(0.5),
                                             BparticlePtCut= cms.double(0.),
                                             isMCatNLO     = cms.untracked.bool(False),
                                             applyPUCorr   = cms.untracked.bool(True),
                                             verbose       = cms.untracked.bool(False),
                                             useClopperPearsonErrors = cms.untracked.bool(True),
                                             saveNTuple    = cms.untracked.bool(True),
                                             PartonLevel   = cms.untracked.bool(False),
                                             DecayChainSelection = cms.untracked.string(options.channel),
                                             OutfileName = cms.string('CorrectionFactors_'+options.channel+'.txt'),
                                             EffbcMCSSVHE = mc_effbc_ssvhem_2011,
                                             EffbcMCSSVHP = mc_effbc_ssvhpt_2011
                                             )

##################################################################
# hack to force NLO weights rule
##################################################################
if options.Sample=="aMCatNLO":
    # process.ZLONLOHistogrammer = cms.untracked.bool(True)     
    # process.MCTruthAnalyzer.isMCatNLO = cms.untracked.bool(True)
    process.calculateAcceptance.isMCatNLO = cms.untracked.bool(True)   

###################################################################
# Definition of analysis sequence
###################################################################
process.AnalysisSequence   = cms.Sequence()

process.calcAcc= cms.Sequence(process.pgenjets*
                              #process.b_heavyflavorfilter*
                              process.JetCleaningSeq*
                              #process.ZLONLOHistogrammer* 
                              #process.MCTruthAnalyzer*
                              process.calculateAcceptance
                              )

##############################################################################
# Patch for sum(E) violation in PYTHIA/MADGRAPH DY samples
# https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/1489.html
##############################################################################
process.load("ZbbAnalysis.Tools.TotalKinematicsFilter_cfi")

if options.Sample=="MadGraph":
    process.AnalysisSequence+=process.totalKinematicsFilter

###################################################################
# switches
###################################################################
if options.Sample=="ZIncl":    
    readFiles.extend(ZjetsInclSrc)
    process.AnalysisSequence+=process.calcAcc
    if options.channel=="Z_m":
        outputRootFileName=cms.string('Acceptance_ZInclToMM_MADGRAPH.root')
    elif options.channel=="Z_e":
        outputRootFileName=cms.string('Acceptance_ZInclToEE_MADGRAPH.root')
elif options.Sample=="MadGraph":    
    readFiles.extend(MadGraphSrc)
    process.AnalysisSequence+=process.calcAcc
    if options.channel=="Z_m":
        outputRootFileName=cms.string('Acceptance_ZbbToMM_MADGRAPH.root')
    elif options.channel=="Z_e":
        outputRootFileName=cms.string('Acceptance_ZbbToEE_MADGRAPH.root')
elif options.Sample=="Sherpa":
    process.AnalysisSequence+=process.calcAcc
    readFiles.extend(SherpaSrc)
    if options.channel=="Z_m":
        outputRootFileName=cms.string('Acceptance_ZbbToMM_SHERPA.root')
    elif options.channel=="Z_e":
        outputRootFileName=cms.string('Acceptance_ZbbToEE_SHERPA.root')
elif options.Sample=="aMCatNLO":
    process.AnalysisSequence+=process.calcAcc
    readFiles.extend(aMCatNLOSrc)
    if options.channel=="Z_m":
        outputRootFileName=cms.string('Acceptance_ZbbToMM_aMCatNLO.root')
    elif options.channel=="Z_e":
        outputRootFileName=cms.string('Acceptance_ZbbToEE_aMCatNLO.root')

###################################################################
# output file
###################################################################
process.TFileService = cms.Service("TFileService", fileName = outputRootFileName)

###################################################################
# partcile list drawer
###################################################################
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                   maxEventsToPrint = cms.untracked.int32(-1),
                                   printVertex = cms.untracked.bool(False),
                                   src = cms.InputTag("genParticles")
                                   )

###################################################################
# check if events in the list
###################################################################
process.checkList = cms.EDFilter("EventListFilter",
                                 Inputfile=cms.untracked.string("./badEvents/listOfBadEvents.txt")
                                 )

process.p = cms.Path(process.AnalysisSequence)
#process.p = cms.Path(process.checkList*process.AnalysisSequence)
#process.p = cms.Path(process.checkList*process.printTree)
#process.p = cms.Path(process.b_heavyflavorfilter)
#process.e = cms.EndPath(process.out)




