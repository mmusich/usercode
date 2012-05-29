import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("TEST")

###################################################################
# Setup 'standard' options
###################################################################
options = VarParsing.VarParsing()
options.register('Sample',
                 "Zl", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string,         # string, int, or float
                 "Sample to be processed (Zl,Zb5f,Zb4f,Zb5fSherpa,Zb5faMCatNLO,Zc,tt,Ztautau,zz,wz)")
options.register('maxEvents',
                 -1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to process (-1 for all)")
options.parseArguments()

###################################################################
# Messages
###################################################################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

###################################################################
# Standard loads
###################################################################
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

###################################################################
# Global tag
###################################################################
#process.GlobalTag.globaltag = 'GR_R_39X_V5::All' # for the DATA
process.GlobalTag.globaltag = 'START42_V13::All'  # for MC

###################################################################
# Source files
###################################################################

# 2010 samples
#DYSource= ['rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/DYJetsToLL_TuneZ2/MergedOutputFile_1_1_hFC.root']
#ZbbSource=['rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/ZbbToLL/MergedOutputFile_1_1_Q6N.root']
#ZccSource=['rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/ZccToLL/MergedOutputFile_1_1_dh4.root']
#TTSource=['rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimDec22ReReco_PAT397_08Apr11_V2/TTJets_TuneZ2/MergedOutputFile_1_1_85V.root']

# new 2010 samples
#DYSource=[
#'rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimMCWinter10_PAT397_22Jun11_V2/DYJetsToLL_TuneZ2_split/MergedOutputFile_1_1_9np.root',
#'rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimMCWinter10_PAT397_22Jun11_V2/DYJetsToLL_TuneZ2_split/MergedOutputFile_2_1_xzP.root',
#'rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimMCWinter10_PAT397_22Jun11_V2/DYJetsToLL_TuneZ2_split/MergedOutputFile_3_1_GNk.root'
#]
#ZbbSource=['rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimMCWinter10_PAT397_22Jun11_V2/ZbbToLL/MergedOutputFile_1_1_QYJ.root']
#ZccSource=['rfio:/castor/cern.ch/cms/store/user/emiglior/zbb/ZbbSkimMCWinter10_PAT397_22Jun11_V2/ZccToLL/MergedOutputFile_1_1_13S.root']

# 2011 samples

DYSource_MADGRAPH_5FS=['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_1_1_lot.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_2_1_epd.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_3_1_2XE.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_4_1_71x.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_5_1_x1F.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_6_1_VYR.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_7_1_spv.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_8_1_uuy.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_9_1_Z8W.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_10_1_2R9.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_11_2_UeP.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_12_1_Bz5.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_13_1_bmU.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_14_1_tYE.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_15_1_dWa.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_16_1_1Yt.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_17_1_i5L.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_18_2_Pk3.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_19_1_Yl8.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_20_2_QIG.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_21_1_RVd.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_22_1_Ews.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_23_1_WVm.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_24_1_C9c.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_25_1_t2C.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_26_2_vp8.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_27_1_0i4.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_28_1_ToC.root',
                       'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_29_1_N0k.root'
                       ]

TTSource=['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/TTJets_TuneZ2/MergedOutputFile_1_2_wnb.root'
          ]

ZbbSource_SHERPA_5FS=['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_1_3_qPY.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_2_3_qTG.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_3_1_cAj.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_4_2_f3e.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_5_2_DHq.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_6_6_SiH.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_7_8_Ks2.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_8_1_sSR.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_9_2_Uxh.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_10_8_I1u.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_11_1_40S.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_12_8_eQA.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_13_2_L1c.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_14_1_8dW.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_15_7_LwJ.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_16_2_9wq.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_17_2_JXi.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_18_8_nAn.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_19_8_lNa.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_20_7_3VH.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_21_1_Q2q.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_22_2_mVR.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_23_8_dWv.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_24_2_LkU.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_25_7_5xi.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_26_8_yek.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_27_4_Jri.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_28_1_rvj.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_29_4_tFa.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_30_1_EAM.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_31_2_oiI.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_32_1_RyR.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_33_2_815.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_34_2_ZmC.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_35_2_V9P.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_36_7_f1x.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_37_7_qBb.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_38_3_JYU.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_39_1_MP5.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_40_6_Yw8.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_41_7_H4y.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_42_6_iuM.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_43_2_7Ql.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_44_1_AAx.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_45_8_vS0.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_46_7_Sou.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_47_1_LMQ.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_48_1_9zI.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_49_1_RQe.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_50_2_0Vw.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_51_2_8sS.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_52_2_91L.root',
                      'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_53_1_P1X.root'
                      ]

ZbbSource_MADGRAPH_5FS=['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_1_1_6Bb.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_2_1_0Hw.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_3_1_BB4.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_4_1_4ce.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_5_1_MkN.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_6_1_k65.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_7_6_y0u.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_8_1_Go1.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_9_1_5xV.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_10_1_SJL.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_11_1_iAO.root',
                        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_12_4_13F.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_12_7_8Wj.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_13_1_WSG.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_14_2_8Cc.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_15_2_g2n.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_16_1_sQ1.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_17_1_ksn.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_18_1_EBS.root'
                        ]

ZbbSource_MADGRAPH_4FS=['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_MADGRAPH_4FS_preProduction/MergedOutputFile_1_1_dgL.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_MADGRAPH_4FS_preProduction/MergedOutputFile_2_1_2De.root'
                        ]

ZbbSource_aMCatNLO_5FS=['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_10_1_YmL.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_11_1_uaE.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_12_1_pro.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_13_1_8jN.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_14_1_RKD.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_15_1_vwd.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_16_1_vK0.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_17_1_TBj.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_18_1_qZa.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_19_1_noY.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_1_1_SlM.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_20_1_pWu.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_21_1_nrV.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_22_1_vSX.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_23_1_spX.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_24_1_Lab.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_25_1_drT.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_26_1_Az1.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_27_1_ZzD.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_28_1_8pX.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_29_1_YG5.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_2_1_zmu.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_3_1_syz.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_4_1_clu.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_5_1_U0a.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_6_1_ymq.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_7_1_DlA.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_8_1_3cm.root',
                        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_9_1_K5m.root'   
                        ]

ZccSource=DYSource_MADGRAPH_5FS

ZZSource=['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZZ_TuneZ2/MergedOutputFile_1_1_JPk.root',
          'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZZ_TuneZ2/MergedOutputFile_2_1_sgc.root',
          'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZZ_TuneZ2/MergedOutputFile_3_1_wIe.root',
          'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZZ_TuneZ2/MergedOutputFile_4_1_rdh.root',
          ]

WZSource=['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/WZ_TuneZ2/MergedOutputFile_1_1_FiG.root',
          'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/WZ_TuneZ2/MergedOutputFile_2_1_G9x.root' 
          ]

readFiles = cms.untracked.vstring()

process.source = cms.Source("PoolSource",
                            fileNames = readFiles,
                            duplicateCheckMode=cms.untracked.string('checkAllFilesOpened')
                            )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
                             
###################################################################
# Trigger Filter 
###################################################################
### Trigger path FIRED in the event
# import HLTrigger.HLTfilters.hltHighLevel_cfi 
# process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone() 
# process.hltFilter.HLTPaths = [ "HLT_DoubleMu3_v*" ]
# process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
# process.hltFilter.throw  = cms.bool(False)

###################################################################
# HLT run range matcher (for data only)
###################################################################

# from ZbbAnalysis.AnalysisStep.hltrunrangematcher_cfi import hltrunrange
# # muons
# from ZbbAnalysis.AnalysisStep.hltrunrangematcher_mu2010_cff import hltpaths_and_runs as  hltpaths_and_runs_mu2010
# process.mu_hltrunrange = hltrunrange.clone(
#     HLTRunRangeList=hltpaths_and_runs_mu2010
# )
# # electrons
# from ZbbAnalysis.AnalysisStep.hltrunrangematcher_ele2010_cff import hltpaths_and_runs as  hltpaths_and_runs_ele2010
# process.ele_hltrunrange = hltrunrange.clone(
#     HLTRunRangeList=hltpaths_and_runs_ele2010
# )

###################################################################
# Z -> tau tau filter
###################################################################
from ZbbAnalysis.AnalysisStep.ZTauTauFilter_cfi import ztautaufilter

process.zttfilter = ztautaufilter.clone(
    src= cms.InputTag("genParticles")
    ) 

###################################################################
# HF filter (for DY - status3+ status2 daughter veto)
###################################################################

from ZbbAnalysis.AnalysisStep.HeavyFlavorFilter_cfi import heavyflavorfilter

process.b_heavyflavorfilter= heavyflavorfilter.clone(
    src= cms.InputTag("genParticles"),
    status2 = cms.bool(True),
    status3 = cms.bool(True),
    bOn = cms.bool(True),
    cOn = cms.bool(False),
    hDaughterVeto = cms.bool(True),
    zDaughterVeto = cms.bool(False),
    ptcut=cms.double(0)
    )

process.notHadr_b_heavyflavorfilter= heavyflavorfilter.clone(
    src= cms.InputTag("genParticles"),
    status2 = cms.bool(True),
    status3 = cms.bool(True),
    bOn = cms.bool(True),
    cOn = cms.bool(False),
    hDaughterVeto = cms.bool(False),
    zDaughterVeto = cms.bool(False),
    ptcut=cms.double(0)
    )

process.c_heavyflavorfilter= heavyflavorfilter.clone(
    src= cms.InputTag("genParticles"),
    status2 = cms.bool(True),
    status3 = cms.bool(True),
    cOn = cms.bool(True),
    bOn = cms.bool(False),
    hDaughterVeto = cms.bool(True),
    zDaughterVeto = cms.bool(False),
    ptcut=cms.double(0)
    )

process.bc_heavyflavorfilter= heavyflavorfilter.clone(
    src= cms.InputTag("genParticles"),
    status2 = cms.bool(True),
    status3 = cms.bool(True),
    cOn = cms.bool(True),
    bOn = cms.bool(True),
    hDaughterVeto = cms.bool(True),
    zDaughterVeto = cms.bool(False),
    ptcut=cms.double(0)
    )

###################################################################
# PAT basic distributions
###################################################################
from ZbbAnalysis.AnalysisStep.PATAnalyzer_cfi import analyzePAT
process.analyzePat =  analyzePAT.clone(
    jetSrc = cms.untracked.InputTag("patJets")
)

###################################################################
# Jet Cleaning
###################################################################
process.cleanPatJets = cms.EDProducer("PATJetCleaner",
                                      src = cms.InputTag("patJets"),
                                      preselection = cms.string('pt > 20.0 && abs(eta) < 2.4'),
                                      checkOverlaps = cms.PSet(ele = cms.PSet(src       = cms.InputTag("offlinePrimaryVertexFromZ:patElectronsFromZ"), #Z daughters
                                                                              #src       = cms.InputTag("goldenElectrons"),
                                                                              algorithm = cms.string("byDeltaR"),
                                                                              preselection        = cms.string(""),
                                                                              deltaR              = cms.double(0.5),
                                                                              checkRecoComponents = cms.bool(False),
                                                                              pairCut             = cms.string(""),
                                                                              requireNoOverlaps = cms.bool(True)
                                                                              ),
                                                               mu = cms.PSet(src       = cms.InputTag("offlinePrimaryVertexFromZ:patMuonsFromZ"), #Z daughters
                                                                             #src       = cms.InputTag("goldenMuons"),
                                                                             algorithm = cms.string("byDeltaR"),
                                                                             preselection        = cms.string(""),
                                                                             deltaR              = cms.double(0.5),
                                                                             checkRecoComponents = cms.bool(False),
                                                                             pairCut             = cms.string(""),
                                                                             requireNoOverlaps = cms.bool(True)
                                                                             ),
                                                               ),
                                      finalCut = cms.string(''),
                                      )

#process.cleanPatJets     = process.cleanPatJets.clone( src = cms.InputTag("patJets") )
#process.cleanPatJetsNoPU = process.cleanPatJets.clone( src = cms.InputTag("patJetsNoPU") )

###################################################################
# Jet number filter... ?
###################################################################
process.jetFilter = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("cleanPatJets"),
                                 minNumber = cms.uint32(2),
                                 )

process.JetCleaningSeq = cms.Sequence( 
    (
    process.cleanPatJets 
    #+ process.cleanPatJets + 
    #+ process.cleanPatJetsNoPU
     )
    #+ process.jetFilter                  # eventually filter on candidate jets 
    )

###################################################################
# PAT basic distributions after cleaning
###################################################################
process.analyzePatAfterCleaning = analyzePAT.clone(                                           
    jetSrc     = cms.untracked.InputTag("cleanPatJets")                                    
    )

###################################################################
# Event content analysis
###################################################################

###################################################################
# MC Truth Analyzer
###################################################################
from ZbbAnalysis.AnalysisStep.MCTruthAnalyzer_cfi import MCTruthAnalyze

process.MCTruthAnalyzer = MCTruthAnalyze.clone(
    isMCatNLO = cms.untracked.bool(False)
    )

####################################
# SFb from DB                      #
#                                  #
# RecoBTag/PerformanceDB V00-04-11 #
#   REQUIRED IN CMSSW_4_2_1        #
#                                  #
####################################
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")
process.load("RecoBTag.PerformanceDB.BTagPerformanceDB1107")

###################################################################
# Import b-efficiency scaling factors
###################################################################
from ZbbAnalysis.AnalysisStep.effbctables_cff import mc_effbc_ssvhem_2011
from ZbbAnalysis.AnalysisStep.effbctables_cff import mc_effbc_ssvhpt_2011
from ZbbAnalysis.AnalysisStep.effbctables_cff import mc_effbc_csvm_2011
from ZbbAnalysis.AnalysisStep.effbctables_cff import mc_effbc_csvt_2011

# other corrections available ('MISTAGSSVHPT', 'MISTAGSSVHEM','MISTAGCSVM','MISTAGCSVT'),
#                             ('BTAGSSVHPT'  , 'BTAGSSVHEM'  ,'BTAGCSVM'  ,'BTAGCSVT')

from ZbbAnalysis.AnalysisStep.ZbbEventContentAnalyzer_cfi import eventcontentanalyze

### all before tagging histograms
process.steps_before_btag = eventcontentanalyze.clone(
    doAllThePlotting = cms.bool(True),
    isMC        = cms.bool(True),
    bTagAlgoWP  = cms.string('SSVHPT'),
    jetSrc      = cms.untracked.InputTag("cleanPatJets"),
    applybEffCalibration = cms.bool(True),
    bEffCalibrationMethod = cms.string('BTAGSSVHPT'),
    bMistagCalibrationMethod = cms.string('MISTAGSSVHPT'),    
    EffbcMCPtRangeList = mc_effbc_ssvhpt_2011,
    applyPUcorrection = cms.bool(True),
    applyLeptonEfficiency = cms.bool(True),
    minBtags =cms.int32(1)          #minimum b-tags in event        
    )

### ssvhem
process.finaldistros_ssvhem = eventcontentanalyze.clone(
    isMC        = cms.bool(True),
    bTagAlgoWP  = cms.string('SSVHEM'),
    jetSrc      = cms.untracked.InputTag("cleanPatJets"),
    applybEffCalibration = cms.bool(True),
    bEffCalibrationMethod = cms.string('BTAGSSVHEM'),
    bMistagCalibrationMethod = cms.string('MISTAGSSVHEM'),    
    EffbcMCPtRangeList = mc_effbc_ssvhem_2011,
    applyPUcorrection = cms.bool(True),
    applyLeptonEfficiency = cms.bool(True),
    minBtags =cms.int32(1)           #minimum b-tags in event              
    )

### ssvhpt
process.finaldistros_ssvhpt = eventcontentanalyze.clone(
    isMC        = cms.bool(True),
    bTagAlgoWP  = cms.string('SSVHPT'),
    jetSrc      = cms.untracked.InputTag("cleanPatJets"),
    applybEffCalibration = cms.bool(True),
    bEffCalibrationMethod = cms.string('BTAGSSVHPT'),
    bMistagCalibrationMethod = cms.string('MISTAGSSVHPT'),    
    EffbcMCPtRangeList = mc_effbc_ssvhpt_2011,
    applyPUcorrection = cms.bool(True),
    applyLeptonEfficiency = cms.bool(True),
    minBtags =cms.int32(1)          #minimum b-tags in event        
    )

### csvm
process.finaldistros_csvm = eventcontentanalyze.clone(
    isMC        = cms.bool(True),
    bTagAlgoWP  = cms.string('CSVM'),
    jetSrc      = cms.untracked.InputTag("cleanPatJets"),
    applybEffCalibration = cms.bool(True),
    bEffCalibrationMethod = cms.string('BTAGCSVM'),
    bMistagCalibrationMethod = cms.string('MISTAGCSVM'),    
    EffbcMCPtRangeList = mc_effbc_csvm_2011,
    applyPUcorrection = cms.bool(True),
    applyLeptonEfficiency = cms.bool(True),
    minBtags =cms.int32(1)           #minimum b-tags in event              
    )

### csvt
process.finaldistros_csvt = eventcontentanalyze.clone(
    isMC        = cms.bool(True),
    bTagAlgoWP  = cms.string('CSVT'),
    jetSrc      = cms.untracked.InputTag("cleanPatJets"),
    applybEffCalibration = cms.bool(True),
    bEffCalibrationMethod = cms.string('BTAGCSVT'),
    bMistagCalibrationMethod = cms.string('MISTAGCSVT'),    
    EffbcMCPtRangeList = mc_effbc_csvt_2011,
    applyPUcorrection = cms.bool(True),
    applyLeptonEfficiency = cms.bool(True),
    minBtags =cms.int32(1)          #minimum b-tags in event        
    )

##################################################################
# hack to force NLO weights rule
##################################################################
if options.Sample=="Zb5faMCatNLO":
    process.finaldistros_ssvhem.isMCatNLO = cms.bool(True)    
    process.finaldistros_ssvhpt.isMCatNLO = cms.bool(True)    
    process.MCTruthAnalyzer.isMCatNLO = cms.untracked.bool(True)
    
# final sequence
process.finaldistros = cms.Sequence(
    process.steps_before_btag+
    process.finaldistros_ssvhem+
    process.finaldistros_ssvhpt+
    process.finaldistros_csvm+
    process.finaldistros_csvt
    )

###################################################################
# Filter on the generated acceptance
###################################################################
from ZbbAnalysis.AnalysisStep.zb_acceptance_filter_cfi import zb_acceptance_filter 

#filter for b
process.zplusbAccfilter = zb_acceptance_filter.clone(
    flavLept  = cms.untracked.string('all'),   # available Zmm,Zee,all
    etaMuMax  = cms.untracked.double(2.1),
    ptMuMin   = cms.untracked.double(17.),
    etaEleMax = cms.untracked.double(2.5),
    ptEleMin  = cms.untracked.double(17.),                                     
    flavJet   = cms.untracked.string('b'),     # available b,c cases
    etaGenJetMax = cms.untracked.double(3.5),
    ptGenJetMin  = cms.untracked.double(15),
    ngoodGenJet  = cms.untracked.int32(1),
    isExclusive  = cms.untracked.bool(False),  # if True -> == ngoodGenJet, if False -> >= ngoodGenJet
    jetSrc = cms.untracked.InputTag("patJets")
    )

# filter for c
process.zpluscAccfilter = zb_acceptance_filter.clone(
    flavLept  = cms.untracked.string('all'),   # available Zmm,Zee,all
    etaMuMax  = cms.untracked.double(2.1),
    ptMuMin   = cms.untracked.double(17.),
    etaEleMax = cms.untracked.double(2.5),
    ptEleMin  = cms.untracked.double(17.),                                     
    flavJet   = cms.untracked.string('c'),     # available b,c cases   
    etaGenJetMax = cms.untracked.double(3.5),
    ptGenJetMin  = cms.untracked.double(15),
    ngoodGenJet  = cms.untracked.int32(1),
    isExclusive  = cms.untracked.bool(False),  # if True -> == ngoodGenJet, if False -> >= ngoodGenJet
    jetSrc = cms.untracked.InputTag("patJets")
    )

###################################################################
# Definition of analysis sequence
###################################################################
process.AnalysisSequence   = cms.Sequence()

process.keepIfB = cms.Sequence(process.b_heavyflavorfilter*(process.analyzePat+
                                                            process.JetCleaningSeq+
                                                            process.MCTruthAnalyzer+
                                                            process.analyzePatAfterCleaning+
                                                            process.finaldistros
                                                            )
                               )

process.keepIfC = cms.Sequence(~process.b_heavyflavorfilter*process.c_heavyflavorfilter*(process.analyzePat+
                                                                                         process.JetCleaningSeq+
                                                                                         process.MCTruthAnalyzer+
                                                                                         process.analyzePatAfterCleaning+
                                                                                         process.finaldistros
                                                                                         )
                               )

process.dropIfBC = cms.Sequence(~process.bc_heavyflavorfilter*(process.analyzePat+
                                                               process.JetCleaningSeq+
                                                               process.MCTruthAnalyzer+
                                                               process.analyzePatAfterCleaning+
                                                               process.finaldistros
                                                               )
                                )

process.keepIfNotHadrB = cms.Sequence(process.notHadr_b_heavyflavorfilter*(process.analyzePat+
                                                                           process.JetCleaningSeq+
                                                                           process.MCTruthAnalyzer+
                                                                           process.analyzePatAfterCleaning+
                                                                           process.finaldistros
                                                                           )
                                      )

## sequence without HeavyFlavourFilter 
process.AnalysisNoFilter = cms.Sequence(process.analyzePat+
                                        process.JetCleaningSeq+
                                        process.MCTruthAnalyzer+
                                        process.analyzePatAfterCleaning+
                                        process.finaldistros
                                        )

##############################################################################
# Patch for sum(E) violation in PYTHIA/MADGRAPH DY samples
# https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/1489.html
##############################################################################
process.load("ZbbAnalysis.Tools.TotalKinematicsFilter_cfi")

if options.Sample=="Zl" or options.Sample=="Zb5f" or options.Sample=="Zc":
    process.AnalysisSequence+=process.totalKinematicsFilter

###################################################################
# switches
###################################################################
if options.Sample=="Zl":    
    readFiles.extend(DYSource_MADGRAPH_5FS)
    outputRootFileName=cms.string('analyzePAT_MC_ZlJets.root')
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.dropIfBC
elif options.Sample=="Zb5f":
    readFiles.extend(ZbbSource_MADGRAPH_5FS)
    outputRootFileName=cms.string('analyzePAT_MC_Zb5fToLL.root')
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.keepIfB
elif options.Sample=="Zb4f":
    print "============================================="
    print "%msg-w: !!!  Beware using bugged sample  !!!!"
    print "============================================="
    readFiles.extend(ZbbSource_MADGRAPH_4FS)
    outputRootFileName=cms.string('analyzePAT_MC_Zb4fToLL.root')
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.keepIfB    
elif options.Sample=="Zc":
    readFiles.extend(ZccSource)
    outputRootFileName=cms.string('analyzePAT_MC_ZcToLL.root')
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.keepIfC
elif options.Sample=="tt":
    readFiles.extend(TTSource)
    outputRootFileName=cms.string('analyzePAT_MC_TTJets.root')
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.AnalysisNoFilter
elif options.Sample=="Zb5fSherpa":
    readFiles.extend(ZbbSource_SHERPA_5FS)
    outputRootFileName=cms.string('analyzePAT_MC_Zb5f_Sherpa.root')
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.keepIfNotHadrB
elif options.Sample=="Zb5faMCatNLO":
    readFiles.extend(ZbbSource_aMCatNLO_5FS)
    outputRootFileName=cms.string('analyzePAT_MC_Zb5f_aMCatNLO.root')
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.AnalysisNoFilter
elif options.Sample=="zz":
    readFiles.extend(ZZSource)
    outputRootFileName=cms.string('analyzePAT_MC_ZZ.root')
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.AnalysisNoFilter
elif options.Sample=="wz":
    readFiles.extend(WZSource)
    outputRootFileName=cms.string('analyzePAT_MC_WZ.root')
    process.AnalysisSequence+=~process.zttfilter
    process.AnalysisSequence+=process.AnalysisNoFilter
elif options.Sample=="Ztautau":
    readFiles.extend(DYSource_MADGRAPH_5FS)
    outputRootFileName=cms.string('analyzePAT_MC_ZtautauJets.root')
    process.AnalysisSequence+=process.zttfilter
    process.AnalysisSequence+=process.AnalysisNoFilter
else :
    print "============================================="
    print "%msg-w: error while exectuing ZbbEventContent"
    print "        Unrecognized sample"
    print "        Please choose among: Zl,Zb5f,Zb4f,Zb5fSherpa"
    print "                             Zb5faMCatNLO,Zc or tt"
    print "============================================="
    outputRootFileName=cms.string('fake.root')

###################################################################
# output file
###################################################################
process.TFileService = cms.Service("TFileService", fileName = outputRootFileName)

###################################################################
# path
###################################################################
process.p = cms.Path(process.AnalysisSequence)



