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
                 "Sample to be processed (Zl,Zb4f,Zb5f,Zc,tt,DY)")
options.register('maxEvents',
                 -1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to process (-1 for all)")
options.register('redoTagging',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "True: re-reun btagging; False: do not re-run btagging")
options.register('isMC',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "True: is a MC sample; False: is a data sample")
options.register('isMuon',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "True: run on muons; False: run on electrons")
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
if options.isMC :
    process.GlobalTag.globaltag = 'START42_V13::All'  # for MC
else :
    process.GlobalTag.globaltag = 'GR_R_42_V20::All'  # for the DATA
    
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

DYSource=['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_1_1_lot.root',
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

Zbb5fSource=['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_1_1_6Bb.root',
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

ZccSource=DYSource

TTSource=['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/TTJets_TuneZ2/MergedOutputFile_1_2_wnb.root']

Zbb4fSource=['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_MADGRAPH_4FS_preProduction/MergedOutputFile_1_1_dgL.root',
             'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_MADGRAPH_4FS_preProduction/MergedOutputFile_2_1_2De.root'
             ]

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
# Jet number filter... (?)
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
# Define btagging sequence
###################################################################
process.btagging = cms.Sequence()

if options.redoTagging:
    
    # right calibrations for JP
    if options.isMC :
        process.GlobalTag.toGet = cms.VPSet(
            cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
                     tag = cms.string("TrackProbabilityCalibration_2D_2011_v1_mc"),
                     connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU")),
            cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
                     tag = cms.string("TrackProbabilityCalibration_3D_2011_v1_mc"),
                     connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU"))
            )
    else :   
        process.GlobalTag.toGet = cms.VPSet(
            cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
                     tag = cms.string("TrackProbabilityCalibration_2D_2011Data_v1_offline"),
                     connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU")),
            cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),    
                     tag = cms.string("TrackProbabilityCalibration_3D_2011Data_v1_offline"),
                     connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU"))
            )

    #define the b-tag squences for offline reconstruction
    process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
    process.MyJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
                                                           process.j2tParametersVX,
                                                           jets = cms.InputTag("cleanPatJets")
                                                           )
    
    process.load("RecoBTag.ImpactParameter.impactParameter_cff")
    process.MyImpactParameterTagInfos = process.impactParameterTagInfos.clone(
        jetTracks = "MyJetTracksAssociatorAtVertex"
        )
    
    # re-run track counting b-tagging algorithm
    process.MyTrackCountingHighEffBJetTags = process.trackCountingHighEffBJetTags.clone(
        tagInfos = cms.VInputTag(cms.InputTag("MyImpactParameterTagInfos"))
        )
    
    process.MyTrackCountingHighPurBJetTags = process.trackCountingHighPurBJetTags.clone(
        tagInfos = cms.VInputTag(cms.InputTag("MyImpactParameterTagInfos"))
        )
    
    process.MyJetProbabilityBJetTags = process.jetProbabilityBJetTags.clone(
        tagInfos = cms.VInputTag(cms.InputTag("MyImpactParameterTagInfos"))
        )
    
    process.btagging += cms.Sequence(process.MyJetTracksAssociatorAtVertex*
                                     (process.MyImpactParameterTagInfos *
                                      (process.MyTrackCountingHighEffBJetTags +
                                       process.MyTrackCountingHighPurBJetTags +
                                       process.MyJetProbabilityBJetTags 
                                       )
                                      )
                                     )

###################################################################
# Event content analysis
###################################################################
from ZbbAnalysis.AnalysisStep.PatBJetTagAnalyzer_cfi import analyzeBJets
process.bjetanalysis_SSVHEM =  analyzeBJets.clone(primaryVertices = cms.InputTag("offlinePrimaryVertexFromZ:ZVertexProducerZVertex:PAT"),
                                                  bTagAlgoWP  = cms.string('SSVHEM'), 
                                                  tracks = cms.InputTag("generalTracks"),
                                                  jets = cms.InputTag("cleanPatJets"),
                                                  isMC = cms.bool(options.isMC),
                                                  redoBTag = cms.bool(options.redoTagging)
                                                  )

process.bjetanalysis_SSVHPT =  analyzeBJets.clone(primaryVertices = cms.InputTag("offlinePrimaryVertexFromZ:ZVertexProducerZVertex:PAT"),
                                                  bTagAlgoWP  = cms.string('SSVHPT'), 
                                                  tracks = cms.InputTag("generalTracks"),
                                                  jets = cms.InputTag("cleanPatJets"),
                                                  isMC = cms.bool(options.isMC),
                                                  redoBTag = cms.bool(options.redoTagging)
                                                  )

process.bjetanalysis_CSVM   =  analyzeBJets.clone(primaryVertices = cms.InputTag("offlinePrimaryVertexFromZ:ZVertexProducerZVertex:PAT"),
                                                  bTagAlgoWP  = cms.string('CSVM'), 
                                                  tracks = cms.InputTag("generalTracks"),
                                                  jets = cms.InputTag("cleanPatJets"),
                                                  isMC = cms.bool(options.isMC),
                                                  redoBTag = cms.bool(options.redoTagging)
                                                  )

process.bjetanalysis_CSVT   =  analyzeBJets.clone(primaryVertices = cms.InputTag("offlinePrimaryVertexFromZ:ZVertexProducerZVertex:PAT"),
                                                  bTagAlgoWP  = cms.string('CSVT'), 
                                                  tracks = cms.InputTag("generalTracks"),
                                                  jets = cms.InputTag("cleanPatJets"),
                                                  isMC = cms.bool(options.isMC),
                                                  redoBTag = cms.bool(options.redoTagging)
                                                  )

process.bjetanalysis = cms.Sequence(
    process.bjetanalysis_SSVHEM+
    process.bjetanalysis_SSVHPT+
    process.bjetanalysis_CSVM+
    process.bjetanalysis_CSVT
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
    jetSrc = cms.untracked.InputTag("cleanPatJets")
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
    jetSrc = cms.untracked.InputTag("cleanPatJets")
    )

###################################################################
# Definition of analysis sequence
###################################################################
process.AnalysisSequence   = cms.Sequence()

process.keepIfB = cms.Sequence(process.b_heavyflavorfilter*(process.JetCleaningSeq+
                                                            process.btagging*
                                                            process.bjetanalysis
                                                            )
                               )

process.keepIfC = cms.Sequence(~process.b_heavyflavorfilter*process.c_heavyflavorfilter*(process.JetCleaningSeq+
                                                                                         process.btagging*
                                                                                         process.bjetanalysis
                                                                                         )
                               )

process.dropIfBC = cms.Sequence(~process.bc_heavyflavorfilter*(process.JetCleaningSeq+
                                                               process.btagging*
                                                               process.bjetanalysis
                                                               )
                                )
## sequence without hf 
process.AnalysisNoFilter = cms.Sequence(process.JetCleaningSeq+
                                        process.btagging*
                                        process.bjetanalysis
                                        )

###################################################################
# switches
###################################################################
if options.isMC :
    if options.Sample=="Zl":    
        readFiles.extend(DYSource)
        outputRootFileName=cms.string('testBTagJet_MC_ZlJets.root')
    elif options.Sample=="DY":    
        readFiles.extend(DYSource)
        outputRootFileName=cms.string('testBTagJet_MC_DYJets.root')
        process.AnalysisSequence+=process.AnalysisNoFilter
    elif options.Sample=="Zb5f":
        readFiles.extend(Zbb5fSource)
        outputRootFileName=cms.string('testBTagJet_MC_Zb5fToLL.root')
        process.AnalysisSequence+=process.keepIfB
    elif options.Sample=="Zc":
        readFiles.extend(ZccSource)
        outputRootFileName=cms.string('testBTagJet_MC_ZcToLL.root')
        process.AnalysisSequence+=process.keepIfC
    elif options.Sample=="tt":
        readFiles.extend(TTSource)
        outputRootFileName=cms.string('testBTagJet_MC_TTJets.root')
        process.AnalysisSequence+=process.AnalysisNoFilter
    elif options.Sample=="Zb4f":
        readFiles.extend(Zbb4fSource)
        outputRootFileName=cms.string('testBTagJet_MC_Zb4fJets.root')
        process.AnalysisSequence+=process.AnalysisNoFilter
    else :
        print "============================================="
        print "%msg-w: error while exectuing ZbbEventContent"
        print "        Unrecognized sample"
        print "        Please choose among DY,Zl,Zb4f,Zb5f,Zc or tt"
        print "============================================="
        outputRootFileName=cms.string('fake.root')
else :
    
    if options.isMuon:    
        readFiles.extend(muonSource)
        process.AnalysisSequence+=process.AnalysisNoFilter
        outputRootFileName=cms.string('testBTagJet_DATA2011_muon.root')
    else:
        readFiles.extend(electronSource)
        process.AnalysisSequence+=process.AnalysisNoFilter
        outputRootFileName=cms.string('testBTagJet_DATA2011_electron.root')
        
###################################################################
# output file
###################################################################

process.TFileService = cms.Service("TFileService", fileName = outputRootFileName)

###################################################################
# path
###################################################################

# process.Out = cms.OutputModule("PoolOutputModule")
# process.Out.fileName = cms.untracked.string('file:/tmp/musich/myzbbPATSkim_42X_withBTag.root')
# process.Out.outputCommands =  cms.untracked.vstring(
#    'keep *',
#    )

process.p = cms.Path(process.AnalysisSequence)
#process.end = cms.EndPath(process.Out)


