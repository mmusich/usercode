import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("TEST")

options = VarParsing.VarParsing()
options.register('Sample',
                 "MadGraph", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string,         # string, int, or float
                 "Sample to be processed (MadGraph,Sherpa,aMCatNLO)")
options.register('maxEvents',
                 -1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to process (-1 for all)")
options.register('isNLO',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "is the MC NLO? (True,False)")
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

SherpaSrc = ['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_1_3_qPY.root',
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

MadGraphSrc = ['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_1_1_6Bb.root',
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

aMCatNLOSrc = ['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_10_1_YmL.root',
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

readFiles = cms.untracked.vstring()

process.source = cms.Source("PoolSource",
                            fileNames = readFiles
                            #firstEvent = cms.untracked.uint32(6391783),
                            #firstRun = cms.untracked.uint32(1)
                            )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

###############################################################
#Producer: create new collections of leptons from HF semileptonic decays
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
process.MCTruthAnalyzer = cms.EDAnalyzer("MCTruthAnalyzer",
                                         JetCollection = cms.InputTag("cleanPatJets"),
                                         genJetSrc = cms.InputTag("patJets:genJets"),
                                         src = cms.InputTag("genParticles"),
                                         ElectronCollection = cms.InputTag("patElectrons"),
                                         MuonCollection = cms.InputTag("patMuons"),
                                         ZmmCollection = cms.InputTag('zMMCand'),
                                         isMCatNLO = cms.untracked.bool(options.isNLO)
                                         )

###################################################################
# Definition of analysis sequence
###################################################################
process.AnalysisSequence   = cms.Sequence()

process.mcTruthAnalyze = cms.Sequence(process.b_heavyflavorfilter*
                                     process.JetCleaningSeq*
                                     process.MCTruthAnalyzer
                                     )                              

###################################################################
# MC Total kinematics filter
###################################################################
process.load("ZbbAnalysis.Tools.TotalKinematicsFilter_cfi")

if options.Sample=="MadGraph":
    process.AnalysisSequence+=process.totalKinematicsFilter

###################################################################
# switches
###################################################################
if options.Sample=="MadGraph":    
    readFiles.extend(MadGraphSrc)
    process.AnalysisSequence+=process.mcTruthAnalyze
    outputRootFileName=cms.string('MCTruth_ZbbToLL_MADGRAPH.root')
elif options.Sample=="Sherpa":
    readFiles.extend(SherpaSrc)
    process.AnalysisSequence+=process.mcTruthAnalyze
    outputRootFileName=cms.string('MCTruth_ZbbToLL_SHERPA.root')
elif options.Sample=="aMCatNLO":
    readFiles.extend(aMCatNLOSrc)
    process.AnalysisSequence+=process.mcTruthAnalyze
    outputRootFileName=cms.string('MCTruth_ZbbToLL_aMCatNLO.root')

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
                                 Inputfile=cms.untracked.string("list.txt")
                                 )

process.p = cms.Path(process.AnalysisSequence) 
#process.p = cms.Path(process.checkList*process.printTree)
#process.p = cms.Path(process.b_heavyflavorfilter)
#process.e = cms.EndPath(process.out)




