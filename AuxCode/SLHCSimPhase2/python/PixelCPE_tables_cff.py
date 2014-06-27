import FWCore.ParameterSet.Config as cms

# taken from https://github.com/cms-sw/cmssw/blob/CMSSW_6_2_X_SLHC/RecoLocalTracker/SiPixelRecHits/src/PixelCPEGeneric.cc

#Pixel CPE for default pixel cell (100um x 150um)
pixelCPE_100x150_default = cms.PSet(
    xerr_barrel_l1_= cms.untracked.vdouble(0.00115, 0.00120, 0.00088),
    xerr_barrel_l1_def_=cms.untracked.double(0.01030),
    yerr_barrel_l1_= cms.untracked.vdouble(0.00375,0.00230,0.00250,0.00250,0.00230,0.00230,0.00210,0.00210,0.00240),
    yerr_barrel_l1_def_=cms.untracked.double(0.00210),
    xerr_barrel_ln_= cms.untracked.vdouble(0.00115, 0.00120, 0.00088),
    xerr_barrel_ln_def_=cms.untracked.double(0.01030),
    yerr_barrel_ln_= cms.untracked.vdouble(0.00375,0.00230,0.00250,0.00250,0.00230,0.00230,0.00210,0.00210,0.00240),
    yerr_barrel_ln_def_=cms.untracked.double(0.00210),
    xerr_endcap_= cms.untracked.vdouble(0.0020, 0.0020),
    xerr_endcap_def_=cms.untracked.double(0.0020),
    yerr_endcap_= cms.untracked.vdouble(0.00210),
    yerr_endcap_def_=cms.untracked.double(0.00075)
    )

#Pixel CPE for upgrade pixel cell (100um x 150um) 
pixelCPE_100x150_upgrade = cms.PSet(
    xerr_barrel_l1_= cms.untracked.vdouble(0.00114,0.00104,0.00214),
    xerr_barrel_l1_def_= cms.untracked.double(0.00425),
    yerr_barrel_l1_= cms.untracked.vdouble(0.00299,0.00203,0.0023,0.00237,0.00233,0.00243,0.00232,0.00259,0.00176),
    yerr_barrel_l1_def_=cms.untracked.double(0.00245),
    xerr_barrel_ln_= cms.untracked.vdouble(0.00114,0.00104,0.00214),
    xerr_barrel_ln_def_=cms.untracked.double(0.00425),
    yerr_barrel_ln_= cms.untracked.vdouble(0.00299,0.00203,0.0023,0.00237,0.00233,0.00243,0.00232,0.00259,0.00176),
    yerr_barrel_ln_def_=cms.untracked.double(0.00245),
    xerr_endcap_= cms.untracked.vdouble(0.00151,0.000813,0.00221),
    xerr_endcap_def_=cms.untracked.double(0.00218),
    yerr_endcap_= cms.untracked.vdouble(0.00261,0.00107,0.00264),
    yerr_endcap_def_=cms.untracked.double(0.00357)
    )


# phase1 100x150 um2  285/285 um 2000e NoAging
pixelCPE_phase1 = cms.PSet(
    xerr_barrel_l1_ = cms.untracked.vdouble(0.00153,0.000721,0.00142),
    xerr_barrel_l1_def_ = cms.untracked.double(0.013),
    yerr_barrel_l1_ = cms.untracked.vdouble(0.00289,0.00122,0.00181,0.00184,0.00189,0.00188,0.0019,0.00189,0.00191),
    yerr_barrel_l1_def_ = cms.untracked.double(0.00192),
    xerr_barrel_ln_ = cms.untracked.vdouble(0.000732,0.000809,0.00177),
    xerr_barrel_ln_def_ = cms.untracked.double(0.0101),
    yerr_barrel_ln_ = cms.untracked.vdouble(0.00288,0.00121,0.0018,0.00183,0.00183,0.00189,0.00187,0.00171,0.035),
    yerr_barrel_ln_def_ = cms.untracked.double(0.134),
    xerr_endcap_ = cms.untracked.vdouble(0.00131,0.000629,0.00423),
    xerr_endcap_def_ = cms.untracked.double(0.00965),
    yerr_endcap_ = cms.untracked.vdouble(0.00259,0.000809,0.00663),
    yerr_endcap_def_ = cms.untracked.double(0.0229)
    )

# phase1 100x150 um2  220/285 um 1000e L1 NoAging
pixelCPE_phase1v2 = cms.PSet(
    xerr_barrel_l1_ = cms.untracked.vdouble(0.00111,0.000555,0.00109),
    xerr_barrel_l1_def_ = cms.untracked.double(0.00875),
    yerr_barrel_l1_ = cms.untracked.vdouble(0.00187,0.000816,0.00127,0.00131,0.00132,0.00136,0.00135,0.00136,0.00136),
    yerr_barrel_l1_def_ = cms.untracked.double(0.00137),
    xerr_barrel_ln_ = cms.untracked.vdouble(0.000628,0.000768,0.00187),
    xerr_barrel_ln_def_ = cms.untracked.double(0.00962),
    yerr_barrel_ln_ = cms.untracked.vdouble(0.00284,0.00106,0.00179,0.00182,0.00178,0.00183,0.00178,0.00179,0.00207),
    yerr_barrel_ln_def_ = cms.untracked.double(0.137),
    xerr_endcap_ = cms.untracked.vdouble(0.00131,0.000629,0.00408),
    xerr_endcap_def_ = cms.untracked.double(0.01),
    yerr_endcap_ = cms.untracked.vdouble(0.00259,0.000808,0.00676),
    yerr_endcap_def_ = cms.untracked.double(0.0241)
    )


# phase1 100x150 um2  285/285 um 2000e Aging 300/fb
pixelCPE_phase1_300fb = cms.PSet(
    xerr_barrel_l1_ = cms.untracked.vdouble(0.00165,0.00093,0.00135),
    xerr_barrel_l1_def_ = cms.untracked.double(0.0153),
    yerr_barrel_l1_ = cms.untracked.vdouble(0.00296,0.00247,0.00355,0.00403,0.00439,0.00451,0.0046,0.00466,0.00475),
    yerr_barrel_l1_def_ = cms.untracked.double(0.00481),
    xerr_barrel_ln_ = cms.untracked.vdouble(0.000788,0.000841,0.00173),
    xerr_barrel_ln_def_ = cms.untracked.double(0.0101),
    yerr_barrel_ln_ = cms.untracked.vdouble(0.00289,0.00132,0.00197,0.00209,0.00219,0.00241,0.00241,0.00213,0.109),
    yerr_barrel_ln_def_ = cms.untracked.double(0.0503),
    xerr_endcap_ = cms.untracked.vdouble(0.00139,0.000664,0.00412),
    xerr_endcap_def_ = cms.untracked.double(0.00863),
    yerr_endcap_ = cms.untracked.vdouble(0.00263,0.00098,0.00712),
    yerr_endcap_def_ = cms.untracked.double(0.00269)
)

# phase1 100x150 um2  285/285 um 2000e Ageing 500/fb
pixelCPE_phase1_500fb = cms.PSet(
    xerr_barrel_l1_ = cms.untracked.vdouble(0.00179,0.00111,0.00137),
    xerr_barrel_l1_def_ = cms.untracked.double(0.0181),
    yerr_barrel_l1_ = cms.untracked.vdouble(0.00323,0.00369,0.00536,0.0062,0.00713,0.00739,0.00765,0.00803,0.00851),
    yerr_barrel_l1_def_ = cms.untracked.double(0.00827),
    xerr_barrel_ln_ = cms.untracked.vdouble(0.00083,0.000869,0.00172),
    xerr_barrel_ln_def_ = cms.untracked.double(0.0102),
    yerr_barrel_ln_ = cms.untracked.vdouble(0.00289,0.00143,0.00215,0.00239,0.0026,0.00297,0.00298,0.00253,0.115),
    yerr_barrel_ln_def_ = cms.untracked.double(0.126),
    xerr_endcap_ = cms.untracked.vdouble(0.00144,0.00069,0.00441),
    xerr_endcap_def_ = cms.untracked.double(0.00757),
    yerr_endcap_ = cms.untracked.vdouble(0.00266,0.00114,0.00766),
    yerr_endcap_def_ = cms.untracked.double(0.0407)
)



#Dummy Pixel CPE (for test purposes)  
pixelCPE_dummy= cms.PSet(
    xerr_barrel_l1_=cms.untracked.vdouble(0.00415,0.0042,0.00388),
    xerr_barrel_l1_def_=cms.untracked.double(0.04030),
    yerr_barrel_l1_=cms.untracked.vdouble(0.01375,0.00830,0.00850,0.00850,0.00830,0.00830,0.00810,0.00810,0.00840),
    yerr_barrel_l1_def_=cms.untracked.double(0.00810),
    xerr_barrel_ln_=cms.untracked.vdouble(0.00415,0.00420,0.00388),
    xerr_barrel_ln_def_=cms.untracked.double(0.04030),
    yerr_barrel_ln_=cms.untracked.vdouble(0.01375,0.00830,0.00850,0.00850,0.00830,0.00830,0.00810,0.00810,0.00840),
    yerr_barrel_ln_def_=cms.untracked.double(0.00810),
    xerr_endcap_=cms.untracked.vdouble(0.0080,0.0080),
    xerr_endcap_def_=cms.untracked.double(0.0080),
    yerr_endcap_=cms.untracked.vdouble(0.00810),
    yerr_endcap_def_=cms.untracked.double(0.00275)
    )

PixelCPE_dict = { 'pixelCPE_100x150_default' : pixelCPE_100x150_default,
                  'pixelCPE_100x150_upgrade' : pixelCPE_100x150_upgrade,
                  'pixelCPE_phase1'          : pixelCPE_phase1         ,
                  'pixelCPE_phase1v2'        : pixelCPE_phase1v2       ,
                  'pixelCPE_phase1_300fb'    : pixelCPE_phase1_300fb   ,
                  'pixelCPE_phase1_500fb'    : pixelCPE_phase1_500fb   ,   
                  'pixelCPE_dummy'           : pixelCPE_dummy
                  }
