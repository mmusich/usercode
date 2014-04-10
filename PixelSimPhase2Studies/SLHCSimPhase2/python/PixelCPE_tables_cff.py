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
    xerr_barrel_ln_ = cms.untracked.vdouble(0.0013,0.000721,0.00179),
    xerr_barrel_ln_def_ = cms.untracked.double(0.00457),
    yerr_barrel_ln_ = cms.untracked.vdouble(0.0124,0.00216,0.00191,0.0018,0.00166,0.00155,0.00197,0.00201,0.00898),
    yerr_barrel_ln_def_ = cms.untracked.double(0.0183),
    xerr_endcap_ = cms.untracked.vdouble(0.003,0.000566,0.00123),
    xerr_endcap_def_ = cms.untracked.double(0.00121),
    yerr_endcap_ = cms.untracked.vdouble(0.00497,0.000983,0.00713),
    yerr_endcap_def_ = cms.untracked.double(0.0049),
    xerr_barrel_l1_ = cms.untracked.vdouble(0.00157,0.000792,0.00137),
    xerr_barrel_l1_def_ = cms.untracked.double(0.00372),
    yerr_barrel_l1_ = cms.untracked.vdouble(0.0032,0.0017,0.00155,0.00165,0.00198,0.00175,0.00177,0.00161,0.00186),
    yerr_barrel_l1_def_ = cms.untracked.double(0.00189)
    )


# phase2v0 100x150 um2  285/285 um 1200e NoAging
pixelCPE_phase2v0 = cms.PSet(
    xerr_barrel_ln_ = cms.untracked.vdouble(0.00121,0.000691,0.00199),
    xerr_barrel_ln_def_ = cms.untracked.double(0.00468),
    yerr_barrel_ln_ = cms.untracked.vdouble(0.0129,0.00217,0.00171,0.00172,0.00185,0.00157,0.00208,0.00214,0.00214),
    yerr_barrel_ln_def_ = cms.untracked.double(0.0164),
    xerr_endcap_ = cms.untracked.vdouble(0.00302,0.000577,0.00125),
    xerr_endcap_def_ = cms.untracked.double(0.0016),
    yerr_endcap_ = cms.untracked.vdouble(0.00495,0.000982,0.00713),
    yerr_endcap_def_ = cms.untracked.double(0.00539),
    xerr_barrel_l1_ = cms.untracked.vdouble(0.00158,0.000714,0.00154),
    xerr_barrel_l1_def_ = cms.untracked.double(0.0064),
    yerr_barrel_l1_ = cms.untracked.vdouble(0.00305,0.00167,0.00181,0.00174,0.00183,0.00185,0.00203,0.00155,0.00198),
    yerr_barrel_l1_def_ = cms.untracked.double(0.00166)
)

# phase2v1 100x150 um2  285/285 um 1200e 300/fb
pixelCPE_phase2v1 = cms.PSet(
    xerr_barrel_ln_ = cms.untracked.vdouble(0.00121,0.0007,0.00194),
    xerr_barrel_ln_def_ = cms.untracked.double(0.00484),
    yerr_barrel_ln_ = cms.untracked.vdouble(0.0128,0.0022,0.00174,0.00199,0.00234,0.0017,0.00223,0.00297,0.0173),
    yerr_barrel_ln_def_ = cms.untracked.double(0.0324),
    xerr_endcap_ = cms.untracked.vdouble(0.00301,0.000647,0.00133),
    xerr_endcap_def_ = cms.untracked.double(0.00192),
    yerr_endcap_ = cms.untracked.vdouble(0.00484,0.00111,0.0138),
    yerr_endcap_def_ = cms.untracked.double(0.00399),
    xerr_barrel_l1_ = cms.untracked.vdouble(0.00164,0.0009,0.00136),
    xerr_barrel_l1_def_ = cms.untracked.double(0.00522),
    yerr_barrel_l1_ = cms.untracked.vdouble(0.00334,0.00237,0.00323,0.00382,0.00364,0.00409,0.00606,0.00723,0.00572),
    yerr_barrel_l1_def_ = cms.untracked.double(0.00608)
)

# phase2v2  50x100 um2  150/285 um 1200e 300/fb
pixelCPE_phase2v2 = cms.PSet(
    xerr_barrel_ln_ = cms.untracked.vdouble(0.00122,0.000716,0.00206),
    xerr_barrel_ln_def_ = cms.untracked.double(0.00627),
    yerr_barrel_ln_ = cms.untracked.vdouble(0.0123,0.00214,0.00189,0.00198,0.00223,0.00242,0.00176,0.00306,0.00378),
    yerr_barrel_ln_def_ = cms.untracked.double(0.0066),
    xerr_endcap_ = cms.untracked.vdouble(0.00229,0.000614,0.00173),
    xerr_endcap_def_ = cms.untracked.double(0.00123),
    yerr_endcap_ = cms.untracked.vdouble(0.00463,0.00108,0.00823),
    yerr_endcap_def_ = cms.untracked.double(0.0175),
    xerr_barrel_l1_ = cms.untracked.vdouble(0.00089,0.000601,0.000665),
    xerr_barrel_l1_def_ = cms.untracked.double(0.00281),
    yerr_barrel_l1_ = cms.untracked.vdouble(0.00229,0.0019,0.00261,0.00292,0.00302,0.00404,0.00385,0.00518,0.00337),
    yerr_barrel_l1_def_ = cms.untracked.double(0.00667)
)

# phase2v3  25x100 um2  100/285 um 1200e 300/fb
pixelCPE_phase2v3 = cms.PSet(
    xerr_barrel_ln_ = cms.untracked.vdouble(0.00121,0.000714,0.002),
    xerr_barrel_ln_def_ = cms.untracked.double(0.00495),
    yerr_barrel_ln_ = cms.untracked.vdouble(0.012,0.00209,0.00185,0.00183,0.00201,0.00208,0.00226,0.00416,0.00909),
    yerr_barrel_ln_def_ = cms.untracked.double(0.00319),
    xerr_endcap_ = cms.untracked.vdouble(0.00223,0.000586,0.00135),
    xerr_endcap_def_ = cms.untracked.double(0.0012),
    yerr_endcap_ = cms.untracked.vdouble(0.00439,0.00105,0.00466),
    yerr_endcap_def_ = cms.untracked.double(0.0193),
    xerr_barrel_l1_ = cms.untracked.vdouble(0.000467,0.000461,0.000506),
    xerr_barrel_l1_def_ = cms.untracked.double(0.00198),
    yerr_barrel_l1_ = cms.untracked.vdouble(0.00236,0.00208,0.00244,0.00259,0.004,0.00354,0.00264,0.224,0.00472),
    yerr_barrel_l1_def_ = cms.untracked.double(0.0454)
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
                  'pixelCPE_phase2v0'        : pixelCPE_phase2v0       ,
                  'pixelCPE_phase2v1'        : pixelCPE_phase2v1       ,
                  'pixelCPE_phase2v2'        : pixelCPE_phase2v2       ,
                  'pixelCPE_phase2v3'        : pixelCPE_phase2v3       ,
                  'pixelCPE_dummy'           : pixelCPE_dummy
                  }
