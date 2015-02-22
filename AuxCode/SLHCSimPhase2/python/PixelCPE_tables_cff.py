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


# phase1 100x150 um2  285/285 um 2000e NoAgeing
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

# phase1 100x150 um2  220/285 um 1000e L1 NoAgeing
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


### study of 300/fb and 500/fb ageing
# phase1 100x150 um2  285/285 um 2000e No Ageing
pixelCPE_phase1_0fb = cms.PSet(
    xerr_barrel_l1_ = cms.untracked.vdouble(0.00123,0.000681,0.00142),
    xerr_barrel_l1_def_ = cms.untracked.double(0.0134),
    yerr_barrel_l1_ = cms.untracked.vdouble(0.00234,0.00115,0.00174,0.00179,0.00183,0.00183,0.00185,0.00183,0.00185),
    yerr_barrel_l1_def_ = cms.untracked.double(0.00185),
    xerr_barrel_ln_ = cms.untracked.vdouble(0.000639,0.000763,0.00179),
    xerr_barrel_ln_def_ = cms.untracked.double(0.0104),
    yerr_barrel_ln_ = cms.untracked.vdouble(0.00263,0.00114,0.00175,0.00176,0.00177,0.00183,0.00183,0.00166,0.111),
    yerr_barrel_ln_def_ = cms.untracked.double(0.0738),
    xerr_endcap_ = cms.untracked.vdouble(0.00129,0.000591,0.00421),
    xerr_endcap_def_ = cms.untracked.double(0.0102),
    yerr_endcap_ = cms.untracked.vdouble(0.00247,0.000765,0.00674),
    yerr_endcap_def_ = cms.untracked.double(0.0469)
)

# phase1 100x150 um2  285/285 um 2000e Ageing CMSSW 300/fb
pixelCPE_phase1_300fb = cms.PSet(
    xerr_barrel_ln_ = cms.untracked.vdouble(0.000742,0.000836,0.00174),
    xerr_barrel_ln_def_ = cms.untracked.double(0.0102),
    yerr_barrel_ln_ = cms.untracked.vdouble(0.00272,0.0013,0.00197,0.00208,0.00219,0.0024,0.00241,0.00213,0.113),
    yerr_barrel_ln_def_ = cms.untracked.double(0.0531),
    xerr_endcap_ = cms.untracked.vdouble(0.00139,0.000664,0.00412),
    xerr_endcap_def_ = cms.untracked.double(0.00863),
    yerr_endcap_ = cms.untracked.vdouble(0.00263,0.00098,0.00712),
    yerr_endcap_def_ = cms.untracked.double(0.00269),
    xerr_barrel_l1_ = cms.untracked.vdouble(0.0015,0.000913,0.00135),
    xerr_barrel_l1_def_ = cms.untracked.double(0.0155),
    yerr_barrel_l1_ = cms.untracked.vdouble(0.00275,0.00247,0.00355,0.00402,0.00445,0.00453,0.00461,0.00467,0.00476),
    yerr_barrel_l1_def_ = cms.untracked.double(0.00481)
    )

# phase1 100x150 um2  285/285 um 2000e Ageing CMSSW 500/fb
pixelCPE_phase1_500fb = cms.PSet(
    xerr_barrel_ln_ = cms.untracked.vdouble(0.000785,0.000864,0.00173),
    xerr_barrel_ln_def_ = cms.untracked.double(0.0102),
    yerr_barrel_ln_ = cms.untracked.vdouble(0.00272,0.00142,0.00216,0.00238,0.0026,0.00298,0.00298,0.00253,0.118),
    yerr_barrel_ln_def_ = cms.untracked.double(0.125),
    xerr_endcap_ = cms.untracked.vdouble(0.00144,0.00069,0.00441),
    xerr_endcap_def_ = cms.untracked.double(0.00757),
    yerr_endcap_ = cms.untracked.vdouble(0.00266,0.00114,0.00766),
    yerr_endcap_def_ = cms.untracked.double(0.0407),
    xerr_barrel_l1_ = cms.untracked.vdouble(0.00167,0.00108,0.00138),
    xerr_barrel_l1_def_ = cms.untracked.double(0.0187),
    yerr_barrel_l1_ = cms.untracked.vdouble(0.00314,0.00369,0.00534,0.0062,0.00732,0.00743,0.00765,0.00809,0.00865),
    yerr_barrel_l1_def_ = cms.untracked.double(0.00829)
    )

# phase1 100x150 um2  285/285 um 2000e Ageing New tuning for 300/fb (2014.07.15)
pixelCPE_phase1_300fb_new = cms.PSet(
    xerr_barrel_ln_ = cms.untracked.vdouble(0.000677,0.000791,0.00177),
    xerr_barrel_ln_def_ = cms.untracked.double(0.0103),
    yerr_barrel_ln_ = cms.untracked.vdouble(0.00264,0.00119,0.00184,0.0019,0.002,0.00221,0.00224,0.00195,0.119),
    yerr_barrel_ln_def_ = cms.untracked.double(0.127),
    xerr_endcap_ = cms.untracked.vdouble(0.00129,0.000591,0.0042),
    xerr_endcap_def_ = cms.untracked.double(0.0108),
    yerr_endcap_ = cms.untracked.vdouble(0.00247,0.000765,0.00677),
    yerr_endcap_def_ = cms.untracked.double(0.0482),
    xerr_barrel_l1_ = cms.untracked.vdouble(0.00126,0.000741,0.00139),
    xerr_barrel_l1_def_ = cms.untracked.double(0.014),
    yerr_barrel_l1_ = cms.untracked.vdouble(0.0025,0.00143,0.00214,0.00231,0.00244,0.00248,0.00251,0.0025,0.00253),
    yerr_barrel_l1_def_ = cms.untracked.double(0.00257)
    )

# phase1 100x150 um2  285/285 um 2000e Ageing New tuning 500/fb  (2014.07.15)
pixelCPE_phase1_500fb_new = cms.PSet(
    xerr_barrel_ln_ = cms.untracked.vdouble(0.000698,0.000802,0.00176),
    xerr_barrel_ln_def_ = cms.untracked.double(0.0103),
    yerr_barrel_ln_ = cms.untracked.vdouble(0.00267,0.00127,0.00194,0.00209,0.00219,0.00233,0.00236,0.00205,0.121),
    yerr_barrel_ln_def_ = cms.untracked.double(0.041),
    xerr_endcap_ = cms.untracked.vdouble(0.00129,0.000591,0.00424),
    xerr_endcap_def_ = cms.untracked.double(0.0103),
    yerr_endcap_ = cms.untracked.vdouble(0.00247,0.000765,0.00675),
    yerr_endcap_def_ = cms.untracked.double(0.026),
    xerr_barrel_l1_ = cms.untracked.vdouble(0.0013,0.000807,0.00137),
    xerr_barrel_l1_def_ = cms.untracked.double(0.0145),
    yerr_barrel_l1_ = cms.untracked.vdouble(0.00238,0.00187,0.00276,0.00308,0.00333,0.00339,0.00345,0.00346,0.00351),
    yerr_barrel_l1_def_ = cms.untracked.double(0.00359)
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
                  'pixelCPE_phase1_0fb'    : pixelCPE_phase1_0fb   ,
                  'pixelCPE_phase1_300fb'    : pixelCPE_phase1_300fb   ,
                  'pixelCPE_phase1_500fb'    : pixelCPE_phase1_500fb   ,
                  'pixelCPE_phase1_300fb_new': pixelCPE_phase1_300fb_new   ,
                  'pixelCPE_phase1_500fb_new': pixelCPE_phase1_500fb_new   ,   
                  'pixelCPE_dummy'           : pixelCPE_dummy
                  }
