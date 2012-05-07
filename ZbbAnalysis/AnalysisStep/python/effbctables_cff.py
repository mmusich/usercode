import FWCore.ParameterSet.Config as cms

## MC EFF b/c SSVHEM
mc_effbc_ssvhem_2011 = cms.VPSet()
mc_effbc_ssvhem_2011.extend([
cms.PSet(ptrange=cms.untracked.vdouble(20,30),effb_barrel=cms.untracked.double(0.439),effb_forward=cms.untracked.double(0.346),effc_barrel=cms.untracked.double(0.109),effc_forward=cms.untracked.double(0.109)),
cms.PSet(ptrange=cms.untracked.vdouble(30,40),effb_barrel=cms.untracked.double(0.527),effb_forward=cms.untracked.double(0.455),effc_barrel=cms.untracked.double(0.143),effc_forward=cms.untracked.double(0.143)),
cms.PSet(ptrange=cms.untracked.vdouble(40,50),effb_barrel=cms.untracked.double(0.617),effb_forward=cms.untracked.double(0.564),effc_barrel=cms.untracked.double(0.144),effc_forward=cms.untracked.double(0.144)),
cms.PSet(ptrange=cms.untracked.vdouble(50,60),effb_barrel=cms.untracked.double(0.696),effb_forward=cms.untracked.double(0.629),effc_barrel=cms.untracked.double(0.206),effc_forward=cms.untracked.double(0.206)),
cms.PSet(ptrange=cms.untracked.vdouble(60,70),effb_barrel=cms.untracked.double(0.707),effb_forward=cms.untracked.double(0.671),effc_barrel=cms.untracked.double(0.215),effc_forward=cms.untracked.double(0.215)),
cms.PSet(ptrange=cms.untracked.vdouble(70,80),effb_barrel=cms.untracked.double(0.71),effb_forward=cms.untracked.double(0.71),effc_barrel=cms.untracked.double(0.224),effc_forward=cms.untracked.double(0.224)),
cms.PSet(ptrange=cms.untracked.vdouble(80,100),effb_barrel=cms.untracked.double(0.725),effb_forward=cms.untracked.double(0.712),effc_barrel=cms.untracked.double(0.211),effc_forward=cms.untracked.double(0.211)),
cms.PSet(ptrange=cms.untracked.vdouble(100,120),effb_barrel=cms.untracked.double(0.735),effb_forward=cms.untracked.double(0.659),effc_barrel=cms.untracked.double(0.249),effc_forward=cms.untracked.double(0.249)),
cms.PSet(ptrange=cms.untracked.vdouble(120,240),effb_barrel=cms.untracked.double(0.735),effb_forward=cms.untracked.double(0.695),effc_barrel=cms.untracked.double(0.242),effc_forward=cms.untracked.double(0.242))
])


## MC EFF b/c SSVHPT
mc_effbc_ssvhpt_2011 = cms.VPSet()
mc_effbc_ssvhpt_2011.extend([
cms.PSet(ptrange=cms.untracked.vdouble(20,30),effb_barrel=cms.untracked.double(0.262),effb_forward=cms.untracked.double(0.193),effc_barrel=cms.untracked.double(0.033),effc_forward=cms.untracked.double(0.033)),
cms.PSet(ptrange=cms.untracked.vdouble(30,40),effb_barrel=cms.untracked.double(0.347),effb_forward=cms.untracked.double(0.272),effc_barrel=cms.untracked.double(0.0615),effc_forward=cms.untracked.double(0.0615)),
cms.PSet(ptrange=cms.untracked.vdouble(40,50),effb_barrel=cms.untracked.double(0.454),effb_forward=cms.untracked.double(0.386),effc_barrel=cms.untracked.double(0.0451),effc_forward=cms.untracked.double(0.0451)),
cms.PSet(ptrange=cms.untracked.vdouble(50,60),effb_barrel=cms.untracked.double(0.527),effb_forward=cms.untracked.double(0.438),effc_barrel=cms.untracked.double(0.096),effc_forward=cms.untracked.double(0.096)),
cms.PSet(ptrange=cms.untracked.vdouble(60,70),effb_barrel=cms.untracked.double(0.559),effb_forward=cms.untracked.double(0.486),effc_barrel=cms.untracked.double(0.0914),effc_forward=cms.untracked.double(0.0914)),
cms.PSet(ptrange=cms.untracked.vdouble(70,80),effb_barrel=cms.untracked.double(0.549),effb_forward=cms.untracked.double(0.538),effc_barrel=cms.untracked.double(0.114),effc_forward=cms.untracked.double(0.114)),
cms.PSet(ptrange=cms.untracked.vdouble(80,100),effb_barrel=cms.untracked.double(0.571),effb_forward=cms.untracked.double(0.526),effc_barrel=cms.untracked.double(0.088),effc_forward=cms.untracked.double(0.088)),
cms.PSet(ptrange=cms.untracked.vdouble(100,120),effb_barrel=cms.untracked.double(0.571),effb_forward=cms.untracked.double(0.537),effc_barrel=cms.untracked.double(0.114),effc_forward=cms.untracked.double(0.114)),
cms.PSet(ptrange=cms.untracked.vdouble(120,240),effb_barrel=cms.untracked.double(0.589),effb_forward=cms.untracked.double(0.539),effc_barrel=cms.untracked.double(0.0838),effc_forward=cms.untracked.double(0.0838))
])



## MC EFF b/c JPT
mc_effbc_jpt_2011 = cms.VPSet()
mc_effbc_jpt_2011.extend([
cms.PSet(ptrange=cms.untracked.vdouble(20,30),effb_barrel=cms.untracked.double(0.3),effb_forward=cms.untracked.double(0.199),effc_barrel=cms.untracked.double(0.0365),effc_forward=cms.untracked.double(0.0365)),
cms.PSet(ptrange=cms.untracked.vdouble(30,40),effb_barrel=cms.untracked.double(0.378),effb_forward=cms.untracked.double(0.285),effc_barrel=cms.untracked.double(0.0615),effc_forward=cms.untracked.double(0.0615)),
cms.PSet(ptrange=cms.untracked.vdouble(40,50),effb_barrel=cms.untracked.double(0.48),effb_forward=cms.untracked.double(0.402),effc_barrel=cms.untracked.double(0.0388),effc_forward=cms.untracked.double(0.0388)),
cms.PSet(ptrange=cms.untracked.vdouble(50,60),effb_barrel=cms.untracked.double(0.53),effb_forward=cms.untracked.double(0.434),effc_barrel=cms.untracked.double(0.0837),effc_forward=cms.untracked.double(0.0837)),
cms.PSet(ptrange=cms.untracked.vdouble(60,70),effb_barrel=cms.untracked.double(0.546),effb_forward=cms.untracked.double(0.447),effc_barrel=cms.untracked.double(0.0538),effc_forward=cms.untracked.double(0.0538)),
cms.PSet(ptrange=cms.untracked.vdouble(70,80),effb_barrel=cms.untracked.double(0.564),effb_forward=cms.untracked.double(0.538),effc_barrel=cms.untracked.double(0.0534),effc_forward=cms.untracked.double(0.0534)),
cms.PSet(ptrange=cms.untracked.vdouble(80,100),effb_barrel=cms.untracked.double(0.579),effb_forward=cms.untracked.double(0.505),effc_barrel=cms.untracked.double(0.0641),effc_forward=cms.untracked.double(0.0641)),
cms.PSet(ptrange=cms.untracked.vdouble(100,120),effb_barrel=cms.untracked.double(0.572),effb_forward=cms.untracked.double(0.526),effc_barrel=cms.untracked.double(0.0976),effc_forward=cms.untracked.double(0.0976)),
cms.PSet(ptrange=cms.untracked.vdouble(120,240),effb_barrel=cms.untracked.double(0.594),effb_forward=cms.untracked.double(0.512),effc_barrel=cms.untracked.double(0.0872),effc_forward=cms.untracked.double(0.0872))
])



## MC EFF b/c CSVM
mc_effbc_csvm_2011 = cms.VPSet()
mc_effbc_csvm_2011.extend([
cms.PSet(ptrange=cms.untracked.vdouble(20,30),effb_barrel=cms.untracked.double(0.565),effb_forward=cms.untracked.double(0.5),effc_barrel=cms.untracked.double(0.152),effc_forward=cms.untracked.double(0.152)),
cms.PSet(ptrange=cms.untracked.vdouble(30,40),effb_barrel=cms.untracked.double(0.631),effb_forward=cms.untracked.double(0.582),effc_barrel=cms.untracked.double(0.203),effc_forward=cms.untracked.double(0.203)),
cms.PSet(ptrange=cms.untracked.vdouble(40,50),effb_barrel=cms.untracked.double(0.687),effb_forward=cms.untracked.double(0.642),effc_barrel=cms.untracked.double(0.153),effc_forward=cms.untracked.double(0.153)),
cms.PSet(ptrange=cms.untracked.vdouble(50,60),effb_barrel=cms.untracked.double(0.749),effb_forward=cms.untracked.double(0.698),effc_barrel=cms.untracked.double(0.256),effc_forward=cms.untracked.double(0.256)),
cms.PSet(ptrange=cms.untracked.vdouble(60,70),effb_barrel=cms.untracked.double(0.742),effb_forward=cms.untracked.double(0.717),effc_barrel=cms.untracked.double(0.251),effc_forward=cms.untracked.double(0.251)),
cms.PSet(ptrange=cms.untracked.vdouble(70,80),effb_barrel=cms.untracked.double(0.746),effb_forward=cms.untracked.double(0.767),effc_barrel=cms.untracked.double(0.221),effc_forward=cms.untracked.double(0.221)),
cms.PSet(ptrange=cms.untracked.vdouble(80,100),effb_barrel=cms.untracked.double(0.768),effb_forward=cms.untracked.double(0.762),effc_barrel=cms.untracked.double(0.208),effc_forward=cms.untracked.double(0.208)),
cms.PSet(ptrange=cms.untracked.vdouble(100,120),effb_barrel=cms.untracked.double(0.794),effb_forward=cms.untracked.double(0.715),effc_barrel=cms.untracked.double(0.25),effc_forward=cms.untracked.double(0.25)),
cms.PSet(ptrange=cms.untracked.vdouble(120,240),effb_barrel=cms.untracked.double(0.809),effb_forward=cms.untracked.double(0.769),effc_barrel=cms.untracked.double(0.229),effc_forward=cms.untracked.double(0.229))
])


## MC EFF b/c CSVT
mc_effbc_csvt_2011 = cms.VPSet()
mc_effbc_csvt_2011.extend([
cms.PSet(ptrange=cms.untracked.vdouble(20,30),effb_barrel=cms.untracked.double(0.411),effb_forward=cms.untracked.double(0.333),effc_barrel=cms.untracked.double(0.0515),effc_forward=cms.untracked.double(0.0515)),
cms.PSet(ptrange=cms.untracked.vdouble(30,40),effb_barrel=cms.untracked.double(0.489),effb_forward=cms.untracked.double(0.425),effc_barrel=cms.untracked.double(0.0864),effc_forward=cms.untracked.double(0.0864)),
cms.PSet(ptrange=cms.untracked.vdouble(40,50),effb_barrel=cms.untracked.double(0.553),effb_forward=cms.untracked.double(0.513),effc_barrel=cms.untracked.double(0.0528),effc_forward=cms.untracked.double(0.0528)),
cms.PSet(ptrange=cms.untracked.vdouble(50,60),effb_barrel=cms.untracked.double(0.622),effb_forward=cms.untracked.double(0.553),effc_barrel=cms.untracked.double(0.0862),effc_forward=cms.untracked.double(0.0862)),
cms.PSet(ptrange=cms.untracked.vdouble(60,70),effb_barrel=cms.untracked.double(0.618),effb_forward=cms.untracked.double(0.567),effc_barrel=cms.untracked.double(0.0605),effc_forward=cms.untracked.double(0.0605)),
cms.PSet(ptrange=cms.untracked.vdouble(70,80),effb_barrel=cms.untracked.double(0.621),effb_forward=cms.untracked.double(0.62),effc_barrel=cms.untracked.double(0.072),effc_forward=cms.untracked.double(0.072)),
cms.PSet(ptrange=cms.untracked.vdouble(80,100),effb_barrel=cms.untracked.double(0.621),effb_forward=cms.untracked.double(0.6),effc_barrel=cms.untracked.double(0.0692),effc_forward=cms.untracked.double(0.0692)),
cms.PSet(ptrange=cms.untracked.vdouble(100,120),effb_barrel=cms.untracked.double(0.614),effb_forward=cms.untracked.double(0.598),effc_barrel=cms.untracked.double(0.0446),effc_forward=cms.untracked.double(0.0446)),
cms.PSet(ptrange=cms.untracked.vdouble(120,240),effb_barrel=cms.untracked.double(0.621),effb_forward=cms.untracked.double(0.543),effc_barrel=cms.untracked.double(0.0579),effc_forward=cms.untracked.double(0.0579))
])
