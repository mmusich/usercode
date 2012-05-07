import FWCore.ParameterSet.Config as cms
## MC EFF b/c SSVHEM
mc_effbc_ssvhem_2011 = cms.VPSet()
mc_effbc_ssvhem_2011.extend([
cms.PSet(ptrange=cms.untracked.vdouble(20,30),effb_barrel=cms.untracked.double(0.403),effb_forward=cms.untracked.double(0.37),effc_barrel=cms.untracked.double(0.112),effc_forward=cms.untracked.double(0.112)),
cms.PSet(ptrange=cms.untracked.vdouble(30,40),effb_barrel=cms.untracked.double(0.499),effb_forward=cms.untracked.double(0.467),effc_barrel=cms.untracked.double(0.148),effc_forward=cms.untracked.double(0.148)),
cms.PSet(ptrange=cms.untracked.vdouble(40,50),effb_barrel=cms.untracked.double(0.604),effb_forward=cms.untracked.double(0.541),effc_barrel=cms.untracked.double(0.137),effc_forward=cms.untracked.double(0.137)),
cms.PSet(ptrange=cms.untracked.vdouble(50,60),effb_barrel=cms.untracked.double(0.706),effb_forward=cms.untracked.double(0.608),effc_barrel=cms.untracked.double(0.223),effc_forward=cms.untracked.double(0.223)),
cms.PSet(ptrange=cms.untracked.vdouble(60,70),effb_barrel=cms.untracked.double(0.767),effb_forward=cms.untracked.double(0.633),effc_barrel=cms.untracked.double(0.208),effc_forward=cms.untracked.double(0.208)),
cms.PSet(ptrange=cms.untracked.vdouble(70,80),effb_barrel=cms.untracked.double(0.689),effb_forward=cms.untracked.double(0.751),effc_barrel=cms.untracked.double(0.235),effc_forward=cms.untracked.double(0.235)),
cms.PSet(ptrange=cms.untracked.vdouble(80,100),effb_barrel=cms.untracked.double(0.673),effb_forward=cms.untracked.double(0.743),effc_barrel=cms.untracked.double(0.202),effc_forward=cms.untracked.double(0.202)),
cms.PSet(ptrange=cms.untracked.vdouble(100,120),effb_barrel=cms.untracked.double(0.731),effb_forward=cms.untracked.double(0.669),effc_barrel=cms.untracked.double(0.239),effc_forward=cms.untracked.double(0.239)),
cms.PSet(ptrange=cms.untracked.vdouble(120,240),effb_barrel=cms.untracked.double(0.749),effb_forward=cms.untracked.double(0.682),effc_barrel=cms.untracked.double(0.255),effc_forward=cms.untracked.double(0.255))
])


## MC EFF b/c SSVHPT
mc_effbc_ssvhpt_2011 = cms.VPSet()
mc_effbc_ssvhpt_2011.extend([
cms.PSet(ptrange=cms.untracked.vdouble(20,30),effb_barrel=cms.untracked.double(0.221),effb_forward=cms.untracked.double(0.205),effc_barrel=cms.untracked.double(0.0277),effc_forward=cms.untracked.double(0.0277)),
cms.PSet(ptrange=cms.untracked.vdouble(30,40),effb_barrel=cms.untracked.double(0.315),effb_forward=cms.untracked.double(0.261),effc_barrel=cms.untracked.double(0.0623),effc_forward=cms.untracked.double(0.0623)),
cms.PSet(ptrange=cms.untracked.vdouble(40,50),effb_barrel=cms.untracked.double(0.472),effb_forward=cms.untracked.double(0.382),effc_barrel=cms.untracked.double(0.0427),effc_forward=cms.untracked.double(0.0427)),
cms.PSet(ptrange=cms.untracked.vdouble(50,60),effb_barrel=cms.untracked.double(0.536),effb_forward=cms.untracked.double(0.432),effc_barrel=cms.untracked.double(0.114),effc_forward=cms.untracked.double(0.114)),
cms.PSet(ptrange=cms.untracked.vdouble(60,70),effb_barrel=cms.untracked.double(0.575),effb_forward=cms.untracked.double(0.489),effc_barrel=cms.untracked.double(0.0818),effc_forward=cms.untracked.double(0.0818)),
cms.PSet(ptrange=cms.untracked.vdouble(70,80),effb_barrel=cms.untracked.double(0.48),effb_forward=cms.untracked.double(0.58),effc_barrel=cms.untracked.double(0.103),effc_forward=cms.untracked.double(0.103)),
cms.PSet(ptrange=cms.untracked.vdouble(80,100),effb_barrel=cms.untracked.double(0.549),effb_forward=cms.untracked.double(0.568),effc_barrel=cms.untracked.double(0.0923),effc_forward=cms.untracked.double(0.0923)),
cms.PSet(ptrange=cms.untracked.vdouble(100,120),effb_barrel=cms.untracked.double(0.592),effb_forward=cms.untracked.double(0.506),effc_barrel=cms.untracked.double(0.092),effc_forward=cms.untracked.double(0.092)),
cms.PSet(ptrange=cms.untracked.vdouble(120,240),effb_barrel=cms.untracked.double(0.615),effb_forward=cms.untracked.double(0.528),effc_barrel=cms.untracked.double(0.0907),effc_forward=cms.untracked.double(0.0907))
])
