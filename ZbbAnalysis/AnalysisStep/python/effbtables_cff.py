import FWCore.ParameterSet.Config as cms

## MC EFF b SSVHEM
mc_effb_ssvhem_2011 = cms.VPSet()
mc_effb_ssvhem_2011.extend([
    cms.PSet(effb=cms.untracked.double(0.138),ptrange=cms.untracked.vdouble(10,20)),
    cms.PSet(effb=cms.untracked.double(0.330),ptrange=cms.untracked.vdouble(20,30)),
    cms.PSet(effb=cms.untracked.double(0.473),ptrange=cms.untracked.vdouble(30,40)),
    cms.PSet(effb=cms.untracked.double(0.572),ptrange=cms.untracked.vdouble(40,50)),
    cms.PSet(effb=cms.untracked.double(0.626),ptrange=cms.untracked.vdouble(50,60)),
    cms.PSet(effb=cms.untracked.double(0.657),ptrange=cms.untracked.vdouble(60,70)),
    cms.PSet(effb=cms.untracked.double(0.681),ptrange=cms.untracked.vdouble(70,80)),
    cms.PSet(effb=cms.untracked.double(0.683),ptrange=cms.untracked.vdouble(80,100)),
    cms.PSet(effb=cms.untracked.double(0.683),ptrange=cms.untracked.vdouble(100,120)),
    cms.PSet(effb=cms.untracked.double(0.652),ptrange=cms.untracked.vdouble(120,999999))
    ])

## MC EFF b SSVHPT
mc_effb_ssvhpt_2011 = cms.VPSet()
mc_effb_ssvhpt_2011.extend([
    cms.PSet(effb=cms.untracked.double(0.018),ptrange=cms.untracked.vdouble(10,20)),
    cms.PSet(effb=cms.untracked.double(0.120),ptrange=cms.untracked.vdouble(20,30)),
    cms.PSet(effb=cms.untracked.double(0.236),ptrange=cms.untracked.vdouble(30,40)),
    cms.PSet(effb=cms.untracked.double(0.333),ptrange=cms.untracked.vdouble(40,50)),
    cms.PSet(effb=cms.untracked.double(0.389),ptrange=cms.untracked.vdouble(50,60)),
    cms.PSet(effb=cms.untracked.double(0.432),ptrange=cms.untracked.vdouble(60,70)),
    cms.PSet(effb=cms.untracked.double(0.463),ptrange=cms.untracked.vdouble(70,80)),
    cms.PSet(effb=cms.untracked.double(0.482),ptrange=cms.untracked.vdouble(80,100)),
    cms.PSet(effb=cms.untracked.double(0.478),ptrange=cms.untracked.vdouble(100,120)),
    cms.PSet(effb=cms.untracked.double(0.457),ptrange=cms.untracked.vdouble(120,999999))
    ])
