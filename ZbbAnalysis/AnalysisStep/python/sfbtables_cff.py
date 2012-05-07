import FWCore.ParameterSet.Config as cms

######################################################################################
# SFb for Fall10 (Winter10?) and 2010 data                                           #
# based on http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2011_114_v2.pdf #
######################################################################################

##SSVHEM
ssvhem_2010 = cms.VPSet()
ssvhem_2010.extend([
    cms.PSet(scale_factor=cms.untracked.double(0.909),ptrange=cms.untracked.vdouble(20,40)),
    cms.PSet(scale_factor=cms.untracked.double(0.909),ptrange=cms.untracked.vdouble(40,50)),
    cms.PSet(scale_factor=cms.untracked.double(0.841),ptrange=cms.untracked.vdouble(50,60)),
    cms.PSet(scale_factor=cms.untracked.double(0.857),ptrange=cms.untracked.vdouble(60,80)),
    cms.PSet(scale_factor=cms.untracked.double(0.753),ptrange=cms.untracked.vdouble(80,140)),
    cms.PSet(scale_factor=cms.untracked.double(0.617),ptrange=cms.untracked.vdouble(140,999999))
    ])

##SSVHPT
ssvhpt_2010 = cms.VPSet()
ssvhpt_2010.extend([
    cms.PSet(scale_factor=cms.untracked.double(0.956),ptrange=cms.untracked.vdouble(20,40)),
    cms.PSet(scale_factor=cms.untracked.double(0.956),ptrange=cms.untracked.vdouble(40,50)),
    cms.PSet(scale_factor=cms.untracked.double(0.927),ptrange=cms.untracked.vdouble(50,60)),
    cms.PSet(scale_factor=cms.untracked.double(0.854),ptrange=cms.untracked.vdouble(60,80)),
    cms.PSet(scale_factor=cms.untracked.double(0.824),ptrange=cms.untracked.vdouble(80,140)),
    cms.PSet(scale_factor=cms.untracked.double(0.713),ptrange=cms.untracked.vdouble(140,999999))
    ])

#########################################################################################################
# SFb for Spring11 (+DA patch) or Summer11 MC and 2011 data                                             #
# based on https://indico.cern.ch/getFile.py/access?contribId=2&resId=1&materialId=slides&confId=130982 #
#########################################################################################################

##SSVHEM
ssvhem_2011 = cms.VPSet()
ssvhem_2011.extend([
    cms.PSet(scale_factor=cms.untracked.double(0.873),ptrange=cms.untracked.vdouble(15,20)),
    cms.PSet(scale_factor=cms.untracked.double(1.069),ptrange=cms.untracked.vdouble(20,30)),
    cms.PSet(scale_factor=cms.untracked.double(1.014),ptrange=cms.untracked.vdouble(30,40)),
    cms.PSet(scale_factor=cms.untracked.double(0.981),ptrange=cms.untracked.vdouble(40,50)),
    cms.PSet(scale_factor=cms.untracked.double(0.977),ptrange=cms.untracked.vdouble(50,60)),
    cms.PSet(scale_factor=cms.untracked.double(0.892),ptrange=cms.untracked.vdouble(60,70)),
    cms.PSet(scale_factor=cms.untracked.double(0.982),ptrange=cms.untracked.vdouble(70,80)),
    cms.PSet(scale_factor=cms.untracked.double(0.966),ptrange=cms.untracked.vdouble(80,100)),
    cms.PSet(scale_factor=cms.untracked.double(0.947),ptrange=cms.untracked.vdouble(100,120)),
    cms.PSet(scale_factor=cms.untracked.double(0.795),ptrange=cms.untracked.vdouble(120,999999))
    ])

##SSVHPT
ssvhpt_2011 = cms.VPSet()
ssvhpt_2011.extend([
    cms.PSet(scale_factor=cms.untracked.double(0.802),ptrange=cms.untracked.vdouble(15,20)),
    cms.PSet(scale_factor=cms.untracked.double(1.000),ptrange=cms.untracked.vdouble(20,30)),
    cms.PSet(scale_factor=cms.untracked.double(0.995),ptrange=cms.untracked.vdouble(30,40)),
    cms.PSet(scale_factor=cms.untracked.double(0.948),ptrange=cms.untracked.vdouble(40,50)),
    cms.PSet(scale_factor=cms.untracked.double(0.912),ptrange=cms.untracked.vdouble(50,60)),
    cms.PSet(scale_factor=cms.untracked.double(0.829),ptrange=cms.untracked.vdouble(60,70)),
    cms.PSet(scale_factor=cms.untracked.double(0.948),ptrange=cms.untracked.vdouble(70,80)),
    cms.PSet(scale_factor=cms.untracked.double(0.900),ptrange=cms.untracked.vdouble(80,100)),
    cms.PSet(scale_factor=cms.untracked.double(0.863),ptrange=cms.untracked.vdouble(100,120)),
    cms.PSet(scale_factor=cms.untracked.double(0.794),ptrange=cms.untracked.vdouble(120,999999))
    ])


