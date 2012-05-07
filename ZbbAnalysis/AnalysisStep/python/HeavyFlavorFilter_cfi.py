##*****USAGE*****
## status2 (3): if True checks HF quarks with status 2 (3)
## hDaughterVeto: if True checks only HF quarks with strings or clusters (i.e. hadrons) as daughters
## zDaughterVeto: if True checks only HF quarks not having Z as daughter (i.e. not initial state partons)
##

import FWCore.ParameterSet.Config as cms

heavyflavorfilter = cms.EDFilter('HeavyFlavorFilter',
                                 src= cms.InputTag("genParticles"),
                                 status2 = cms.bool(True),
                                 status3 = cms.bool(True),
                                 bOn = cms.bool(True),
                                 cOn = cms.bool(True),
                                 hDaughterVeto = cms.bool(False),
                                 zDaughterVeto = cms.bool(False),
                                 ptcut=cms.double(15)
                                 )


