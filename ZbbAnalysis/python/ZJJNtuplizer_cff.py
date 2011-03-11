import FWCore.ParameterSet.Config as cms
import copy

#   ______                  _       _     _          
#  |___  /                 (_)     | |   | |         
#     / /  __   ____ _ _ __ _  __ _| |__ | | ___ ___ 
#    / /   \ \ / / _` | '__| |/ _` | '_ \| |/ _ | __|
#   / /__   \ V / (_| | |  | | (_| | |_) | |  __|__ \
#  /_____|   \_/ \__,_|_|  |_|\__,_|_.__/|_|\___|___/


zed =  cms.EDProducer(
    "CandViewNtpProducer",
    src=cms.InputTag("zMMCand"),
    lazyParser=cms.untracked.bool(True),
    prefix=cms.untracked.string("Zmm"),
    eventInfo=cms.untracked.bool(True),
    variables = cms.VPSet(
    ## Candidate Mass
    cms.PSet(
    tag = cms.untracked.string("Mass"),
    quantity = cms.untracked.string("mass")
    ),
    ## Candidate pt
    cms.PSet(
    tag = cms.untracked.string("Pt"),
    quantity = cms.untracked.string("pt")
    ),
    ## Candidate eta
    cms.PSet(
    tag = cms.untracked.string("Eta"),
    quantity = cms.untracked.string("eta")
    ),
    ## Candidate y
    cms.PSet(
    tag = cms.untracked.string("Y"),
    quantity = cms.untracked.string("y")
    ),
    ## Candidate phi
    cms.PSet(
    tag = cms.untracked.string("Phi"),
    quantity = cms.untracked.string("phi")
    ),
    ## Candidate fourMomentum
    cms.PSet(
    tag = cms.untracked.string("energy"),
    quantity = cms.untracked.string("energy")
    ),
    cms.PSet(
    tag = cms.untracked.string("px"),
    quantity = cms.untracked.string("px")
    ),
    cms.PSet(
    tag = cms.untracked.string("py"),
    quantity = cms.untracked.string("py")
    ),
    cms.PSet(
    tag = cms.untracked.string("pz"),
    quantity = cms.untracked.string("pz")
    ),
    )
    )

#   ______  _       _                                             
#  |___  / (_)     | |                                            
#     / /   _ _ __ | |_  ___    _ __ ___  _   _   _ __ ___  _   _ 
#    / /   | | '_ \| __|/ _ \  | '_ ` _ \| | | | | '_ ` _ \| | | |
#   / /__  | | | | | |_| (_) | | | | | | | |_| | | | | | | | |_| |
#  /_____| |_|_| |_|\__|\___/  |_| |_| |_|\__,_| |_| |_| |_|\__,_|
                                                                

zmm = (

    ## leptons pt
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Pt"),
    quantity = cms.untracked.string("daughter(0).pt ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Pt"),
    quantity = cms.untracked.string("daughter(1).pt")
    ),
    ## leptons eta
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Eta"),
    quantity = cms.untracked.string("daughter(0).eta ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Eta"),
    quantity = cms.untracked.string("daughter(1).eta ")
    ),
    ## leptons phi
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Phi"),
    quantity = cms.untracked.string("daughter(0).phi ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Phi"),
    quantity = cms.untracked.string("daughter(1).phi ")
    ),
    ## leptons charge
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Q"),
    quantity = cms.untracked.string("daughter(0).charge")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Q"),
    quantity = cms.untracked.string("daughter(1).charge")
    ),

    ## leptons fourMomentum
    cms.PSet(
    tag = cms.untracked.string("LeptDau1energy"),
    quantity = cms.untracked.string("daughter(0).energy")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1px"),
    quantity = cms.untracked.string("daughter(0).px")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1py"),
    quantity = cms.untracked.string("daughter(0).py")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1pz"),
    quantity = cms.untracked.string("daughter(0).pz")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2energy"),
    quantity = cms.untracked.string("daughter(1).energy")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2px"),
    quantity = cms.untracked.string("daughter(1).px")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2py"),
    quantity = cms.untracked.string("daughter(1).py")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2pz"),
    quantity = cms.untracked.string("daughter(1).pz")
    ),
    
    ## Trigger Matching
    
    cms.PSet(
    tag = cms.untracked.string("TriggerMatchDau1"),
    quantity = cms.untracked.string("? !daughter(0).masterClone.triggerObjectMatchesByPath('HLT_Mu9').empty || !daughter(0).masterClone.triggerObjectMatchesByPath('HLT_Mu11').empty || !daughter(0).masterClone.triggerObjectMatchesByPath('HLT_Mu15_v1').empty?1:0")
    ),
    
    cms.PSet(
    tag = cms.untracked.string("TriggerMatchDau2"),
    quantity = cms.untracked.string("? !daughter(1).masterClone.triggerObjectMatchesByPath('HLT_Mu9').empty || !daughter(1).masterClone.triggerObjectMatchesByPath('HLT_Mu11').empty || !daughter(1).masterClone.triggerObjectMatchesByPath('HLT_Mu15_v1').empty?1:0")
    ),
        
    cms.PSet(
    tag = cms.untracked.string("LeptDau1GlobalMuonBit"),
    quantity = cms.untracked.string("daughter(0).isGlobalMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2GlobalMuonBit"),
    quantity = cms.untracked.string("daughter(1).isGlobalMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1StandAloneBit"),
    quantity = cms.untracked.string("daughter(0).isStandAloneMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2StandAloneBit"),
    quantity = cms.untracked.string("daughter(1).isStandAloneMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1TrackerMuonBit"),
    quantity = cms.untracked.string("daughter(0).isTrackerMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2TrackerMuonBit"),
    quantity = cms.untracked.string("daughter(1).isTrackerMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofMuonHits"),
    quantity = cms.untracked.string("?daughter(0).isGlobalMuon?daughter(0).masterClone.globalTrack.hitPattern.numberOfValidMuonHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2NofMuonHits"),
    quantity = cms.untracked.string("?daughter(1).isGlobalMuon?daughter(1).masterClone.globalTrack.hitPattern.numberOfValidMuonHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofStripHits"),
    quantity = cms.untracked.string("?daughter(0).isGlobalMuon?daughter(0).masterClone.globalTrack.hitPattern.numberOfValidStripHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2NofStripHits"),
    quantity = cms.untracked.string("?daughter(1).isGlobalMuon?daughter(1).masterClone.globalTrack.hitPattern.numberOfValidStripHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofPixelHits"),
    quantity = cms.untracked.string("?daughter(0).isGlobalMuon?daughter(0).masterClone.globalTrack.hitPattern.numberOfValidPixelHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2NofPixelHits"),
    quantity = cms.untracked.string("?daughter(1).isGlobalMuon?daughter(1).masterClone.globalTrack.hitPattern.numberOfValidPixelHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NormChi2"),
    quantity = cms.untracked.string("?daughter(0).isGlobalMuon?daughter(0).masterClone.globalTrack.normalizedChi2: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2NormChi2"),
    quantity = cms.untracked.string("?daughter(1).isGlobalMuon?daughter(1).masterClone.globalTrack.normalizedChi2: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofChambers"),
    quantity = cms.untracked.string("?daughter(0).isGlobalMuon?daughter(0).masterClone.numberOfChambers: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2NofChambers"),
    quantity = cms.untracked.string("?daughter(1).isGlobalMuon?daughter(1).masterClone.numberOfChambers: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1dB"),
    quantity = cms.untracked.string("daughter(0).masterClone.dB")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2dB"),
    quantity = cms.untracked.string("daughter(1).masterClone.dB")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1TrkIso"),
    quantity = cms.untracked.string("daughter(0).masterClone.trackIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2TrkIso"),
    quantity = cms.untracked.string("daughter(1).masterClone.trackIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1EcalIso"),
    quantity = cms.untracked.string("daughter(0).masterClone.ecalIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2EcalIso"),
    quantity = cms.untracked.string("daughter(1).masterClone.ecalIso")
    ),
      cms.PSet(
    tag = cms.untracked.string("LeptDau1HcalIso"),
    quantity = cms.untracked.string("daughter(0).masterClone.hcalIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2HcalIso"),
    quantity = cms.untracked.string("daughter(1).masterClone.hcalIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1CombRelIso"),
    quantity = cms.untracked.string("(daughter(0).masterClone.hcalIso + daughter(0).masterClone.ecalIso + daughter(0).masterClone.trackIso )/ daughter(0).pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2CombrelIso"),
    quantity = cms.untracked.string("(daughter(1).masterClone.hcalIso + daughter(1).masterClone.ecalIso + daughter(1).masterClone.trackIso )/ daughter(1).pt")
    ),

    )

#   ______  _       _                       
#  |___  / (_)     | |                      
#     / /   _ _ __ | |_  ___     ___    ___ 
#    / /   | | '_ \| __|/ _ \   / _ \  / _ \
#   / /__  | | | | | |_| (_) | |  __/ |  __/
#  /_____| |_|_| |_|\__|\___/   \___|  \___|
                                          
zee =(

    ## leptons pt
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Pt"),
    quantity = cms.untracked.string("daughter(0).pt ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Pt"),
    quantity = cms.untracked.string("daughter(1).pt")
    ),
    ## leptons eta
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Eta"),
    quantity = cms.untracked.string("daughter(0).eta ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Eta"),
    quantity = cms.untracked.string("daughter(1).eta ")
    ),
    ## leptons phi
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Phi"),
    quantity = cms.untracked.string("daughter(0).phi ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Phi"),
    quantity = cms.untracked.string("daughter(1).phi ")
    ),
    ## leptons charge
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Q"),
    quantity = cms.untracked.string("daughter(0).charge")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Q"),
    quantity = cms.untracked.string("daughter(1).charge")
    ),

    ## leptons fourMomentum
    cms.PSet(
    tag = cms.untracked.string("LeptDau1energy"),
    quantity = cms.untracked.string("daughter(0).energy")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1px"),
    quantity = cms.untracked.string("daughter(0).px")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1py"),
    quantity = cms.untracked.string("daughter(0).py")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1pz"),
    quantity = cms.untracked.string("daughter(0).pz")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2energy"),
    quantity = cms.untracked.string("daughter(1).energy")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2px"),
    quantity = cms.untracked.string("daughter(1).px")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2py"),
    quantity = cms.untracked.string("daughter(1).py")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2pz"),
    quantity = cms.untracked.string("daughter(1).pz")
    ),

    ## Electron variables
    
    cms.PSet(
    tag = cms.untracked.string("LeptDau1isEB"),
    quantity = cms.untracked.string("daughter(0).masterClone.isEB")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2isEB"),
    quantity = cms.untracked.string("daughter(1).masterClone.isEB")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1isEE"),
    quantity = cms.untracked.string("daughter(0).masterClone.isEE")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2isEE"),
    quantity = cms.untracked.string("daughter(1).masterClone.isEE")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1ecaliso"),
    quantity = cms.untracked.string("daughter(0).masterClone.dr03EcalRecHitSumEt")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2ecaliso"),
    quantity = cms.untracked.string("daughter(1).masterClone.dr03EcalRecHitSumEt")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1hcaliso"),
    quantity = cms.untracked.string("daughter(0).masterClone.dr03HcalTowerSumEt")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2hcaliso"),
    quantity = cms.untracked.string("daughter(1).masterClone.dr03HcalTowerSumEt")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1trackiso"),
    quantity = cms.untracked.string("daughter(0).masterClone.dr03TkSumPt")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2trackiso"),
    quantity = cms.untracked.string("daughter(1).masterClone.dr03TkSumPt")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1classification"),
    quantity = cms.untracked.string("daughter(0).masterClone.classification")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2classification"),
    quantity = cms.untracked.string("daughter(1).masterClone.classification")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1sihih"),
    quantity = cms.untracked.string("daughter(0).masterClone.sigmaIetaIeta")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2sihih"),
    quantity = cms.untracked.string("daughter(1).masterClone.sigmaIetaIeta")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1dfi"),
    quantity = cms.untracked.string("daughter(0).masterClone.deltaPhiSuperClusterTrackAtVtx")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2dfi"),
    quantity = cms.untracked.string("daughter(1).masterClone.deltaPhiSuperClusterTrackAtVtx")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1dhi"),
    quantity = cms.untracked.string("daughter(0).masterClone.deltaEtaSuperClusterTrackAtVtx")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2dhi"),
    quantity = cms.untracked.string("daughter(1).masterClone.deltaEtaSuperClusterTrackAtVtx")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1hoe"),
    quantity = cms.untracked.string("daughter(0).masterClone.hcalOverEcal")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2hoe"),
    quantity = cms.untracked.string("daughter(1).masterClone.hcalOverEcal")
    ),

    #    FIXME!!!! To be added Electron ID and Trigger Match for Electrons
    
    #     cms.PSet(
    #     tag = cms.untracked.string("EleDau1VBTF80CombID"),
    #     quantity = cms.untracked.string("daughter(0).masterClone.electronID(\"eidVBTFCom80\")")
    #     ),
    #     cms.PSet(
    #     tag = cms.untracked.string("EleDau2VBTF80CombID"),
    #     quantity = cms.untracked.string("daughter(1).masterClone.electronID(\"eidVBTFCom80\")")
    #     )
    #     cms.PSet(
    #     tag = cms.untracked.string("eidRobustTight"),
    #     quantity = cms.untracked.string("?isElectronIDAvailable(\"eidRobustTight\")?electronID(\"eidRobustTight\"):-1")
    #     ),
    )

#        _      _       
#       | |    | |      
#       | | ___| |_ ___ 
#   _   | |/ _ \ __/ __|
#  | |__| |  __/ |_\__ \
#   \____/ \___|\__|___/
                      

jets = (
    
#   _           _                     _             
#  | |         | |                   (_)            
#  | |__ ______| |_  __ _  __ _  __ _ _ _ __   __ _ 
#  | '_ \______| __|/ _` |/ _` |/ _` | | '_ \ / _` |
#  | |_) |     | |_| (_| | (_| | (_| | | | | | (_| |
#  |_.__/       \__|\__,_|\__, |\__, |_|_| |_|\__, |
#                          __/ | __/ |         __/ |
#                         |___/ |___/         |___/ 

    cms.PSet(
    tag = cms.untracked.string("JetCSV"),
    quantity = cms.untracked.string("bDiscriminator(\"combinedSecondaryVertexBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetCSVMVA"),
    quantity = cms.untracked.string("bDiscriminator(\"combinedSecondaryVertexMVABJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetJProb"),
    quantity = cms.untracked.string("bDiscriminator(\"jetProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetJbProb"),
    quantity = cms.untracked.string(" bDiscriminator(\"jetBProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetSSVHE"),
    quantity = cms.untracked.string("bDiscriminator(\"simpleSecondaryVertexHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetSSVHP"),
    quantity = cms.untracked.string(" bDiscriminator(\"simpleSecondaryVertexHighPurBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetElPt"),
    quantity = cms.untracked.string("bDiscriminator(\"softElectronByPtBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetElIp"),
    quantity = cms.untracked.string("bDiscriminator(\"softElectronByIP3dBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetMu"),
    quantity = cms.untracked.string("bDiscriminator(\"softMuonBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetMuPt"),
    quantity = cms.untracked.string("bDiscriminator(\"softMuonByPtBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetMuIp"),
    quantity = cms.untracked.string("bDiscriminator(\"softMuonByIP3dBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetTKHE"),
    quantity = cms.untracked.string("bDiscriminator(\"trackCountingHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetTKHP"),
    quantity = cms.untracked.string("bDiscriminator(\"trackCountingHighPurBJetTags\") ")
    ),  
    )

#  _______      __
# |  __ \ \    / /
# | |__) \ \  / / 
# |  ___/ \ \/ /  
# | |      \  /   
# |_|       \/    

PrimaryVerticesNtuplizer = cms.EDProducer(
    "CandViewNtpProducer",
    src = cms.InputTag("offlinePrimaryVertices"),
    lazyParser = cms.untracked.bool(True),
    prefix = cms.untracked.string("PV"),
    eventInfo = cms.untracked.bool(True),
    variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("x"),
    quantity = cms.untracked.string("x")
    ),
    cms.PSet(
    tag = cms.untracked.string("y"),
    quantity = cms.untracked.string("y")
    ),
    cms.PSet(
    tag = cms.untracked.string("z"),
    quantity = cms.untracked.string("z")
    ),
    cms.PSet(
    tag = cms.untracked.string("xError"),
    quantity = cms.untracked.string("xError")
    ),
    cms.PSet(
    tag = cms.untracked.string("yError"),
    quantity = cms.untracked.string("yError")
    ),
    cms.PSet(
    tag = cms.untracked.string("zError"),
    quantity = cms.untracked.string("zError")
    ),
    cms.PSet(
    tag = cms.untracked.string("Chi2"),
    quantity = cms.untracked.string("chi2")
    ),
    cms.PSet(
    tag = cms.untracked.string("ndof"),
    quantity = cms.untracked.string("ndof")
    ),
    cms.PSet(
    tag = cms.untracked.string("isFake"),
    quantity = cms.untracked.string("isFake")
    ),
    )
    )

#  ____   _____ 
# |  _ \ / ____|
# | |_) | (___  
# |  _ < \___ \ 
# | |_) |____) |
# |____/|_____/ 
              
BeamSpotNtuplizer = cms.EDProducer(
    "CandViewNtpProducer",
    src = cms.InputTag("offlineBeamSpot"),
    lazyParser = cms.untracked.bool(True),
    prefix = cms.untracked.string("BS"),
    eventInfo = cms.untracked.bool(True),
    variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("x"),
    quantity = cms.untracked.string("x0")
    ),
    cms.PSet(
    tag = cms.untracked.string("y"),
    quantity = cms.untracked.string("y0")
    ),
    cms.PSet(
    tag = cms.untracked.string("z"),
    quantity = cms.untracked.string("z0")
    ),
    cms.PSet(
    tag = cms.untracked.string("xError"),
    quantity = cms.untracked.string("x0Error")
    ),
    cms.PSet(
    tag = cms.untracked.string("yError"),
    quantity = cms.untracked.string("y0Error")
    ),
    cms.PSet(
    tag = cms.untracked.string("zError"),
    quantity = cms.untracked.string("z0Error")
    ),
    cms.PSet(
    tag = cms.untracked.string("type"),
    quantity = cms.untracked.string("type")
    ),
    )
    )


#  __  __ ______ _______ 
# |  \/  |  ____|__   __|
# | \  / | |__     | |   
# | |\/| |  __|    | |   
# | |  | | |____   | |   
# |_|  |_|______|  |_|   

METNtuplizer =  cms.EDProducer(
    "CandViewNtpProducer",
    src=cms.InputTag("pfMet"),
    lazyParser=cms.untracked.bool(True),
    prefix=cms.untracked.string("pfMet"),
    eventInfo=cms.untracked.bool(False),
    variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("Pt"),
    quantity = cms.untracked.string("pt")
    ),
    )
    )

#   _____                _                 
#  |  __ \              | |                
#  | |__) |_ __ ___   __| |_   _  ___  ___ 
#  |  ___/| '__/ _ \ / _` | | | |/ __|/ _ \
#  | |    | | | (_) | (_| | |_| | (__|  __/
#  |_|    |_|  \___/ \__,_|\__,_|\___|\___|
                                         
ZMMNtuplizer = copy.deepcopy(zed)
ZMMNtuplizer.variables += zmm

ZEENtuplizer = copy.deepcopy(zed)
ZEENtuplizer.src = cms.InputTag("zEECand")
ZEENtuplizer.prefix = cms.untracked.string("Zee")
ZEENtuplizer.variables += zee

JetNtuplizer = copy.deepcopy(zed)
JetNtuplizer.src = cms.InputTag("cleanPatJets")
JetNtuplizer.prefix = cms.untracked.string("Jet")
JetNtuplizer.variables += jets

edmNtuplesOut = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('ZJJCandNtuples.root'),
    outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_ZMMNtuplizer_*_*",
    "keep *_ZEENtuplizer_*_*",
    "keep *_JetNtuplizer_*_*",
  #  "keep *_PrimaryVerticesNtuplizer_*_*",
  #  "keep *_BeamSpotNtuplizer_*_*",
    "keep *_TriggerResults_*_*",
    "keep *_offlinePrimaryVertices_*_*",
    "keep *_offlineBeamSpot_*_*",
    "keep *_METNtuplizer_*_*"
    )
    )
