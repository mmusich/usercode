import ROOT

# CMSSW modules
from DataFormats.FWLite import Handle

# local modules
import eventHelperFunctions as ehf

#####################
class zbbPlots:
#####################
    
    def __init__(self,isMC=False):
#########################################
        self.root_output_file = ROOT.TFile("zbbPlots.root","RECREATE")
        self.dir = self.root_output_file.mkdir("eventSelection")

        self.isMC = isMC

        self.counterLabels = ('TotalEventCounter','AfterPVFilterCounter', 'AfterNSFilterCounter', 'AfterPATCounter', 'AfterCandidatesCounter', 'AfterJetsCounter')
        self.muonLabel="patMuonsWithTrigger"
        self.electronLabel ="patElectronsWithTrigger"
        self.metLabel="patMETs"
        self.jetLabel="cleanPatJetsPF"
        self.zllLabel="zMMCand"
        self.triggerLabel="patTriggerEvent"
        

    def beginJob(self):
#########################################
        #self.root_output_file.cd()
        self.dir.cd()

        # lumi block plots
        self.h1_eventCounters = ROOT.TH1F("EventCounter","Event Counters",6,-0.5,5.5)
        for theCounterLabel in self.counterLabels:
            self.h1_eventCounters.GetXaxis().SetBinLabel(self.counterLabels.index(theCounterLabel)+1,theCounterLabel)

        # muon plots
        self.h1_muonPt  = ROOT.TH1F("muonPt", "muon pt; p_{T} (GeV)",    100,  0.,300.)
        self.h1_muonEta = ROOT.TH1F("muonEta","muon #eta; #eta", 100, -3.,  3.)
        self.h1_muonPhi = ROOT.TH1F("muonPhi","muon #phi; #phi (rad)", 100, -3.14, 3.14)


        # trigger plots

        self.h1_triggerSelection = ROOT.TH1F("triggerSelection","triggerSelection",2,0,2)
        self.h1_triggerBit = ROOT.TH1F("triggerBits","trigger bits",20,0,20)
        for theTriggerLabel in ehf.triggers:
            self.h1_triggerBit.GetXaxis().SetBinLabel(ehf.triggers.index(theTriggerLabel)+1,theTriggerLabel)

    
        # create handle outside of loop        
        self.counterHandle      = Handle("edm::MergeableCounter")
        self.muonHandle         = Handle("std::vector<pat::Muon>")
        self.electronHandle     = Handle("std::vector<pat::Electron>")
        self.jetHandle          = Handle("std::vector<pat::Jet>")
        self.metHandle          = Handle("std::vector<pat::MET>")
        self.zllCandidateHandle = Handle("std::vector<reco::CompositeCandidate>")
        self.triggerEventHandle = Handle("pat::TriggerEvent")

    def analyzeEvent(self,event):
#########################################        
        # get the muon product
        try:
            event.getByLabel (self.muonLabel, self.muonHandle)
            if self.muonHandle.isValid():
                muonsCollection = self.muonHandle.product()
                for theMuon in muonsCollection :
                    self.h1_muonPt.Fill( theMuon.pt() )
                    self.h1_muonEta.Fill( theMuon.eta() )
                    self.h1_muonPhi.Fill( theMuon.phi() )
        except:
            print "Muon not found"
            
        # get the electron product
        try:
            event.getByLabel (self.electronLabel, self.electronHandle)
            if self.electronHandle.isValid():
                electronsCollection =  self.electronHandle.product()
        except:
            print "Electron not found"
                
        # get the jet product
        try:
            event.getByLabel (self.jetLabel, self.jetHandle)
            if self.jetHandle.isValid():
                jetCollection = self.jetHandle.product()
        except:
            print "Jet not found"

        # get the met product
        try:
            event.getByLabel (self.metLabel, self.metHandle)
            if self.metHandle.isValid():
                metCollection = self.metHandle.product()
        except:
            print "MET not found"

        # get the zll candidate product
        try:
            event.getByLabel (self.zllLabel, self.zllCandidateHandle)
            if self.zllCandidateHandle.isValid():
                zllCandidatesCollection = self.zllCandidateHandle.product()
        except:
            print "Zll not found"

        # get the trigger info product
        try:
            event.getByLabel (self.triggerLabel, self.triggerEventHandle)
            if self.triggerEventHandle.isValid():                 
                triggerEvent = self.triggerEventHandle.product()
                selTriggers = ehf.selectedTriggers(triggerEvent,ehf.triggers)
                for trigger,triggered in enumerate(selTriggers):
                    if triggered : self.h1_triggerBit.Fill(trigger)

                #### list triggers fired    
                #ehf.listTheTriggers(triggerEvent)    

                self.h1_triggerSelection.Fill(ehf.isTriggerOK(triggerEvent,True,event.eventAuxiliary().run()))
                #self.h1_triggerSelection.Fill(ehf.isTriggerOK(triggerEvent))

        except:        
            print "Trigger not found"


    def analyzeLumiBlock(self,lumi_block):
#########################################
        for theCounterLabel in self.counterLabels:
            try:
                lumi_block.getByLabel (theCounterLabel, self.counterHandle)
                if self.counterHandle.isValid():
                    self.h1_eventCounters.AddBinContent(self.counterLabels.index(theCounterLabel)+1,self.counterHandle.product().value)
            except:
                print "Lumi block not found"
            
    def endJob(self):
#########################################        
        self.root_output_file.cd()
        self.root_output_file.Write()
        if not self.root_output_file is None:
            self.root_output_file.Close()
