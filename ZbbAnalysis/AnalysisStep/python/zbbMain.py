#!/usr/bin/env python

# generic  python modules
from optparse import OptionParser
import inspect
  
# CMSSW modules
from DataFormats.FWLite import Events, Lumis
from PhysicsTools.PythonAnalysis import *
from ROOT import *

# prepare the FWLite autoloading mechanism
gSystem.Load("libFWCoreFWLite.so")
AutoLibraryLoader.enable()

# enable support for files > 2 GB
#gSystem.Load("libIOPoolTFileAdaptor")
#ui = TFileAdaptorUI()

# local modules
from fileList_Muon2010 import input_files
from zbbPlots import zbbPlots

# load analyzer module
ROOT.gSystem.Load('libZbbAnalysisAnalysisStep.so')

def main():
#########################
# optparse for pedestrian at http://www.alexonlinux.com/pythons-optparse-for-human-beings
 
    desc="""This is a description of %prog.""" 
    parser = OptionParser(description=desc,version='%prog version 0.1')
    parser.add_option('--mc',  help='is Monte Carlo',  dest='isMC',default=False, action='store_true')
    (opts, args) = parser.parse_args()

    zbb_plots = zbbPlots(opts.isMC)
    zbb_plots.beginJob()
    
    patsimpleanalyzer = ROOT.patba.PatSimpleAnalyzer()
    
    #print inspect.getmembers(patsimpleanalyzer)
    patsimpleanalyzer.beginJob()
    
    events     = Events(input_files)
    # loop over events
    for i, event in enumerate(events):
        if i%1000==0:
            print "zbbMain processing evt: ", i
            print "runNumber: ", event.eventAuxiliary().run()
     
        zbb_plots.analyzeEvent(event)
        patsimpleanalyzer.analyzeEventOnly(event.object())

    lumiBlocks = Lumis(input_files)
    # loop over lumi blocks
    for i, lumi in enumerate(lumiBlocks):
        if i%100==0: print "zbbMain processing lumi block ", i
        zbb_plots.analyzeLumiBlock(lumi)

    patsimpleanalyzer.endJob()
    zbb_plots.endJob()        

if __name__ == "__main__":
    main()

