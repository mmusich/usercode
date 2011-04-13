#!/usr/bin/env python

# generic  python modules
from optparse import OptionParser
  
# CMSSW modules
from DataFormats.FWLite import Events, Lumis

# local modules
from fileList_Muon2010 import input_files
from zbbPlots import zbbPlots

def main():
#########################
# optparse for pedestrian at http://www.alexonlinux.com/pythons-optparse-for-human-beings
 
    desc="""This is a description of %prog.""" 
    parser = OptionParser(description=desc,version='%prog version 0.1')
    parser.add_option('--mc',  help='is Monte Carlo',  dest='isMC',default=False, action='store_true')
    (opts, args) = parser.parse_args()

    zbb_plots = zbbPlots(opts.isMC)
    zbb_plots.beginJob()
    
    events     = Events(input_files)
    # loop over events
    for i, event in enumerate(events):
        if i%1000==0:
            print "zbbMain processing evt: ", i
            print "runNumber: ", event.eventAuxiliary().run()
     
        zbb_plots.analyzeEvent(event)

    lumiBlocks = Lumis(input_files)
    # loop over lumi blocks
    for i, lumi in enumerate(lumiBlocks):
        if i%100==0: print "zbbMain processing lumi block ", i
        zbb_plots.analyzeLumiBlock(lumi)

    zbb_plots.endJob()        
    

if __name__ == "__main__":
    main()

