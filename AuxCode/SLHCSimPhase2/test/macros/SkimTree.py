#!/usr/bin/env python
import ROOT
from optparse import OptionParser

###############
def main():
###############

# command line options
    parser = OptionParser()
    parser.add_option("-f", "--file",  
                      action="store", type="string", dest="input_root_filename",
                      help="input root file")

    (options, args) = parser.parse_args()

# input
    chain = ROOT.TChain("PixelNtuple")
    chain.Add(options.input_root_filename)
    chain.GetEntry(0)
    nentries = chain.GetEntries()
    print  "+++++ No. of entries in the input tree: ", nentries

# define selection
#    string_eta_cut = "fabs(asinh(gz/sqrt(gx*gx+gy*gy)))<0.2&&subid==1&&layer==1&&sqrt(gx*gx+gy*gy)<2.9"
#    eta_cut = "_BPixL1_inner_eta_lt_02"
    string_eta_cut = "subid==1&&layer==1&&sqrt(gx*gx+gy*gy)<2.9"
    eta_cut = "_BPixL1_inner"

# output
    output_root_filename = options.input_root_filename.split('.root')[0]+eta_cut+'.root'
    print output_root_filename
    fileTreeOut = ROOT.TFile(output_root_filename,"RECREATE")

    fileTreeOut.cd()
    newtree = ROOT.TTree()
    newtree= chain.CopyTree(string_eta_cut,"",nentries)
    newtree.AutoSave()

    nentries = newtree.GetEntries()
    print "+++++ No. of entries in the output tree: ", nentries 
    
    fileTreeOut.Write()
    fileTreeOut.Close()
    
##################################
if __name__ == "__main__":        
    main()
