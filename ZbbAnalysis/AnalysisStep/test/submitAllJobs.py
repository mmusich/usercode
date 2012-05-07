#!/usr/bin/env python

import os,sys
import string, re
# generic  python modules
from optparse import OptionParser
from time import gmtime, localtime, strftime

desc="""This is a description of %prog."""
parser = OptionParser(description=desc,version='%prog version 0.1')
parser.add_option('-s','--submit',  help='job submitted', dest='submit', action='store_true', default=False)
(opts, args) = parser.parse_args()

#########################
# matrix of needed imputs
#########################

channels  = ["DATA2010_electron","DATA2010_muon","DATA2011_electron","DATA2011_muon","MC_ZbToLL","MC_ZcToLL","MC_ZlJets","MC_TTJets"]
selection = ["isMuon","isMuon","isMuon","isMuon","Sample","Sample","Sample","Sample"]
samples   = ["False","True","False","True","Zb","Zc","Zl","tt"]

#########################
# modify the template 
#########################
def produceLSFTemplateFile(channel,index):
    fin = open("ZbbAnalysis_TEMPL.lsf")
    if "MC_" in channel:
        pset_cfg = "ZbbEventContentAnalyzer_MC_2011_cfg.py"
    else :
        if "2010" in channel :
            pset_cfg = "ZbbEventContentAnalyzer_DATA_2010_cfg.py"
        else :
            pset_cfg = "ZbbEventContentAnalyzer_DATA_2011_cfg.py" 
    pset_lsf      = "ZbbAnalysis_" + channel + ".lsf"
    outfile_root  = "analyzePAT_"  + channel + ".root"
    fout = open(pset_lsf,"w")
    for line in fin.readlines():
        if  line.find("JOBNAMETEMPLATE")!=-1:
            line=line.replace("JOBNAMETEMPLATE",channel)
        if line.find("OUTFILETEMPLATE")!=-1:
            line=line.replace("OUTFILETEMPLATE",outfile_root)    
        if  line.find("CFGTEMPLATE")!=-1:
            line=line.replace("CFGTEMPLATE",pset_cfg)
        if  line.find("SELTEMPLATE")!=-1:
            line=line.replace("SELTEMPLATE",selection[index])
        if  line.find("SAMPLETEMPLATE")!=-1:
            line=line.replace("SAMPLETEMPLATE",samples[index])    
        fout.write(line)             
    print pset_lsf + " has been written.\n"

#########################
# produce lsf files
#########################
for i in range(len(channels)):
    produceLSFTemplateFile(channels[i],i)

for i in range(len(channels)):
    submitcommand1 = "chmod +x " +  "ZbbAnalysis_" + channels[i] + ".lsf"
    child1  = os.system(submitcommand1)
    if opts.submit:
        submitcommand2 = "bsub -o tmp.tmp -q 8nh " +  "ZbbAnalysis_" + channels[i] + ".lsf "
        child2  = os.system(submitcommand2)
    else :
        print "job not submitted -- running in dry mode"


