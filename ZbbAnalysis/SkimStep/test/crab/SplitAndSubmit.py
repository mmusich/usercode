#!/usr/bin/env python

import os
from optparse import OptionParser

#############################################################
def SplitByNumberOfLines(inputFile,maxNumLines,outputName,outputDir):
#############################################################
    output_header="""import FWCore.ParameterSet.Config as cms
readFiles.extend( [
"""

    openInputFile = open(inputFile, 'r')

    numLines=0
    outPart=1    
    outputCffName=outputName.replace("N",str(outPart))+"_cff.py"
    openOutputFile=open(os.path.join(outputDir,outputCffName),'w')
    openOutputFile.write(output_header)

    localFilePath="rfio:"

# read-in file 
    while 1:
        line = openInputFile.readline()
        if not line:
            break # end of input file
        numLines+=1
        if numLines < maxNumLines:
            lineOut="'"+localFilePath+line[:-1]+"'"
            if not numLines == (maxNumLines-1):
                lineOut+=",\n"
            openOutputFile.write(lineOut)            
        else:
# close previous file
            openOutputFile.write("]) \n")
            outputPATName=outputName.replace("N",str(outPart))+".root"
            lineOut = "zbbPAT = "+outputPATName
            openOutputFile.write(lineOut)            
            openOutputFile.close()

# open next file
            outPart+=1
            outputCffName=outputName.replace("N",str(outPart))+"_cff.py"
            openOutputFile=open(os.path.join(outputDir,outputCffName),'w')
            openOutputFile.write(output_header)
# reset counter
            numLines=0    

# close previous file
    openOutputFile.write("]) \n")
    outputPATName=outputName.replace("N",str(outPart))+".root"
    lineOut = "zbbPAT = "+outputPATName
    openOutputFile.write(lineOut)      
    openOutputFile.close()

    openInputFile.close()
    return outPart


##############################################
def main():
##############################################
    desc="""This is a description of %prog.""" 
    parser = OptionParser(description=desc,version='%prog version 0.1')
    parser.add_option('-f','--file',help='file list (plain text)', dest='filelist', action='store')
    parser.add_option('-j','--jobname',help='additional job name, eg Mu2010A', dest='jobname', action='store')
    (opts, args) = parser.parse_args()

    if opts.filelist is None:
        print "-f not specified\n"
        parser.print_help()
        exit(-1)
    else:                
        input_list = opts.filelist


# CMSSW section
    Input_CMSSW_PATH    = os.environ.get('CMSSW_BASE')

# directory to store the cff and the cfg
    PATskimming_dir = os.path.join(Input_CMSSW_PATH,"src/ZbbAnalysis/SkimStep/test/PATskimming")

    PATskimming_python_dir = os.path.join(PATskimming_dir,"python")
    if not os.path.exists(PATskimming_python_dir):
        os.makedirs(PATskimming_python_dir)

    PATskimming_test_dir = os.path.join(PATskimming_dir,"test")
    if not os.path.exists(PATskimming_test_dir):
        os.makedirs(PATskimming_test_dir)

    Template_ConfigFile = os.path.join(PATskimming_dir,"mergeEDMfiles_cfg.py")

# split section
    if opts.jobname is None:
        Output_BaseName = 'zbbPATSkim_Merge_N'
    else:
        Output_BaseName = 'zbbPATSkim_Merge_'+opts.jobname+'_N'

    NumOfJobs = SplitByNumberOfLines(input_list,255,Output_BaseName,PATskimming_python_dir)    
    for aJob in range(NumOfJobs):
        fin = open(Template_ConfigFile)
        cfg_list = fin.readlines()
        cfg_list.append("process.load(\"ZbbAnalysis.SkimStep.test.PATskimming."+Output_BaseName.replace("N",str(aJob+1))+"_cff\")")
        
# write the cfg file 
        outputCfgName=Output_BaseName.replace("N",str(aJob+1))+".py"
        theConfigFile = os.path.join(PATskimming_test_dir,outputCfgName)
        fout = open(theConfigFile, "w")
        fout.writelines(cfg_list)





    
#Output_EDA = 

#    bsubfile.append("bsub" % (queue, director, jobnumber, waiter, jobnumber))

#
if __name__ == "__main__":        
    main()
