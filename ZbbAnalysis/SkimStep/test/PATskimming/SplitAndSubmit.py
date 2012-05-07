#!/usr/bin/env python

import os
from optparse import OptionParser

#############################################################
def SplitByNumberOfLines(inputFile,maxNumLines,outputName,outputDir):
#############################################################
    output_header="""import FWCore.ParameterSet.Config as cms
zbbPAT = cms.untracked.string(\'\')
readFiles = cms.untracked.vstring()
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
            lineOut = "zbbPAT = '"+outputPATName+"'"
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
    lineOut = "zbbPAT = '"+outputPATName+"'"
    openOutputFile.write(lineOut)      
    openOutputFile.close()

    openInputFile.close()
    return outPart

############################################
def createTheLSFFile(job_name,cfg_dir,lsf_dir,inputCfgName,OUTDIR):
############################################

    Input_CMSSW_PATH    = os.environ.get('CMSSW_BASE')
    
    output_LSF_name=job_name+".lsf"
    fout=open(os.path.join(lsf_dir,output_LSF_name),'w')
    
    log_dir = os.path.join(cfg_dir,"log")
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
        
    fout.write("#!/bin/sh \n") 
    fout.write("#BSUB -L /bin/sh\n")       
    fout.write("#BSUB -J "+job_name+"\n")
    fout.write("#BSUB -o "+os.path.join(log_dir,job_name+".log")+"\n")
    fout.write("#BSUB -q 8nh\n")
    fout.write("JobName="+job_name+" \n")
    fout.write("OUT_DIR="+OUTDIR+" \n")
    fout.write("LXBATCH_DIR=`pwd` \n") 
    fout.write("cd "+os.path.join(Input_CMSSW_PATH,"src")+" \n")
    fout.write("eval `scram runtime -sh` \n")
    fout.write("cd $LXBATCH_DIR \n") 
    fout.write("cmsRun "+os.path.join(cfg_dir,inputCfgName)+" \n")
    fout.write("ls -lh . \n")
    fout.write("for RootOutputFile in $(ls *root ); do rfcp ${RootOutputFile}  ${OUT_DIR}/${RootOutputFile} ; done \n")
    
    fout.close()

##############################################
def submit(job_id,job_name,lsf_dir):
##############################################    
        print "submit job", job_id
        the_jobname = job_name+".lsf"
        submitcommand1 = "chmod u+x " + os.path.join(lsf_dir,the_jobname)
        child1  = os.system(submitcommand1)
        submitcommand2 = "bsub < "+os.path.join(lsf_dir,the_jobname)
        child2  = os.system(submitcommand2)

##############################################
def main():
##############################################
    desc="""This is a description of %prog.""" 
    parser = OptionParser(description=desc,version='%prog version 0.1')
    parser.add_option('-f','--file',help='file list (plain text)', dest='filelist', action='store')
    parser.add_option('-j','--jobname',help='additional job name, eg Mu2010A', dest='jobname', action='store')
    parser.add_option('-s','--submit',  help='job submitted', dest='submit', action='store_true', default=False)
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

    PATskimming_python_dir = os.path.join(Input_CMSSW_PATH,"src/ZbbAnalysis/SkimStep/python")
    if not os.path.exists(PATskimming_python_dir):
        os.makedirs(PATskimming_python_dir)

    PATskimming_test_dir = os.path.join(PATskimming_dir,"test")
    if not os.path.exists(PATskimming_test_dir):
        os.makedirs(PATskimming_test_dir)

    # directory to store the LSF to be submitted
    LSF_dir = os.path.join(PATskimming_test_dir,"LSF")
    if not os.path.exists(LSF_dir):
        os.makedirs(LSF_dir)     

    Template_ConfigFile = os.path.join(PATskimming_dir,"mergeEDMfiles_cfg.py")

# split section
    if opts.jobname is None:
        Output_BaseName = 'zbbPATSkim_Merge_N'
    else:
        Output_BaseName = 'zbbPATSkim_Merge_'+opts.jobname+'_N'

    OUTDIR = "/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2"

    NumOfJobs = SplitByNumberOfLines(input_list,6,Output_BaseName,PATskimming_python_dir)    
    for aJob in range(NumOfJobs):
        fin = open(Template_ConfigFile)
        cfg_list = []
        cfg_list.append("import FWCore.ParameterSet.Config as cms \n")
        cfg_list.append("process = cms.Process(\"MergeEDMFiles\")\n")
        cfg_list.append("import ZbbAnalysis.SkimStep."+Output_BaseName.replace("N",str(aJob+1))+"_cff as myfiles \n")
        cfg_list+=fin.readlines()

# write the cfg file 
        outputCfgName=Output_BaseName.replace("N",str(aJob+1))+"_cfg.py"
        outputJobName=Output_BaseName.replace("N",str(aJob+1))
        theConfigFile = os.path.join(PATskimming_test_dir,outputCfgName)
        fout = open(theConfigFile, "w")
        fout.writelines(cfg_list)

        createTheLSFFile(outputJobName,PATskimming_test_dir,LSF_dir,outputCfgName,OUTDIR)
        if opts.submit:
            submit(aJob,outputJobName,LSF_dir)

#Output_EDA = 

#    bsubfile.append("bsub" % (queue, director, jobnumber, waiter, jobnumber))

#
if __name__ == "__main__":        
    main()
