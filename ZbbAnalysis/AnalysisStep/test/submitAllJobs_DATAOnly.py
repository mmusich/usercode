#!/usr/bin/env python

import datetime,time
import os,sys
import string, re
# generic  python modules
from optparse import OptionParser

def split(sequence, size):
##########################    
# aux generator function to split lists
# based on http://sandrotosi.blogspot.com/2011/04/python-group-list-in-sub-lists-of-n.html
# about generators see also http://stackoverflow.com/questions/231767/the-python-yield-keyword-explained
##########################
    for i in xrange(0, len(sequence), size):
        yield sequence[i:i+size] 

#############
class Job:
#############

    def __init__(self, job_id, channel, isData, selection, sample, bsel, pusel, tnpsel, CMSSW_dir ,the_dir):
##############################################################
        self.job_id=job_id
# FIXME add comments!        
        self.job_name=channel
        self.isData=isData
        self.selection=selection
        self.sample=sample
        self.bsel= bsel
        self.pusel=pusel
        self.tnpsel=tnpsel        

        self.the_dir=the_dir
        self.CMSSW_dir=CMSSW_dir

        self.output_full_name=self.getOutputBaseName()+"_"+str(self.job_id)

        self.cfg_dir=None
        self.outputCfgName=None
# LSF variables        
        self.LSF_dir=None
        self.output_LSF_name=None

        self.lfn_list=list()      

        #self.OUTDIR = "/castor/cern.ch/user/m/musich/ZbbAnalysis/test01Sept" # TODO: write a setter method
        #self.OUTDIR = self.createCASTORout()

    def __del__(self):
###############################
        del self.lfn_list

    def setCASTORout(self,theCASTORdir):    
###############################
        self.OUTDIR = theCASTORdir
          
    def getOutputBaseName(self):
########################    
        return "analyzePAT_"+self.job_name
        
    def createTheCfgFile(self,lfn):
###############################
        
# write the cfg file 
        self.cfg_dir = os.path.join(self.the_dir,"cfg")
        if not os.path.exists(self.cfg_dir):
            os.makedirs(self.cfg_dir)

        self.outputCfgName=self.output_full_name+"_cfg.py"
        fout=open(os.path.join(self.cfg_dir,self.outputCfgName),'w')

# decide which template according to data/mc
        if self.isData:
            template_cfg_file = os.path.join(self.the_dir,"ZbbEventContentAnalyzer_DATA_TEMPL_cfg.py")
        else:
            template_cfg_file = os.path.join(self.the_dir,"ZbbEventContentAnalyzer_MC_TEMPL_cfg.py")

        fin = open(template_cfg_file)
        for line in fin.readlines():
            if self.isData:
                if line.find("OUTLISTTEMPLATE")!=-1:
                    line=line.replace("OUTLISTTEMPLATE",self.job_name+"_"+str(self.job_id))
            else:                    
                if line.find("PUTEMPLATE")!=-1:
                    line=line.replace("PUTEMPLATE",self.pusel)
                if line.find("TNPTEMPLATE")!=-1:
                    line=line.replace("TNPTEMPLATE",self.tnpsel)
                if line.find("BSELTEMPLATE")!=-1:
                    line=line.replace("BSELTEMPLATE",self.bsel)    
                    
            if line.find("FILESOURCETEMPLATE")!=-1:
                lfn_with_quotes = map(lambda x: "\'"+x+"\'",lfn)                   
                #print "["+",".join(lfn_with_quotes)+"]"
                line=line.replace("FILESOURCETEMPLATE","["+",".join(lfn_with_quotes)+"]") 
            if line.find("OUTFILETEMPLATE")!=-1:
                line=line.replace("OUTFILETEMPLATE",self.output_full_name+".root")     
            fout.write(line)    

        fout.close()

    def createTheLSFFile(self):
###############################

# directory to store the LSF to be submitted
        self.LSF_dir = os.path.join(self.the_dir,"LSF")
        if not os.path.exists(self.LSF_dir):
            os.makedirs(self.LSF_dir)

        self.output_LSF_name=self.output_full_name+".lsf"
        fout=open(os.path.join(self.LSF_dir,self.output_LSF_name),'w')
    
        job_name = self.output_full_name

        log_dir = os.path.join(self.the_dir,"log")
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)

        fout.write("#!/bin/sh \n") 
        fout.write("#BSUB -L /bin/sh\n")       
        fout.write("#BSUB -J "+job_name+"\n")
        fout.write("#BSUB -o "+os.path.join(log_dir,job_name+".log")+"\n")
        fout.write("#BSUB -q 8nh\n")
        fout.write("JobName="+job_name+" \n")
        fout.write("OUT_DIR="+self.OUTDIR+" \n")
        fout.write("LXBATCH_DIR=`pwd` \n") 
        fout.write("cd "+os.path.join(self.CMSSW_dir,"src")+" \n")
        fout.write("eval `scram runtime -sh` \n")
        fout.write("cd $LXBATCH_DIR \n") 
        fout.write("cmsRun "+os.path.join(self.cfg_dir,self.outputCfgName)+" "+self.selection+"="+self.sample+" \n")
        fout.write("ls -lh . \n")
        fout.write("for RootOutputFile in $(ls *root ); do rfcp ${RootOutputFile}  ${OUT_DIR}/${RootOutputFile} ; done \n")
        fout.write("for TxtOutputFile in $(ls *txt ); do rfcp ${TxtOutputFile}  ${OUT_DIR}/${TxtOutputFile} ; done \n")

        fout.close()

    def getOutputFileName(self):
############################################
        return os.path.join(self.OUTDIR,self.output_full_name+".root")
        
    def submit(self):
###############################        
        print "submit job", self.job_id
        job_name = self.output_full_name
        submitcommand1 = "chmod u+x " + os.path.join(self.LSF_dir,self.output_LSF_name)
        child1  = os.system(submitcommand1)
        submitcommand2 = "bsub < "+os.path.join(self.LSF_dir,self.output_LSF_name)
        child2  = os.system(submitcommand2)

##############################################
def main():
##############################################

    desc="""This is a description of %prog."""
    parser = OptionParser(description=desc,version='%prog version 0.1')
    parser.add_option('-s','--submit',  help='job submitted', dest='submit', action='store_true', default=False)
    parser.add_option('-j','--jobname', help='task name', dest='taskname', action='store', default='')
    (opts, args) = parser.parse_args()

    now = datetime.datetime.now()
    t = now.strftime("test_%Y_%m_%d_%H_%M_%S_DATA")
    t+=opts.taskname
    castordir = os.path.join("/castor/cern.ch/user/m/musich/ZbbAnalysis",t)
    mkdir = "rfmkdir "+ castordir
    osmkdir  = os.system(mkdir)
    rfchmod = "rfchmod +775 "+ castordir
    osrfchmod = os.system(rfchmod)
    
    electronSrc=[
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_1_1_4IW.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_2_1_Uce.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_3_1_Qhk.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_4_1_Nuq.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_5_1_Khn.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_6_1_dg5.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_7_1_bJK.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_8_1_cbt.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_9_2_LrS.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_10_3_pL3.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_11_1_zYu.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_12_2_efg.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_13_1_jwC.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_14_2_UEa.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_15_2_kEJ.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_16_2_vl4.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_17_1_bNs.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_18_2_YvL.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_19_1_Vmj.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_20_1_MkM.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_21_1_pZf.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_22_1_iQD.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_23_1_K2F.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_24_1_w6G.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_25_1_X2o.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleElectron/MergedOutputFile_26_1_KuX.root'
        ] 
    
    muonSrc=[
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_1_1_Uui.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_2_1_lUC.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_3_1_iNw.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_4_1_Z8m.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_5_1_44d.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_6_1_vyd.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_7_1_717.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_8_1_6bu.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_9_1_2K8.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_10_1_lBg.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_11_1_s3U.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_12_1_SfB.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_13_2_hi7.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_14_2_7Yq.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_15_2_ffA.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_16_3_JAL.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_17_1_Osd.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_18_2_T0M.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_19_2_jlN.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_20_2_zac.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_21_1_WYP.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_22_3_3PF.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_23_1_1xg.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_24_1_1GM.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_25_2_s8x.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_26_2_ZdI.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_27_1_WbB.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_28_2_x4D.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_29_3_ld4.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_30_2_XwM.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_31_2_LRz.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_32_4_je7.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_33_3_LKF.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_34_1_KNH.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_35_2_QAN.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_36_2_lzV.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_37_2_z2G.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DoubleMu/MergedOutputFile_38_1_neq.root'
        ]
    
    channels  =   ["DATA2011_electron","DATA2011_muon"]
    srcFiles  =   [electronSrc, muonSrc]
    isData    =   [True,  True] 
    selections=   ["isMuon","isMuon"]
    samples   =   ["False", "True"]     #false = ele,  false= mu
    bsels     =   ["DUMMY", "DUMMY"]
    pusels    =   ["DUMMY", "DUMMY"]
    tnpsels   =   ["DUMMY", "DUMMY"]


# CMSSW section
    input_CMSSW_BASE = os.environ.get('CMSSW_BASE')
    AnalysisStep_dir = os.path.join(input_CMSSW_BASE,"src/ZbbAnalysis/AnalysisStep/test")

    for iChannel in range(len(channels)):

# for hadd script
        scripts_dir = os.path.join(AnalysisStep_dir,"scripts")
        if not os.path.exists(scripts_dir):
            os.makedirs(scripts_dir)
        hadd_script_file = os.path.join(scripts_dir,channels[iChannel])
        fout = open(hadd_script_file,'w')

        output_file_list1=list()      
        output_file_list2=list()
        output_file_list2.append("hadd ")

            
        for jobN,theSrcFiles in enumerate(split(srcFiles[iChannel],8)):            
            aJob = Job(jobN,channels[iChannel],isData[iChannel],selections[iChannel],samples[iChannel],bsels[iChannel],pusels[iChannel],tnpsels[iChannel],input_CMSSW_BASE,AnalysisStep_dir)
            aJob.setCASTORout(castordir)
            aJob.createTheCfgFile(theSrcFiles)
            aJob.createTheLSFFile()

            output_file_list1.append("xrdcp root://castorcms/"+aJob.getOutputFileName()+" . \n")
            if jobN == 0:
                output_file_list2.append(aJob.getOutputBaseName()+".root ")
            output_file_list2.append(os.path.split(aJob.getOutputFileName())[1]+" ")    
   
            if opts.submit:
                aJob.submit()
            del aJob

        fout.writelines(output_file_list1)
        fout.writelines(output_file_list2)

        fout.close()
        del output_file_list1
        

if __name__ == "__main__":        
    main()


   
