#!/usr/bin/env python
# TODO: 
# 1. pass the number of files in each group as argument

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

    def __init__(self, job_id, channel, isData, selection, sample, bsel, pusel, tnpsel,channelsel, CMSSW_dir ,the_dir):
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
        self.channelsel=channelsel

        self.the_dir=the_dir
        self.CMSSW_dir=CMSSW_dir

        self.output_full_name=self.getOutputBaseName()+"_"+str(self.job_id)

        self.cfg_dir=None
        self.outputCfgName=None
# LSF variables        
        self.LSF_dir=None
        self.output_LSF_name=None

        self.lfn_list=list()      

    def __del__(self):
###############################
        del self.lfn_list

    def setCASTORout(self,theCASTORdir):    
###############################
        self.OUTDIR = theCASTORdir
          
    def getOutputBaseName(self):
########################    
        return "acceptanceCalculator_"+self.job_name
        
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
# FIXME: in case DATA  is used add a WARNING and EXIT!
        else:
            template_cfg_file = os.path.join(self.the_dir,"AcceptanceCalculator_MC_TEMPL_cfg.py")

        fin = open(template_cfg_file)
        for line in fin.readlines():
            if self.isData:
                if line.find("OUTLISTTEMPLATE")!=-1:
                    line=line.replace("OUTLISTTEMPLATE",self.job_name+"_"+str(self.job_id))
            else:
                if line.find("OUTLISTTEMPLATE")!=-1:
                    line=line.replace("OUTLISTTEMPLATE",self.job_name+"_"+str(self.job_id))
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
            if line.find("CHANNELTEMPLATE")!=-1:
                line=line.replace("CHANNELTEMPLATE",self.channelsel)    
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
    t = now.strftime("test_%Y_%m_%d_%H_%M_%S_MC")
    t+=opts.taskname
    castordir = os.path.join("/castor/cern.ch/user/m/musich/ZbbAnalysis/zbbAcceptance",t)
    mkdir = "rfmkdir "+ castordir
    osmkdir  = os.system(mkdir)
    rfchmod = "rfchmod +775 "+ castordir
    osrfchmod = os.system(rfchmod)
    

    zbbSrc_SHERPA_5FS=[
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_1_3_KkT.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_2_3_JT1.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_3_3_cnb.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_4_2_lKb.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_5_3_0PV.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_6_1_3qS.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_7_1_Tkk.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_8_3_qeG.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_9_3_rqB.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_10_2_179.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_11_3_gDi.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_12_2_rdZ.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_13_3_BmO.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_14_3_OB6.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_15_2_mON.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_16_2_GKJ.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_17_3_nUD.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_18_2_D8S.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_19_3_s4T.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_20_2_sSs.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_21_2_A4l.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_22_3_pDS.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_23_2_Dbl.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_24_3_H9s.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_25_1_vvi.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_26_1_xGf.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_27_2_PRJ.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_28_2_xSE.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_29_1_P0y.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_30_2_GAI.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_31_2_F5k.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_32_2_ObC.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_33_1_Esr.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_34_2_q7S.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_35_2_hRQ.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_36_2_rti.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_37_2_iht.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched_forAcceptanceStudies/MergedOutputFile_38_1_Vk4.root'
        ]

    zbbSrc_MADGRAPH_5FS=[
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_1_2_PX9.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_2_2_1yW.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_3_1_w2Z.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_4_1_qIM.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_5_1_oz4.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_6_2_Fv1.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_7_2_oVO.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_8_1_1c1.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_9_2_S0I.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_10_2_IOo.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_11_2_zfe.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_12_2_0i3.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_13_2_oFF.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_14_2_WxD.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_15_2_V9P.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_16_2_ySu.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_17_1_Et3.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_18_2_jtw.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_19_1_jVD.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_20_2_IR4.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_21_2_Tu2.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_22_2_dGJ.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_23_1_pmz.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_24_2_HqO.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_25_1_hxw.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_26_1_bgq.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_27_1_GPa.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_28_1_CR9.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_29_1_cCs.root',     
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_30_1_k7E.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_31_1_tqe.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_32_2_WpE.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_33_1_nMC.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_34_1_uBP.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_35_1_UZE.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_36_1_rbC.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_37_1_dqJ.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_38_1_MNR.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched_forAcceptanceStudies_V2/MergedOutputFile_39_1_0o2.root'
        ]

    zbbSrc_aMCATNLO_5FS=[
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_1_1_rKV.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_2_1_mWo.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_3_1_8Wp.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_4_1_k6w.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_5_1_XfD.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_6_1_cbz.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_7_1_Kk4.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_8_1_rnF.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_9_1_Ld9.root',      
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_10_1_GM3.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_11_1_luU.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_12_1_4dS.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_13_1_YZB.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_14_1_5S6.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_15_1_Lri.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_16_1_HJ8.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_17_3_4bl.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_18_1_IiM.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_19_1_vkw.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_20_1_Pa2.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_21_1_3if.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_22_3_VME.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_23_1_uE6.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_24_3_C5Y.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_25_1_dPO.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_26_1_DRt.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_27_1_rYT.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_28_3_moU.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_29_3_oTS.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_30_3_l6d.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_31_3_Wvo.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_32_1_xeR.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_33_1_bDC.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_34_1_kwV.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_35_3_6a5.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_36_1_091.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_37_1_1g9.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_38_3_Zfb.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_39_3_Kj4.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_40_1_R9N.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_41_1_hQS.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_42_1_o4V.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_43_3_i9j.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_44_1_STU.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_45_1_j0U.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_46_1_omL.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_47_1_TO0.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_48_1_TXs.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_49_1_MG6.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_50_1_OU8.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_51_1_0YC.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_52_1_Vey.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_53_1_s5q.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_54_1_Ctt.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_55_1_kwD.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_56_1_Sjo.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_57_1_arq.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_58_1_2oT.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_59_1_VzW.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_60_1_0Cy.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_61_1_h4c.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO_forAcceptanceStudies/MergedOutputFile_62_1_L6h.root'
        ]

    ZjetsInclSrc = [
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_100_1_mjJ.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_104_1_CSZ.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_109_1_lcq.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_111_1_4rL.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_113_5_e0d.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_117_5_syP.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_118_5_Ifc.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_11_1_P45.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_121_1_wWm.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_122_1_0Ln.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_128_1_ODY.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_130_4_Mib.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_146_4_2pN.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_151_2_e6d.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_158_4_Ewk.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_161_6_ira.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_162_4_LzI.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_1_1_ug7.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_28_4_5aJ.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_2_1_Bw5.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_39_4_hQG.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_48_5_VjS.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_49_5_5ii.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_50_5_ysU.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_51_5_ycw.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_53_5_W9Z.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_60_2_mdT.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_72_5_Fbk.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_73_5_TeW.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_76_5_GjH.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_85_1_yQv.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_8_4_ksk.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_forAcceptanceStudies_V2/MergedOutputFile_93_1_ZJU.root'
        ]
    
# assign the  datasets
    zbbSrc         = zbbSrc_MADGRAPH_5FS
    zbbSherpaSrc   = zbbSrc_SHERPA_5FS
    zbbaMCatNLOSrc = zbbSrc_aMCATNLO_5FS

    channels  = ["MC_ZeeInclusive_PU","MC_Zeeb5fToLLSherpa_PU","MC_Zeeb5fToLL_PU","MC_ZmmInclusive_PU","MC_Zmmb5fToLLSherpa_PU","MC_Zmmb5fToLL_PU"]

    srcFiles = [ZjetsInclSrc,zbbSherpaSrc,zbbSrc,ZjetsInclSrc,zbbSherpaSrc,zbbSrc]

    isData =   [False,False,False,False,False,False]
    
    selections  = ["Sample","Sample","Sample","Sample","Sample","Sample"]
 
    samples     = ["ZIncl","Sherpa","MadGraph","ZIncl","Sherpa","MadGraph"]

    pusels      = ["True","True","True","True","True","True"]

    tnpsels     = ["False","False","False","False","False","False"]
    
    bsels       = ["False","False","False","False","False","False"]

    decays      = ["Z_e","Z_e","Z_e","Z_m","Z_m","Z_m"]
    
# CMSSW section
    input_CMSSW_BASE = os.environ.get('CMSSW_BASE')
    AnalysisStep_dir = os.path.join(input_CMSSW_BASE,"src/ZbbAnalysis/AnalysisStep/test/AcceptanceStudies/")

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

            
        for jobN,theSrcFiles in enumerate(split(srcFiles[iChannel],5)):            
            aJob = Job(jobN,channels[iChannel],isData[iChannel],selections[iChannel],samples[iChannel],bsels[iChannel],pusels[iChannel],tnpsels[iChannel],decays[iChannel],input_CMSSW_BASE,AnalysisStep_dir)
            aJob.setCASTORout(castordir)
            aJob.createTheCfgFile(theSrcFiles)
            aJob.createTheLSFFile()

            output_file_list1.append("rfcp "+aJob.getOutputFileName()+" . \n")
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


   
