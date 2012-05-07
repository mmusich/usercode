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
    castordir = os.path.join("/castor/cern.ch/user/m/musich/ZbbAnalysis",t)
    mkdir = "rfmkdir "+ castordir
    osmkdir  = os.system(mkdir)
    rfchmod = "rfchmod +775 "+ castordir
    osrfchmod = os.system(rfchmod)
    
    dySrc_MADGRAPH_5FS=[
        ## ok for inclusive studies
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_1_1_lot.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_2_1_epd.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_3_1_2XE.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_4_1_71x.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_5_1_x1F.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_6_1_VYR.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_7_1_spv.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_8_1_uuy.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_9_1_Z8W.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_10_1_2R9.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_11_2_UeP.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_12_1_Bz5.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_13_1_bmU.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_14_1_tYE.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_15_1_dWa.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_16_1_1Yt.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_17_1_i5L.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_18_2_Pk3.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_19_1_Yl8.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_20_2_QIG.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_21_1_RVd.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_22_1_Ews.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_23_1_WVm.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_24_1_C9c.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_25_1_t2C.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_26_2_vp8.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_27_1_0i4.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_28_1_ToC.root',
        # 'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS/MergedOutputFile_29_1_N0k.root'
        
        ## ok for template studies
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_1.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_10.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_11.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_12.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_13.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_14.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_15.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_16.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_17.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_18.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_19.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_2.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_20.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_21.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_22.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_23.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_24.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_25.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_26.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_27.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_28.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_29.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_3.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_30.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_31.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_32.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_33.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_34.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_35.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_36.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_37.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_38.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_39.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_4.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_40.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_41.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_42.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_43.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_44.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_45.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_46.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_47.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_48.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_49.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_5.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_50.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_51.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_6.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_7.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_8.root',
        'rfio:/castor/cern.ch/user/m/musich/ZbbAnalysis/DY_forTemplateStudies/merged2/zbbPATSkim_Merge_Merged_9.root'
        ]

    ttSrc=[
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/TTJets_TuneZ2/MergedOutputFile_1_2_wnb.root'
        ]

    zbbSrc_SHERPA_5FS=[
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_1_3_qPY.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_2_3_qTG.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_3_1_cAj.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_4_2_f3e.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_5_2_DHq.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_6_6_SiH.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_7_8_Ks2.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_8_1_sSR.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_9_2_Uxh.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_10_8_I1u.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_11_1_40S.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_12_8_eQA.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_13_2_L1c.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_14_1_8dW.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_15_7_LwJ.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_16_2_9wq.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_17_2_JXi.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_18_8_nAn.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_19_8_lNa.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_20_7_3VH.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_21_1_Q2q.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_22_2_mVR.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_23_8_dWv.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_24_2_LkU.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_25_7_5xi.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_26_8_yek.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_27_4_Jri.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_28_1_rvj.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_29_4_tFa.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_30_1_EAM.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_31_2_oiI.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_32_1_RyR.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_33_2_815.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_34_2_ZmC.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_35_2_V9P.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_36_7_f1x.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_37_7_qBb.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_38_3_JYU.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_39_1_MP5.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_40_6_Yw8.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_41_7_H4y.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_42_6_iuM.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_43_2_7Ql.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_44_1_AAx.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_45_8_vS0.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_46_7_Sou.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_47_1_LMQ.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_48_1_9zI.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_49_1_RQe.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_50_2_0Vw.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_51_2_8sS.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_52_2_91L.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_SHERPA_5FS_bEnriched/MergedOutputFile_53_1_P1X.root'
        ]

    zbbSrc_MADGRAPH_5FS=[
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_1_1_6Bb.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_2_1_0Hw.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_3_1_BB4.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_4_1_4ce.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_5_1_MkN.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_6_1_k65.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_7_6_y0u.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_8_1_Go1.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_9_1_5xV.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_10_1_SJL.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_11_1_iAO.root',
        #'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_12_4_13F.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_12_7_8Wj.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_13_1_WSG.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_14_2_8Cc.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_15_2_g2n.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_16_1_sQ1.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_17_1_ksn.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/DYJetsToLL_TuneZ2_MADGRAPH_5FS_bEnriched/MergedOutputFile_18_1_EBS.root'
        ]

    zbbSrc_aMCATNLO_5FS=[
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_10_1_YmL.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_11_1_uaE.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_12_1_pro.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_13_1_8jN.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_14_1_RKD.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_15_1_vwd.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_16_1_vK0.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_17_1_TBj.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_18_1_qZa.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_19_1_noY.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_1_1_SlM.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_20_1_pWu.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_21_1_nrV.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_22_1_vSX.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_23_1_spX.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_24_1_Lab.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_25_1_drT.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_26_1_Az1.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_27_1_ZzD.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_28_1_8pX.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_29_1_YG5.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_2_1_zmu.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_3_1_syz.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_4_1_clu.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_5_1_U0a.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_6_1_ymq.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_7_1_DlA.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_8_1_3cm.root',
        'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZbbToLL_aMCatNLO/MergedOutputFile_9_1_K5m.root'   
        ]

    zzSrc=['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZZ_TuneZ2/MergedOutputFile_1_1_JPk.root',
           'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZZ_TuneZ2/MergedOutputFile_2_1_sgc.root',
           'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZZ_TuneZ2/MergedOutputFile_3_1_wIe.root',
           'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/ZZ_TuneZ2/MergedOutputFile_4_1_rdh.root',
           ]

    wzSrc=['rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/WZ_TuneZ2/MergedOutputFile_1_1_FiG.root',
           'rfio:/castor/cern.ch/user/e/emiglior/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/WZ_TuneZ2/MergedOutputFile_2_1_G9x.root' 
           ]

# assign the  datasets
    dySrc          = dySrc_MADGRAPH_5FS
    zbbSrc         = zbbSrc_MADGRAPH_5FS
    zbbSherpaSrc   = zbbSrc_SHERPA_5FS
    zbbaMCatNLOSrc = zbbSrc_aMCATNLO_5FS

    channels  = ["MC_Zb5fToLLaMCatNLO_OOB","MC_Zb5fToLLSherpa_OOB","MC_Zb5fToLL_OOB","MC_ZcToLL_OOB","MC_ZlJets_OOB","MC_TTJets_OOB","MC_ZtautauJets_OOB","MC_ZZ_OOB","MC_ZW_OOB", 
                 "MC_Zb5fToLLaMCatNLO_PU" ,"MC_Zb5fToLLSherpa_PU" ,"MC_Zb5fToLL_PU" ,"MC_ZcToLL_PU" ,"MC_ZlJets_PU" ,"MC_TTJets_PU" ,"MC_ZtautauJets_PU", "MC_ZZ_PU", "MC_ZW_PU", 
                 "MC_Zb5fToLLaMCatNLO_TNP","MC_Zb5fToLLSherpa_TNP","MC_Zb5fToLL_TNP","MC_ZcToLL_TNP","MC_ZlJets_TNP","MC_TTJets_TNP","MC_ZtautauJets_TNP","MC_ZZ_TNP","MC_ZW_TNP",
                 "MC_Zb5fToLLaMCatNLO_All","MC_Zb5fToLLSherpa_All","MC_Zb5fToLL_All","MC_ZcToLL_All","MC_ZlJets_All","MC_TTJets_All","MC_ZtautauJets_All","MC_ZZ_All","MC_ZW_All"]

    srcFiles  =  [zbbaMCatNLOSrc,zbbSherpaSrc,zbbSrc,dySrc,dySrc,ttSrc,dySrc,zzSrc,wzSrc,
                  zbbaMCatNLOSrc,zbbSherpaSrc,zbbSrc,dySrc,dySrc,ttSrc,dySrc,zzSrc,wzSrc,
                  zbbaMCatNLOSrc,zbbSherpaSrc,zbbSrc,dySrc,dySrc,ttSrc,dySrc,zzSrc,wzSrc,
                  zbbaMCatNLOSrc,zbbSherpaSrc,zbbSrc,dySrc,dySrc,ttSrc,dySrc,zzSrc,wzSrc]

    isData =   [False,False,False,False,False,False,False,False,False,
                False,False,False,False,False,False,False,False,False,
                False,False,False,False,False,False,False,False,False,
                False,False,False,False,False,False,False,False,False]
    
    selections = ["Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample",
                  "Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample",
                  "Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample",
                  "Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample"]

    samples     = ["Zb5faMCatNLO","Zb5fSherpa","Zb5f","Zc","Zl","tt","Ztautau","zz","wz",
                   "Zb5faMCatNLO","Zb5fSherpa","Zb5f","Zc","Zl","tt","Ztautau","zz","wz",
                   "Zb5faMCatNLO","Zb5fSherpa","Zb5f","Zc","Zl","tt","Ztautau","zz","wz",
                   "Zb5faMCatNLO","Zb5fSherpa","Zb5f","Zc","Zl","tt","Ztautau","zz","wz"]

    pusels      = ["False","False","False","False","False","False","False","False","False",
                   "True" ,"True" ,"True" ,"True" ,"True" ,"True","True","True","True",
                   "True" ,"True" ,"True" ,"True" ,"True" ,"True","True","True","True",
                   "True" ,"True" ,"True" ,"True" ,"True" ,"True","True","True","True"]
                                                                                                                                                      
    tnpsels     = ["False","False","False","False","False","False","False", "False","False",
                   "False","False","False","False","False","False","False","False","False",
                   "True" ,"True" ,"True" ,"True" ,"True" ,"True","True","True","True",
                   "True" ,"True" ,"True" ,"True" ,"True" ,"True","True","True","True"]
                                                                                                                                                  
    bsels       = ["False","False","False","False","False","False","False", "False","False",
                   "False","False","False","False","False","False","False","False","False",
                   "False","False","False","False","False","False","False","False","False",
                   "True" ,"True" ,"True" ,"True" ,"True" ,"True","True","True","True"]

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

            
        for jobN,theSrcFiles in enumerate(split(srcFiles[iChannel],3)):            
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


   
