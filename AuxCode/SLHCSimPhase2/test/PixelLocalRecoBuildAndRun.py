#!/usr/bin/env python

import datetime,time
import os,sys
import string, re
import subprocess
import ConfigParser
from optparse import OptionParser,OptionGroup

##### method to parse the input file ################################

def ConfigSectionMap(config, section):
    the_dict = {}
    options = config.options(section)
    for option in options:
        try:
            the_dict[option] = config.get(section, option)
            if the_dict[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            the_dict[option] = None
    return the_dict

###### method to define global variable instead of export #############

def set_global_var():
    global USER
    global HOME
    global LSF_DIR
    global LOG_DIR
    
    global SCRAM_ARCH
    global CMSSW_VER
    global LAUNCH_BASE

    USER = os.environ.get('USER')
    HOME = os.environ.get('HOME')
    LAUNCH_BASE = os.environ.get('CMSSW_BASE')
    LSF_DIR = os.path.join(os.getcwd(),"LSF")
    LOG_DIR = os.path.join(os.getcwd(),"log")
    SCRAM_ARCH = "slc6_amd64_gcc472"
    CMSSW_VER="CMSSW_6_2_0_SLHC17_patch1"
    
###### method to create recursively directories on EOS  #############
    
def mkdir_eos(out_path):
    newpath='/'
    for dir in out_path.split('/'):
        newpath=os.path.join(newpath,dir)
        # do not issue mkdir from very top of the tree
        if newpath.find('SLHCSimPhase2') > 0:
            p = subprocess.Popen(["cmsMkdir",newpath], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out, err) = p.communicate()
            p.wait()

    # now check that the directory exists
    p = subprocess.Popen(["cmsLs",out_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = p.communicate()
    p.wait()
    if p.returncode !=0:
        print out

###########################################################################
class Job:
    """Main class to create and submit LSF jobs"""
###########################################################################

    def __init__(self, job_id, maxevents, pu, ageing, pixelrocrows_l0l1, pixelroccols_l0l1, pixelrocrows_l2l3, pixelroccols_l2l3, pixeleperadc, pixmaxadc, chanthr, seedthr, clusthr, bpixthr, bpixl0l1thickness, bpixl2l3thickness, myseed, islocal, queue, task_name):
############################################################################################################################
        
        # store the job-ID (since it is created in a for loop)
        self.job_id=job_id
        self.task_name=task_name
        self.queue=queue
        
        # max event used in this job
        self.maxevents=maxevents

        ## FIXME: always check that these are specified
        
        # parameters of the pixel digitizer 
        self.pixelrocrows_l0l1=pixelrocrows_l0l1
        self.pixelrocrows_l2l3=pixelrocrows_l2l3
        self.pixelroccols_l0l1=pixelroccols_l0l1
        self.pixelroccols_l2l3=pixelroccols_l2l3

        self.pixeleperadc=pixeleperadc 
        self.pixmaxadc=pixmaxadc

        # parameter of the pixel clusterizer
        self.chanthr=chanthr
        self.seedthr=seedthr
        self.clusthr=clusthr
        
        self.bpixl0l1thickness=bpixl0l1thickness
        self.bpixl2l3thickness=bpixl2l3thickness

        self.bpixthr=bpixthr
        self.ageing=ageing
        self.pu=pu
        self.myseed=myseed
        self.islocal=islocal
        self.launch_dir=LAUNCH_BASE

        #        self.out_dir=os.path.join("/store/caf/user",USER,"SLHCSimPhase2/phase2/out62XSLHC17patch1/32bit/TestBricked",\
        #        self.out_dir=os.path.join("/store/caf/user",USER,"SLHCSimPhase2/phase2/out62XSLHC17patch1nn/DataCompression",\
        self.out_dir=os.path.join(self.launch_dir,"src/results/SLHCSimPhase2/phase2/out62XSLHC17patch1/32bit/OccupancyStudy",\
                                      "PixelROCRows_"+pixelrocrows_l0l1+"_"+pixelrocrows_l2l3,\
                                      "PixelROCCols_"+pixelroccols_l0l1+"_"+pixelroccols_l2l3,\
                                      "BPIXThick_"+bpixl0l1thickness+"_"+bpixl2l3thickness,\
                                      "BPixThr_"+bpixthr,"eToADC_"+pixeleperadc,"MaxADC_"+pixmaxadc,\
                                      "ChanThr_"+chanthr,"SeedThr_"+seedthr,"ClusThr_"+clusthr,\
                                      "PU_"+self.pu)
    
        if(self.job_id==1):
            os.makedirs(self.out_dir)
           
            
        self.task_basename = self.task_name + "_pixelCPE_age" + self.ageing + "_pu" + self.pu + "_PixelROCRows" + self.pixelrocrows_l0l1 + "_" +  self.pixelrocrows_l2l3 +\
            "_PixelROCCols" + self.pixelroccols_l0l1 + "_" +  self.pixelroccols_l2l3 +\
            "_BPIXThick" + self.bpixl0l1thickness + "_" + self.bpixl2l3thickness +\
            "_BPixThr" + self.bpixthr + 'eToADC_' + self.pixeleperadc + 'MaxADC_' + self.pixmaxadc + 'ChanThr_'+ self.chanthr + 'SeedThr_' + self.seedthr + "ClusThr_" + self.clusthr 

        self.job_basename  = self.task_basename + "_seed" +str(self.myseed)

        self.cfg_dir=None
        self.outputPSetName=None

        # LSF variables        
        self.output_LSF_name=None

###############################
    def createTheLSFFile(self):
###############################

        # directory to store the LSF to be submitted
        self.lsf_dir = LSF_DIR
        if not os.path.exists(self.lsf_dir):
            os.makedirs(self.lsf_dir)

        jobs_dir = os.path.join(LSF_DIR,"jobs")
        if not os.path.exists(jobs_dir):
            os.makedirs(jobs_dir)

        self.output_LSF_name=self.job_basename+".lsf"
        fout=open(os.path.join(self.lsf_dir,'jobs',self.output_LSF_name),'w')
    
        if not os.path.exists(LOG_DIR):
            os.makedirs(LOG_DIR)

        fout.write("#!/bin/sh \n") 

        fout.write("#BSUB -L /bin/sh \n")       
        fout.write("#BSUB -J "+self.job_basename+"\n")

        fout.write("#BSUB -oo "+os.path.join(LOG_DIR,self.job_basename)+".log \n") # LXBATCH     
        fout.write("#BSUB -q "+self.queue+" \n")                                   # LXBATCH
     
        fout.write("### Auto-Generated Script by PixelCPEBuildAndRun.py ### \n")
        fout.write("JobName="+self.job_basename+" \n")
        fout.write("OUT_DIR="+self.out_dir+" \n")
        fout.write("islocal="+str(self.islocal)+" \n")
        fout.write("maxevents="+str(self.maxevents)+" \n")
        fout.write("pixelroccols_l0l1="+self.pixelroccols_l0l1+" \n")
        fout.write("pixelrocrows_l0l1="+self.pixelrocrows_l0l1+" \n")
        fout.write("pixelroccols_l2l3="+self.pixelroccols_l2l3+" \n")
        fout.write("pixelrocrows_l2l3="+self.pixelrocrows_l2l3+" \n")
        fout.write("ageing="+self.ageing+" \n")
        fout.write("puscenario="+self.pu+" \n")
        fout.write("bpixthr="+self.bpixthr+" \n")
        fout.write("pixeleperadc="+self.pixeleperadc+" \n") 
        fout.write("pixmaxadc="+self.pixmaxadc+" \n")
        fout.write("chanthr="+self.chanthr+" \n")
        fout.write("seedthr="+self.seedthr+" \n")
        fout.write("clusthr="+self.clusthr+" \n")
        fout.write("myseed="+str(self.myseed)+" \n")
        
        fout.write("if [ ! \"$LSB_JOBID\" = \"\" ]; then \n")
        fout.write("echo \"I AM IN BATCH\" \n")

        fout.write("export HOME=$WORKDIR \n") # LXBATCH
        fout.write("cd \n")
        fout.write("fi \n")
                   
        fout.write("export SCRAM_ARCH=slc6_amd64_gcc472 \n")
        fout.write("# Setup variables   \n")
                                      
        fout.write("cmssw_ver="+CMSSW_VER+" \n")
        fout.write("# Install and Compile CMSSW on batch node  \n")
        fout.write("scram p CMSSW $cmssw_ver  \n")
        fout.write("cd ${cmssw_ver}/src \n")
        fout.write("eval `scram r -sh` \n")

        if(self.islocal):
            fout.write("echo \"I AM IN LOCAL MODE\" \n")
            fout.write("export PKG_DIR="+self.launch_dir+"/src/AuxCode/SLHCSimPhase2/test \n")
        else:
            fout.write("echo \"I AM NOT IN LOCAL MODE\" \n")
            fout.write("export PKG_DIR=${CMSSW_BASE}/src/AuxCode/SLHCSimPhase2/test \n")

        # implement in the LSF script E.Brownson's recipe for changing the size of the pixels / part #1
        fout.write("# Eric Brownson's recipe to change the size of the pixels \n")
        fout.write("### 1: checkout CMSSW patches \n")

        fout.write("if [ ! \"$LSB_JOBID\" = \"\" ]; then \n")
        fout.write("echo \"!!!!!!!!!!!!!!RUNNING IN BATCH!!!!!!!!!!!!! \" \n")
        fout.write("cd \n")
        fout.write("# git config needed to avoid \n")
        fout.write("# error: SSL certificate problem: unable to get local issuer certificate while accessing \n")
        fout.write("ln -fs "+os.path.join(HOME,".gitconfig")+ " .\n")        
        fout.write("# since /cmshome/emiglior is mounted on the batch node \n")
        fout.write("eval `ssh-agent -s` \n")
        fout.write("ssh-add "+os.path.join(HOME,".ssh","id_rsa")+"\n")
        fout.write("# checkpoint: test that ssh connection is ok \n")
        fout.write("ssh -T git@github.com \n")
        fout.write("cd $HOME/${cmssw_ver}/src  \n")
        fout.write("fi \n")
        fout.write("git cms-addpkg CalibTracker/SiPixelESProducers \n")
        fout.write("git cms-addpkg CondFormats/SiPixelObjects \n")
        fout.write("git cms-addpkg CondTools/SiPixel \n")
        fout.write("git cms-addpkg Geometry/TrackerCommonData \n")
        fout.write("git cms-addpkg SLHCUpgradeSimulations/Geometry \n")
        fout.write("git cms-addpkg SLHCUpgradeSimulations/Configuration \n")
        fout.write("git cms-addpkg RecoLocalTracker/SiPixelRecHits \n")
        fout.write("git cms-addpkg RecoLocalTracker/SiPixelClusterizer \n")
        fout.write("git cms-addpkg DPGAnalysis/SiStripTools \n")
        fout.write("git cms-addpkg DataFormats/SiPixelCluster \n")
        fout.write("echo \"After git cms-addpkg\" \n")
        fout.write("pwd \n")
        fout.write("ls -l . \n")
        fout.write("git pull https://github.com/mmusich/cmssw ChangePitch_on620_SLHC17_patch1_32bit \n")
        fout.write("### 1 ended  \n")

        fout.write("git clone -b 620_slhc17_patch1_phase2 git://github.com/emiglior/usercode.git \n")
        fout.write("mv usercode/AuxCode .\n")

        ###### please make sure to delete this line afterwards!!!!!! #######
        #fout.write("mv usercode/RecoLocalTracker .\n")        
        # for the moment we ignore this (to be used to change the matching window)
        #fout.write("mv usercode/SimTracker .\n")
        fout.write("rm -fr usercode \n")
        fout.write("git cms-checkdeps -a \n")

        fout.write("# compile \n")
        if(self.islocal):
            fout.write("cp "+self.launch_dir+"/src/AuxCode/SLHCSimPhase2/plugins/StdPixelHitNtuplizer.cc ./AuxCode/SLHCSimPhase2/plugins/StdPixelHitNtuplizer.cc \n")
 
        fout.write("scram b -j 8 USER_CXXFLAGS=\"-DPHASE2\" \n")
        fout.write("eval `scram r -sh` \n")

        # implement in the LSF script E.Brownson's recipe for changing the size of the pixels / part #2
        fout.write("# Eric Brownson's recipe to change the size of the pixels \n")
        fout.write("### 2: modify the topology \n")
        fout.write("# trackerStructureTopology_template.xml   -> BPIX is changed \n")
        fout.write("sed -e \"s%PIXELROCROWS_L0L1%"+self.pixelrocrows_l0l1+"%g\" -e \"s%PIXELROCCOLS_L0L1%"+self.pixelroccols_l0l1+"%g\" -e \"s%PIXELROCROWS_L2L3%"+self.pixelrocrows_l2l3+"%g\" -e \"s%PIXELROCCOLS_L2L3%"+self.pixelroccols_l2l3+"%g\" ${PKG_DIR}/trackerStructureTopology_template.xml > Geometry/TrackerCommonData/data/PhaseII/Pixel10D/trackerStructureTopology.xml \n")
        fout.write("# Run CMSSW to complete the recipe for changing the size of the pixels \n")

        # recipe for phase II tracking
        fout.write("cmsRun AuxCode/SLHCSimPhase2/test/writeFile_phase2Pixel10D_cfg.py \n")
        fout.write("mv PixelSkimmedGeometry_phase2Pixel10D.txt ${CMSSW_BASE}/src/SLHCUpgradeSimulations/Geometry/data/PhaseII/Pixel10D/PixelSkimmedGeometry.txt \n")


        fout.write("### 2 ended  \n")

        # implement the recipe for changing the bpix sensor thickness from A. Tricomi
        fout.write("# A Tricomi's recipe to change the sensors thickness \n")
        fout.write("sed -e \"s%BPIXLAYER01THICKNESS%"+self.bpixl0l1thickness+"%g\" ${PKG_DIR}/pixbarladderfull0_template.xml > Geometry/TrackerCommonData/data/PhaseI/pixbarladderfull0.xml \n")
        fout.write("sed -e \"s%BPIXLAYER01THICKNESS%"+self.bpixl0l1thickness+"%g\" ${PKG_DIR}/pixbarladderfull1_template.xml > Geometry/TrackerCommonData/data/PhaseI/pixbarladderfull1.xml \n")
        fout.write("sed -e \"s%BPIXLAYER23THICKNESS%"+self.bpixl2l3thickness+"%g\" ${PKG_DIR}/pixbarladderfull2_template.xml > Geometry/TrackerCommonData/data/PhaseI/pixbarladderfull2.xml \n")
        fout.write("sed -e \"s%BPIXLAYER23THICKNESS%"+self.bpixl2l3thickness+"%g\" ${PKG_DIR}/pixbarladderfull3_template.xml > Geometry/TrackerCommonData/data/PhaseI/pixbarladderfull3.xml \n")
        
        fout.write("# Run CMSSW for GEN-NTUPLE steps \n")
        fout.write("cd "+os.path.join("AuxCode","SLHCSimPhase2","test")+"\n")
#        fout.write("edmConfigDump ${PKG_DIR}/OneNuM_GEN_TO_RECO_PixelLocalRecoStudy_cfg.py >> pset_dumped.py \n")
#        fout.write("cmsRun ${PKG_DIR}/OneNuM_GEN_TO_RECO_PixelLocalRecoStudy_cfg.py  maxEvents=${maxevents} PixElePerADC=${pixeleperadc} PixMaxADC=${pixmaxadc} BPixThr=${bpixthr} PUScenario=${puscenario} AgeingScenario=${ageing} MySeed=${myseed} ChannelThreshold=${chanthr} SeedThreshold=${seedthr} ClusterThreshold=${clusthr} \n")
        fout.write("cmsRun ${PKG_DIR}/TenMuE_0_200_GEN_TO_RECO_PixelLocalRecoStudy_cfg.py maxEvents=${maxevents} PixElePerADC=${pixeleperadc} PixMaxADC=${pixmaxadc} BPixThr=${bpixthr} PUScenario=${puscenario} AgeingScenario=${ageing} MySeed=${myseed} ChannelThreshold=${chanthr} SeedThreshold=${seedthr} ClusterThreshold=${clusthr} \n")

        fout.write("ls -lh . \n")
#        fout.write("cmsStage -f pset_dumped.py ${OUT_DIR}/pset_dumped.py \n")
        #fout.write("cd Brownson \n")
        #fout.write("make \n")
        #fout.write("ln -fs ../stdgrechitfullph1g_ntuple.root . \n")
        #fout.write("./res \n")        
        fout.write(" # retrieve the outputs \n")
        fout.write("for RootOutputFile in $(ls *root |grep ntuple); do cmsStage -f  ${RootOutputFile}  ${OUT_DIR}/${RootOutputFile} ; done \n")
        fout.write("for RootOutputFile in $(ls *root); \n")
        fout.write("do \n")
        fout.write("events=`edmEventSize -v $RootOutputFile | grep Events | awk '{print $4}'` \n")
        fout.write("echo $RootOutputFile $events >> "+os.path.join(LAUNCH_BASE,"src","AuxCode","SLHCSimPhase2","test","eventsCount_"+ self.task_basename +".txt")+" \n")
        fout.write("done \n")
        #fout.write("for EpsOutputFile in $(ls *eps ); do cmsStage -f ${EpsOutputFile}  ${OUT_DIR}/${EpsOutputFile} ; done \n")
        fout.close() 

############################################
    def submit(self): 
############################################
        os.system("chmod u+x " + os.path.join(self.lsf_dir,'jobs',self.output_LSF_name))
        os.system("bsub < "+os.path.join(self.lsf_dir,'jobs',self.output_LSF_name)) #LXBATCH

#################
def main():            
### MAIN LOOP ###

    desc="""This is a description of %prog."""
    parser = OptionParser(description=desc,version='%prog version 0.1') #,epilog="**** Use option -i for overring all parameters from config ****")
    
    group = OptionGroup(parser,"Job configuration options",
                       "You can specify several parameters of the task"
                       )
    
    group.add_option('-j','--jobname', help='task name', dest='jobname', action='store', default='myjob')
    group.add_option('-s','--submit',  help='job submitted', dest='submit', action='store_true', default=False)
    group.add_option('-q','--queue',help='lxbatch queue for submission', dest='queue',action='store',default='cmscaf1nd')
    group.add_option('-l','--local', help='reads local branch',dest='localmode',action='store_true', default=False)
    group.add_option('-n','--numberofevents', help='number of events', dest='numberofevents', action='store', default='1')
    group.add_option('-N','--jobsInTask', help='number of jobs in this task', dest='jobsInTask', action='store',default='500')
    parser.add_option_group(group)

    group = OptionGroup(parser, "Pixel layout options",
                        "You can specift several parameters of the pixel layout")
    
    group.add_option('--ROCRows01',help='ROC Rows L0L1 (default 80 -> du=100 um)', dest='rocrows_l0l1', action='store', default='80')
    group.add_option('--ROCCols01',help='ROC Cols L0L1 (default 52 -> dv=150 um)', dest='roccols_l0l1', action='store', default='52')
    group.add_option('--Layer01Thick',help='BPix L0L1 sensor thickness', dest='layer01thick', action='store', default='0.285')

    group.add_option('--ROCRows23',help='ROC Rows L2L3 (default 80 -> du=100 um)', dest='rocrows_l2l3', action='store', default='80')
    group.add_option('--ROCCols23',help='ROC Cols L2L3(default 52 -> dv=150 um)', dest='roccols_l2l3', action='store', default='52')
    group.add_option('--Layer23Thick',help='BPix L2L3 sensor thickness', dest='layer23thick', action='store', default='0.285')

    group.add_option('-T','--BPixThr',help='BPix Threshold', dest='bpixthr', action='store', default='2000')

    group.add_option('--PixElePerADC',help='Pix ele per ADC', dest='pixeleperadc', action='store', default='135')
    group.add_option('--PixMaxADC'   ,help='Pix max ADC',     dest='pixmaxadc',  action='store', default='255')
        
    group.add_option('--SeedThr'   ,help='Cluster seed threshold',    dest='seedthr',  action='store', default='1000')
    group.add_option('--ChanThr'   ,help='Cluster channel threshold', dest='chanthr',  action='store', default='1000')
    group.add_option('--ClusThr'   ,help='Cluster channel threshold', dest='clusthr',  action='store', default='1000')

    group.add_option('-p','--pileup',help='set pileup',dest='pu',    action='store',default='NoPU')
    group.add_option('-a','--ageing',help='set ageing',dest='ageing',action='store',default='NoAgeing')
    parser.add_option_group(group)

    parser.add_option('-i','--input',help='set input configuration (overrides default)',dest='inputconfig',action='store',default=None)
    (opts, args) = parser.parse_args()

    # initialize needed input 
    mRocRows_l0l1 = None
    mRocCols_l0l1 = None
    mL0L1Thick = None

    mRocRows_l2l3 = None
    mRocCols_l2l3 = None
    mL2L3Thick = None

    mBPixThr = None
    mPixElePerADC = None
    mPixMaxADC    = None
    mSeedThr = None
    mClusThr = None
    mChanThr = None
    mAgeing  = None
    mPileUp  = None
    mQueue   = None
    mJobsInTask=None
 
    ConfigFile = opts.inputconfig
    
    if ConfigFile is not None:

        print "********************************************************"
        print "*         Parsing from input file:", ConfigFile,"    "
        
        config = ConfigParser.ConfigParser()
        config.read(ConfigFile)
        
        mRocRows_l0l1    = ConfigSectionMap(config,"PixelConfiguration")['rocrows_l0l1']   
        mRocCols_l0l1    = ConfigSectionMap(config,"PixelConfiguration")['roccols_l0l1']   
        mL0L1Thick       = ConfigSectionMap(config,"PixelConfiguration")['layer01thickness']
        mRocRows_l2l3    = ConfigSectionMap(config,"PixelConfiguration")['rocrows_l2l3']   
        mRocCols_l2l3    = ConfigSectionMap(config,"PixelConfiguration")['roccols_l2l3']   
        mL2L3Thick       = ConfigSectionMap(config,"PixelConfiguration")['layer23thickness']
        mBPixThr    = ConfigSectionMap(config,"PixelConfiguration")['bpixthr']
        mPixElePerADC = ConfigSectionMap(config,"PixelConfiguration")['pixeleperadc']
        mPixMaxADC    = ConfigSectionMap(config,"PixelConfiguration")['pixmaxadc']   
        mSeedThr    = ConfigSectionMap(config,"PixelConfiguration")['seedthr']
        mClusThr    = ConfigSectionMap(config,"PixelConfiguration")['clusthr']   
        mChanThr    = ConfigSectionMap(config,"PixelConfiguration")['chanthr']               
        mAgeing     = ConfigSectionMap(config,"PixelConfiguration")['ageing']    
        mPileUp     = ConfigSectionMap(config,"PixelConfiguration")['pileup']
        mNOfEvents  = ConfigSectionMap(config,"JobConfiguration")['numberofevents']
        mJobsInTask = ConfigSectionMap(config,"JobConfiguration")['numberofjobs']
        mQueue      = ConfigSectionMap(config,"JobConfiguration")['queue']
    
    else :

        print "********************************************************"
        print "*             Parsing from command line                *"
        print "********************************************************"
        
        mRocRows_l0l1    = opts.rocrows_l0l1
        mRocCols_l0l1    = opts.roccols_l0l1
        mL0L1Thick    = opts.layer01thick
        mRocRows_l2l3    = opts.rocrows_l2l3
        mRocCols_l2l3    = opts.roccols_l2l3
        mL2L3Thick    = opts.layer23thick
        mBPixThr    = opts.bpixthr
        mPixElePerADC = opts.pixeleperadc
        mPixMaxADC    = opts.pixmaxadc
        mSeedThr    = opts.seedthr
        mClusThr    = opts.clusthr
        mChanThr    = opts.chanthr
        mAgeing     = opts.ageing
        mPileUp     = opts.pu
        mNOfEvents  = opts.numberofevents
        mJobsInTask = opts.jobsInTask
        mQueue      = opts.queue
     
    # # check that chosen pixel size matches what is currently available in the trackerStructureTopology
    # # https://twiki.cern.ch/twiki/bin/view/CMS/ExamplePhaseI#Changing_the_Pixel_Size
    # if int(mRocRows) % 80:
    #     print 'illegal value for PixelROCRows' 
    # exit

    # if int(mRocCols) % 52:
    #     print "illegal value for PixelROCCols"
    # exit

    # Set global variables
    set_global_var()
 
    print "********************************************************"
    print "*                 Configuration info                   *"
    print "********************************************************"
    print "  Launching this script from : ",os.getcwd()
    print "- submitted                  : ",opts.submit
    print "- isLocal version            : ",opts.localmode
    print "- Jobs in Task               : ",mJobsInTask
    print "- Events/Job                 : ",mNOfEvents
    print "- Queue                      : ",mQueue
    print "- Jobname                    : ",opts.jobname
    print "- e per ADC                  : ",mPixElePerADC
    print "- Max ADC                    : ",mPixMaxADC
    print "- Clusterizer Channel Thresh : ",mChanThr
    print "- Clusterizer Seed Thresh    : ",mSeedThr
    print "- Clusterizer Cluster Thresh : ",mClusThr
    print "- ROCRows                    : ",mRocRows_l0l1,mRocRows_l2l3
    print "- ROCCols                    : ",mRocCols_l0l1,mRocCols_l2l3
    print "- L0 Thickness               : ",mL0L1Thick,mL2L3Thick
    print "- BPIX Digitizer Threshold   : ",mBPixThr
    print "- Ageing Scenario            : ",mAgeing
    print "- PileUp Scenario            : ",mPileUp
    print "- Total events to run        : ",int(mNOfEvents)*int(mJobsInTask)     

    nEvents=int(mNOfEvents)
    
    jobIndex=0

    for theseed in range(1,int(mJobsInTask)+1):

        ajob=Job(theseed, nEvents, mPileUp, mAgeing, mRocRows_l0l1, mRocCols_l0l1, mRocRows_l2l3, mRocCols_l2l3, mPixElePerADC, mPixMaxADC, mChanThr, mSeedThr, mClusThr, mBPixThr, mL0L1Thick, mL2L3Thick, theseed, opts.localmode,mQueue,opts.jobname)
        ajob.createTheLSFFile()        

        out_dir = ajob.out_dir # save for later usage
    
        if opts.submit:
            ajob.submit()
            del ajob

        jobIndex+=1       
            
        print "- Output will be saved in   :",out_dir
    print "********************************************************"

    #############################################
    # prepare the script for the harvesting step
    #############################################

    link_name=opts.jobname+"_PixelROCRows_"+mRocRows_l0l1+"_"+mRocRows_l2l3\
                          +"_PixelROCCols_"+mRocCols_l0l1+"_"+mRocCols_l2l3\
                          +"_Thick_"+mL0L1Thick+"_"+mL2L3Thick+"_BPixThr_"+mBPixThr

    harvestingname = LSF_DIR + "/jobs/PixelLocalRecoNtuple_"+link_name +".sh"
    fout=open(harvestingname,"w")

    fout.write("#!/bin/bash \n")
    fout.write("OUT_DIR="+out_dir+" \n")
    fout.write("cd "+os.path.join(LAUNCH_BASE,"src","AuxCode","SLHCSimPhase2","test")+"\n")
    fout.write("eval `scram r -sh` \n")
    fout.write("mkdir -p /tmp/$USER/"+link_name+" \n")
    fout.write("for inputfile in `cmsLs "+out_dir+" |grep seed | grep root`; do \n")
    fout.write("   namebase=`echo $inputfile |awk '{split($0,b,\"/\"); print b[21]}'` \n")
    fout.write("   cmsStage -f $OUT_DIR/$namebase /tmp/$USER/"+link_name+" \n")
    fout.write("# Uncomment next line to clean up EOS space \n")
    fout.write("#  cmsRm $OUT_DIR/$namebase \n")
    fout.write("done \n")
    fout.write("cd /tmp/$USER/"+link_name+" \n")
    fout.write("hadd /tmp/$USER/pixellocalreco_ntuple_"+link_name+".root /tmp/$USER/"+link_name+"/*_ntuple_*.root \n")
    fout.write("cmsStage -f /tmp/$USER/pixellocalreco_ntuple_"+link_name+".root $OUT_DIR \n")
    os.system("chmod u+x "+harvestingname)
    
if __name__ == "__main__":        
    main()

    


