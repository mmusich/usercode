#!/usr/bin/env python

import datetime,time
import os,sys
import string, re
import subprocess
import ConfigParser
from optparse import OptionParser


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

def set_global_var(sample):
    global USER
    global HOME
    global PBS_DIR
    global LOG_DIR
    
    global SCRAM_ARCH
    global CMSSW_VER

    global GENSIM_FILE

    USER = os.environ.get('USER')
    HOME = os.environ.get('HOME')
    PBS_DIR = os.path.join(os.getcwd(),"PBS") 
    LOG_DIR = os.path.join(os.getcwd(),"log")
    SCRAM_ARCH = "slc5_amd64_gcc472"
    CMSSW_VER="CMSSW_6_2_0_SLHC10"
    
    if (sample=="TTbar") | (sample=="ttbar") | (sample=="TTBar") :
        GENSIM_FILE = "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/612_slhc8/Extended2017/TTbar/step1_TTtoAnything_14TeV_pythia6_15k_evts.root"
    elif (sample=="MinBias") | (sample=="minbias") :
        GENSIM_FILE = "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/612_slhc8/Extended2017/MinBias/step1_MinBias_TuneZ2star_14TeV_pythia6_15k_evts.root"
    elif (sample=="IsoMuons") | (sample=="muons") | (sample=="Muons") :
        GENSIM_FILE = "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/612_slhc8/Extended2017/ParticleGun/step1_FourMuPartGun_100kEvents.root"
    else :
        print "unrecongnize input sample, using default (=TTbar)"
        GENSIM_FILE = "root://eoscms//eos/cms/store/caf/user/emiglior/SLHCSimPhase2/612_slhc8/Extended2017/TTbar/step1_TTtoAnything_14TeV_pythia6_15k_evts.root"

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
    """Main class to create and submit PBS jobs"""
###########################################################################

    def __init__(self, job_id,firstevent,maxevents, sample, pu, ageing, pixelrocrows, pixelroccols, bpixthr, bpixl0thickness, pixelcpe):
############################################################################################################################
        
        # store the job-ID (since it is created in a for loop)
        self.job_id=job_id
        
        # first/max event used in this job
        self.firstevent=firstevent
        self.maxevents=maxevents

        ## FIXME: always check that these are specified
        self.sample=sample
        self.pu=pu
        
        # parameters of the pixel digitizer 
        self.pixelrocrows=pixelrocrows
        self.pixelroccols=pixelroccols
        self.bpixl0thickness=bpixl0thickness
        
        self.bpixthr=bpixthr
        self.ageing=ageing

        # parameters for RECO
        self.pixelcpe=pixelcpe

# >>>>>>>>> BA
#                self.out_dir=os.path.join("/lustre/cms/store/user",USER,"SLHCSimPhase2/out","sample_"+sample,"pu_"+pu,"PixelROCRows_"+pixelrocrows+"_PixelROCCols_"+pixelroccols,"L0Thick_"+self.bpixl0thickness, "BPixThr_"+bpixthr)
        self.out_dir=os.path.join("/store/caf/user",USER,"TMP","SLHCSimPhase2/out","sample_"+sample,"pu_"+pu,"PixelROCRows_"+pixelrocrows+"_PixelROCCols_"+pixelroccols,"L0Thick_"+self.bpixl0thickness, "BPixThr_"+bpixthr)
# <<<<<<<<< LXBATCH
#        os.system("mkdir -p "+self.out_dir)
        mkdir_eos(self.out_dir)        

        self.job_basename= "step_digitodqm_" +self.sample+ "_pu" + self.pu + "_age" + self.ageing + "_PixelROCRows" + self.pixelrocrows + "_PixelROCCols" + self.pixelroccols + "_L0Thick" + self.bpixl0thickness + "_BPixThr" + self.bpixthr + "_event" + str(self.firstevent)
        
        self.cfg_dir=None
        self.outputPSetName=None

# PBS variables        
        self.output_PBS_name=None

###############################
    def createThePBSFile(self):
###############################

# directory to store the PBS to be submitted
        self.pbs_dir = PBS_DIR
        if not os.path.exists(self.pbs_dir):
            os.makedirs(self.pbs_dir)

        jobs_dir = os.path.join(PBS_DIR,"jobs")
        if not os.path.exists(jobs_dir):
            os.makedirs(jobs_dir)

        self.output_PBS_name=self.job_basename+".pbs"
        fout=open(os.path.join(self.pbs_dir,'jobs',self.output_PBS_name),'w')
    
        if not os.path.exists(LOG_DIR):
            os.makedirs(LOG_DIR)

        fout.write("#!/bin/sh \n") 
        fout.write("#BSUB -L /bin/sh \n")       
        fout.write("#BSUB -J "+self.job_basename+"\n")
        fout.write("#BSUB -oo "+os.path.join(LOG_DIR,self.job_basename)+".log \n") # LXBATCH
        fout.write("#BSUB -q cmscaf1nd \n")                                        # LXBATCH
        fout.write("#BSUB -R \"rusage[mem=5000]\"\n")
        fout.write("### Auto-Generated Script by LoopCMSSWBuildAndRun.py ### \n")
        fout.write("JobName="+self.job_basename+" \n")
        fout.write("outfilename="+self.job_basename+".root"+" \n")
        fout.write("OUT_DIR="+self.out_dir+" \n")
        fout.write("firstevent="+str(self.firstevent)+" \n")
        fout.write("maxevents="+str(self.maxevents)+" \n")
        fout.write("pixelroccols="+self.pixelroccols+" \n")
        fout.write("pixelrocrows="+self.pixelrocrows+" \n")
        fout.write("pixelcpe="+self.pixelcpe+" \n")
        fout.write("puscenario="+self.pu+" \n")
        fout.write("ageing="+self.ageing+" \n")
        fout.write("bpixthr="+self.bpixthr+" \n")
        fout.write("inputgensimfilename="+GENSIM_FILE+" \n")

        fout.write("if [ ! \"$LSB_JOBID\" = \"\" ]; then \n")
        fout.write("echo \"I AM IN BATCH\" \n")
        fout.write("export HOME=$WORKDIR \n") # LXBATCH
        fout.write("cd \n")
        fout.write("fi \n")

        
        fout.write("export SCRAM_ARCH=slc5_amd64_gcc472 \n")
        fout.write("# Setup variables   \n")

        fout.write("cmssw_ver="+CMSSW_VER+" \n")
        fout.write("# Install and Compile CMSSW on batch node  \n")
        fout.write("scram p CMSSW $cmssw_ver  \n")
        fout.write("cd ${cmssw_ver}/src \n")
        fout.write("eval `scram r -sh` \n")

        # implement in the PBS script E.Brownson's recipe for changing the size of the pixels / part #1
        fout.write("# Eric Brownson's recipe to change the size of the pixels \n")
        fout.write("### 1: checkout CMSSW patches \n")

        fout.write("if [ ! \"$LSB_JOBID\" = \"\" ]; then \n")
        fout.write("cd \n")
        fout.write("# git config needed to avoid \n")
        fout.write("# error: SSL certificate problem: unable to get local issuer certificate while accessing \n")
        fout.write("ln -fs "+os.path.join(HOME,".gitconfig")+ " .\n")        
        fout.write("# since /cmshome/emiglior is mounted on the batch node \n")
        fout.write("eval `ssh-agent -s` \n")
        fout.write("ssh-add "+os.path.join(HOME,".ssh","id_rsa")+"\n")
        fout.write("# checkpoint: test that ssh connection is ok \n")
        fout.write("ssh -T git@github.com \n")
        fout.write("git clone https://github.com/fwyzard/cms-git-tools \n")
        fout.write("PATH=$PWD/cms-git-tools:$PATH \n")
        fout.write("git cms-init -y --ssh \n")
        fout.write("cd $HOME/${cmssw_ver}/src  \n")
        fout.write("fi \n")
        fout.write("git cms-addpkg CalibTracker/SiPixelESProducers \n")
        fout.write("git cms-addpkg CondFormats/SiPixelObjects \n")
        fout.write("git cms-addpkg CondTools/SiPixel \n")
        fout.write("git cms-addpkg Geometry/TrackerCommonData \n")
        fout.write("git cms-addpkg SLHCUpgradeSimulations/Geometry \n")
        fout.write("git cms-addpkg SLHCUpgradeSimulations/Configuration \n")
        fout.write("echo \"After git cms-addpkg\" \n")
        fout.write("pwd \n")
        fout.write("ls -l . \n")
        fout.write("git pull https://github.com/mmusich/cmssw ChangePitch_on62X \n")
        fout.write("### 1 ended  \n")
        
        fout.write("git clone -b 620_slhc10_phase1 git://github.com/emiglior/usercode.git \n")
        fout.write("mv usercode/AuxCode .\n")
        fout.write("mv usercode/RecoLocalTracker .\n")
        fout.write("rm -fr usercode \n")
        fout.write("git cms-checkdeps -a \n")

        fout.write("# compile \n")
        fout.write("scram b -j 8 \n")
        fout.write("eval `scram r -sh` \n")
        
        # implement in the PBS script E.Brownson's recipe for changing the size of the pixels / part #2
        fout.write("# Eric Brownson's recipe to change the size of the pixels \n")
        fout.write("### 2: modify the topology \n")
        fout.write("# trackerStructureTopology_template_L0.xml   -> L0    BPIX is changed \n")
        fout.write("sed -e \"s%PIXELROCROWS%"+self.pixelrocrows+"%g\" -e \"s%PIXELROCCOLS%"+self.pixelroccols+"%g\" AuxCode/SLHCSimPhase2/test/trackerStructureTopology_template_L0.xml > Geometry/TrackerCommonData/data/PhaseI/trackerStructureTopology.xml \n")
        fout.write("# Run CMSSW to complete the recipe for changing the size of the pixels \n")

        # recipe for phase I tracking  
        fout.write("cmsRun SLHCUpgradeSimulations/Geometry/test/writeFile_phase1_cfg.py \n")
        fout.write("mv PixelSkimmedGeometry_phase1.txt ${CMSSW_BASE}/src/SLHCUpgradeSimulations/Geometry/data/PhaseI \n")
        
        # recipe for phase II tracking
        # fout.write("cmsRun SLHCUpgradeSimulations/Geometry/test/writeFile_phase2BE_cfg.py \n")
        # fout.write("mv PixelSkimmedGeometry_phase2BE.txt ${CMSSW_BASE}/src/SLHCUpgradeSimulations/Geometry/data/PhaseII/BarrelEndcap/PixelSkimmedGeometry.txt \n")
        
        fout.write("### 2 ended  \n")

        # implement the recipe for changing the bpix sensor thickness from A. Tricomi
        fout.write("# A Tricomi's recipe to change the sensors thickness \n")
        fout.write("sed -e \"s%BPIXLAYER0THICKNESS%"+self.bpixl0thickness+"%g\" AuxCode/SLHCSimPhase2/test/pixbarladderfull0_template.xml > Geometry/TrackerCommonData/data/PhaseI/pixbarladderfull0.xml \n")

        fout.write("# Run CMSSW for DIGI-to-DQM steps \n")
        fout.write("cd "+os.path.join("AuxCode","SLHCSimPhase2","test")+"\n")  
        fout.write("cmsRun step_digitodqmvalidation_PUandAge.py maxEvents=${maxevents} firstEvent=${firstevent} BPixThr=${bpixthr} PixelCPE=${pixelcpe} InputFileName=${inputgensimfilename} OutFileName=${outfilename} PUScenario=${puscenario} AgeingScenario=${ageing} \n")
        fout.write("cmsStage step_digitodqmvalidation_PUandAge.py ${OUT_DIR}/step_digitodqmvalidation_PUandAge.py \n")
        fout.write("ls -lh . \n")
        fout.write(" # retrieve the outputs \n")
#        fout.write("for RootOutputFile in $(ls *root ); do cp  ${RootOutputFile}  ${OUT_DIR}/${RootOutputFile} ; done \n")
        fout.write("for RootOutputFile in $(ls *root ); do cmsStage  ${RootOutputFile}  ${OUT_DIR}/${RootOutputFile} ; done \n")
        fout.write("#cp ${CMSSW_BASE}/src/SLHCUpgradeSimulations/Geometry/data/PhaseI/PixelSkimmedGeometry_phase1.txt ${OUT_DIR} \n")
        fout.write("#cp ${CMSSW_BASE}/src/Geometry/TrackerCommonData/data/PhaseI/trackerStructureTopology.xml ${OUT_DIR} \n")
        fout.close()

############################################
    def submit(self):
############################################
        os.system("chmod u+x " + os.path.join(self.pbs_dir,'jobs',self.output_PBS_name))
        os.system("bsub < "+os.path.join(self.pbs_dir,'jobs',self.output_PBS_name)) #LXBATCH

#################
def main():            
### MAIN LOOP ###

    desc="""This is a description of %prog."""
    parser = OptionParser(description=desc,version='%prog version 0.1')
    parser.add_option('-s','--submit',          help='job submitted',           dest='submit', action='store_true', default=False)
    parser.add_option('-n','--numberofjobs',    help='number of splitted jobs', dest='numberofjobs', action='store',  default=1)
    parser.add_option('-j','--jobname',         help='task name',               dest='jobname', action='store', default='myjob')
    parser.add_option('-r','--ROCRows',         help='ROC Rows (default 80 -> du=100 um)', dest='rocrows', action='store', default='80')
    parser.add_option('-c','--ROCCols',         help='ROC Cols (default 52 -> dv=150 um)', dest='roccols', action='store', default='52')
    parser.add_option('-t','--Layer0Thick',     help='BPix L0 sensor thickness',dest='layer0thick', action='store', default='0.285')
    parser.add_option('-T','--BPixThr',         help='BPix Threshold',          dest='bpixthr', action='store', default='2000')
    parser.add_option('-p','--pileup',          help='set pileup',              dest='pu',action='store',default='NoPU')
    parser.add_option('-S','--sample',          help='set sample name',         dest='sample',action='store',default='TTbar')
    parser.add_option('-a','--ageing',          help='set ageing',              dest='ageing',action='store',default='NoAgeing')
    parser.add_option('-i','--input',           help='set input configuration (overrides default)',dest='inputconfig',action='store',default=None)
    (opts, args) = parser.parse_args()

    # initialize needed input 
    mRocRows  = None
    mRocCols  = None
    mBPixThr  = None
    mL0Thick  = None
    mAgeing   = None
    mSample   = None
    mPileUp   = None
    mPixelCPE = None

    ConfigFile = opts.inputconfig
    
    if ConfigFile is not None:

        print "********************************************************"
        print "*         Parsing from input file:", ConfigFile,"    "
        
        config = ConfigParser.ConfigParser()
        config.read(ConfigFile)
        
        mRocRows  = ConfigSectionMap(config,"PixelConfiguration")['rocrows']   
        mRocCols  = ConfigSectionMap(config,"PixelConfiguration")['roccols']   
        mBPixThr  = ConfigSectionMap(config,"PixelConfiguration")['bpixthr']   
        mL0Thick  = ConfigSectionMap(config,"PixelConfiguration")['layer0thickness']
        mAgeing   = ConfigSectionMap(config,"PixelConfiguration")['ageing']    
        mSample   = ConfigSectionMap(config,"SampleConfiguration")['sample']   
        mPileUp   = ConfigSectionMap(config,"SampleConfiguration")['pileup']
        mPixelCPE = ConfigSectionMap(config,"SampleConfiguration")['pixelcpe']   

    else :

        print "********************************************************"
        print "*             Parsing from command line                *"
        print "********************************************************"
        
        mRocRows  = opts.rocrows
        mRocCols  = opts.roccols
        mBPixThr  = opts.bpixthr
        mL0Thick  = opts.layer0thick
        mAgeing   = opts.ageing
        mSample   = opts.sample
        mPileUp   = opts.pu
        mPixelCPE = 'pixelCPE_100x150_upgrade'

# check that chosen pixel size matches what is currently available in the trackerStructureTopology
# https://twiki.cern.ch/twiki/bin/view/CMS/ExamplePhaseI#Changing_the_Pixel_Size
    if int(mRocRows) % 80:
        print 'illegal value for PixelROCRows' 
    exit

    if int(mRocCols) % 52:
        print "illegal value for PixelROCCols"
    exit

    # Set global variables
    set_global_var(mSample)
    
    print "********************************************************"
    print "*                 Configuration info                   *"
    print "********************************************************"
    print "  Launching this script from : ",os.getcwd()
    print "- submitted                  : ",opts.submit
    print "- Jobname                    : ",opts.jobname
    print "- Sample                     : ",mSample
    print "- Input generated sample     : ",GENSIM_FILE
    print "- PileUp Scenario            : ",mPileUp
    print "- ROCRows                    : ",mRocRows
    print "- ROCCols                    : ",mRocCols
    print "- L0 Thickness               : ",mL0Thick
    print "- Clusterizer Threshold      : ",mBPixThr
    print "- Ageing Scenario            : ",mAgeing
    
    # Setup CMSSW variables
#     os.system("source /swcms_slc5/CMSSW/cmsset_default.sh") # LXBATCH
#     os.chdir(os.path.join(HOME,"SLHCSimPhase2",CMSSW_VER,"src"))
    os.chdir(os.path.join(HOME,"MyWorkSpace/public","TMP","SLHCSimPhase2",CMSSW_VER,"src"))
    os.system("eval `scram r -sh`")

    # Split and submit
    child_edm = subprocess.Popen(["edmEventSize","-v",GENSIM_FILE],stdout=subprocess.PIPE)
    (out,err) = child_edm.communicate()

    ### uncomment next to debug the script on 50 events
    nEvents=10 # this line should be commented for running on the full GEN-SIM sample
    #nEvents = int((out.split("\n")[1]).split()[3])
    
    eventsPerJob = nEvents/int(opts.numberofjobs)

    print "********************************************************"
    print "*                 Job submission pattern               *"
    print "********************************************************"    
    print "- Total events to run       : ",nEvents," in ",opts.numberofjobs," jobs" 
    print "- Total events/job          : ",eventsPerJob
    
    remainder=nEvents
    firstEvent=1
    jobIndex=0
    
    #prepare the list of the DQM files for the harvesting
    DQMFileList=""
    
    out_dir = None

    #############################################
    # loop on the jobs
    #############################################
    
    while remainder>0:

        print "- Job n. ",jobIndex," will process events from: ",firstEvent," to ",firstEvent+eventsPerJob-1
        
        ajob=Job(opts.jobname, firstEvent, eventsPerJob, mSample, mPileUp, mAgeing, mRocRows, mRocCols, mBPixThr, mL0Thick, mPixelCPE)
        ajob.createThePBSFile()

        dqmoutput=ajob.job_basename+".root"
        dqmoutput.replace("digitodqm_","digitodqm_inDQM")        
        # this is needed for the script doing the harvesting
        DQMFileList+="file:"+os.path.join(ajob.out_dir,dqmoutput)+","

        out_dir = ajob.out_dir # save for later usage
        
        if opts.submit:
            ajob.submit()
            del ajob

        remainder -= eventsPerJob
        firstEvent += eventsPerJob
        jobIndex+=1
        
    #############################################
    # link the output folder
    #############################################
    
    link_name="sample_"+mSample+"_pu"+mPileUp+"_PixelROCRows_"+mRocRows+"_PixelROCCols_"+mRocCols+"_BPixThr_"+mBPixThr+"_L0Thick"+mL0Thick
    linkthedir="ln -fs "+out_dir+" "+os.path.join(LOG_DIR,link_name)     
    os.system(linkthedir)    

    print "- Output will be saved in   :",out_dir
    print "********************************************************"

    #############################################
    # prepare the script for the harvesting step
    #############################################
    
    harvestingname = PBS_DIR + "/jobs/"+opts.jobname+"_sample_"+mSample+"_pu"+mPileUp+"_PixelRocRows"+mRocRows+"_PixelROCCols_"+mRocCols+"_BPixThr"+mBPixThr+"_L0Thick"+mL0Thick+".sh"
    fout=open(harvestingname,"w")
    fout.write("#!/bin/sh \n")
#     fout.write("source /swcms_slc5/CMSSW/cmsset_default.sh \n") # LXBATCH
    fout.write("cmssw_ver="+CMSSW_VER+" \n")
#    fout.write("cd "+os.path.join(HOME,"SLHCSimPhase2","${cmssw_ver}","src")+"\n") # LXBTACH
    fout.write("cd "+os.path.join(HOME,"MyWorkSpace/public","TMP","SLHCSimPhase2","${cmssw_ver}","src")+"\n")
    fout.write("eval `scram r -sh`\n")
    fout.write("DQMFileList="+DQMFileList[:-1]+" \n")
    fout.write("cmsDriver.py step4 --geometry Extended2017 --customise SLHCUpgradeSimulations/Configuration/phase1TkCustoms.customise,AuxCode/SLHCSimPhase2/TkOnlyValidationCustoms.customise_tkonly --conditions auto:upgrade2017 --mc -s HARVESTING:validationHarvesting+dqmHarvesting --filein $DQMFileList --fileout file:step4_sample_"+mSample+"_pu"+mPileUp+"_PixelRocRows"+mRocRows+"_PixelROCCols_"+mRocCols+"_BPixThr"+mBPixThr+"_L0Thick"+mL0Thick+".root > step4_sample_"+mSample+"_pu"+mPileUp+"_PixelRocRows"+mRocRows+"_PixelROCCols_"+mRocCols+"_BPixThr"+mBPixThr+"_L0Thick"+mL0Thick+".log \n")
    fout.write("mv DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root AuxCode/SLHCSimPhase2/test/step4_sample_"+mSample+"_pu"+mPileUp+"_PixelRocRows"+mRocRows+"_PixelROCCols_"+mRocCols+"_BPixThr"+ mBPixThr+"_L0Thick"+mL0Thick+".root")

 
    fout.close()

if __name__ == "__main__":        
    main()

    


