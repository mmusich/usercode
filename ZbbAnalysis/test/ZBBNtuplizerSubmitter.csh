#!/bin/tcsh

echo  "Job started at "
\date;

if ($#argv < 1) then
    echo "************** Argument Error: at least 1 arg. required **************"
    echo "*   Usage:                                                           *"
    echo "*     ./ZBBNtuplizerSubmitter.csh <TagFile.dat> <options>            *"
    echo "*                                                                    *"
    echo "*   Options:                                                         *"
    echo "*   --dryRun             do not submit the job to lxbatch            *"
    echo "**********************************************************************"
    exit 1
endif

set inputsource=$1
set options=$2

setenv CMSSW_VER 4_1_2
setenv CMSSW_DIR ${CMSSW_BASE}/src/ZZAnalysis/Examples/test 
setenv LXBATCH_DIR `pwd`

cd $CMSSW_DIR
eval `scramv1 runtime -csh`
cd $LXBATCH_DIR

cp ${CMSSW_DIR}/${inputsource} .

set jobname=`tail ${inputsource} | grep jobname | awk '{print $2}'`
set globaltag=`tail ${inputsource} | grep globaltag | awk '{print $2}'`
set infile=`tail ${inputsource} | grep infile | awk '{print $2}'`
set plotoutfile=`tail ${inputsource} | grep plotoutfile | awk '{print $2}'`
set outfile=`tail ${inputsource} | grep outfile | awk '{print $2}'`

cp ${CMSSW_DIR}/ZBBNtuplizer_TEMPL_cfg.py .
cat ZBBNtuplizer_TEMPL_cfg.py | sed "s?GLOBALTAGTEMPLATE?${globaltag}?g" | sed "s?INPUTFILETEMPLATE?${infile}?g" |sed "s?PLOTOUTFILETEMPLATE?${plotoutfile}?g" |sed "s?OUTFILETEMPLATE?${outfile}?g" >! ${jobname}_cfg.py

cp ${CMSSW_DIR}/ZBBNtuplizer_TEMPL.lsf .
cat ZBBNtuplizer_TEMPL.lsf | sed  "s?JOBNAMETEMPLATE?${jobname}?g" | sed "s?OUTFILETEMPLATE?${outfile}?g" >! ${jobname}.lsf

if (! -d  ${CMSSW_DIR}/submittedCfg) then
    mkdir ${CMSSW_DIR}/submittedCfg
    cp ${jobname}_cfg.py  ${CMSSW_DIR}/submittedCfg
    cp ${jobname}.lsf     ${CMSSW_DIR}/submittedCfg
else
    echo "${CMSSW_DIR}/submittedCfg already exists"
    cp ${jobname}_cfg.py  ${CMSSW_DIR}/submittedCfg
    cp ${jobname}.lsf     ${CMSSW_DIR}/submittedCfg
endif

if(${options} != "--dryRun") then
 echo "cmsRun ${jobname}_cfg.py"
 cmsRun ${jobname}_cfg.py >& ${jobname}.out;

 echo "Content of working directory is: "
 \ls -lrt

 if (! -d  ${CMSSW_DIR}/NtuplizerResults) then
     mkdir ${CMSSW_DIR}/NtuplizerResults
     rfcp ${outfile} $CASTOR_HOME/Zbb2010Ntuples/v3
     cp ${plotoutfile} ${CMSSW_DIR}/NtuplizerResults 
     cp ${jobname}.out ${CMSSW_DIR}/NtuplizerResults 
 else
     rfcp ${outfile}   
     rfcp ${outfile} $CASTOR_HOME/Zbb2010Ntuples/v3
     cp ${plotoutfile} ${CMSSW_DIR}/NtuplizerResults 
     cp ${jobname}.out ${CMSSW_DIR}/NtuplizerResults 
 endif
endif

echo  "Job ended at "
\date;

exit 0
