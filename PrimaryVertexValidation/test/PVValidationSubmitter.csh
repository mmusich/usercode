#!/bin/tcsh

echo  "Job started at "
\date;

if ($#argv < 1) then
    echo "************** Argument Error: at least 1 arg. required **************"
    echo "*   Usage:                                                           *"
    echo "*     ./PVValidationSubmitter.csh <TagFile.dat> <options>            *"
    echo "*                                                                    *"
    echo "*   Options:                                                         *"
    echo "*   --dryRun             do not submit the job to lxbatch            *"
    echo "**********************************************************************"
    exit 1
endif

set inputsource=$1
set options=$2


setenv CMSSW_VER 4_1_2
setenv CMSSW_DIR ${CMSSW_BASE}/src/Alignment/OfflineValidation/test
setenv LXBATCH_DIR `pwd`

cd $CMSSW_DIR
eval `scramv1 runtime -csh`
cd $LXBATCH_DIR

cp ${CMSSW_DIR}/${inputsource} .

set jobname=`tail ${inputsource} | grep jobname | awk '{print $2}'`
set datasetpath=`tail ${inputsource} | grep datasetpath | awk '{print $2}'`
set globaltag=`tail ${inputsource} | grep globaltag | awk '{print $2}'`
set alignobj=`tail ${inputsource} | grep alignobj | awk '{print $2}'`
set taggeom=`tail ${inputsource} | grep taggeom | awk '{print $2}'`
set apeobj=`tail ${inputsource} | grep apeobj | awk '{print $2}'`
set tagape=`tail ${inputsource} | grep tagape | awk '{print $2}'`
set validationtype=`tail ${inputsource} | grep validationtype | awk '{print $2}'`
set tracktype=`tail ${inputsource} | grep tracktype | awk '{print $2}'`
set outfile=`tail ${inputsource} | grep outfile | awk '{print $2}'`

echo ${validationtype}

if(${validationtype} == "PrimaryVertexValidation") then
   echo "1 Vertex Validation"
   cp ${CMSSW_DIR}/PVValidation_TEMPL_cfg.py .
   cat PVValidation_TEMPL_cfg.py | sed "s?DATASETTEMPLATE?${datasetpath}?g" | sed "s?GLOBALTAGTEMPLATE?${globaltag}?g" | sed "s?ALIGNOBJTEMPLATE?${alignobj}?g" | sed "s?GEOMTAGTEMPLATE?${taggeom}?g" | sed "s?APEOBJTEMPLATE?${apeobj}?g" | sed "s?ERRORTAGTEMPLATE?${tagape}?g"  | sed "s?VALIDATIONMODULETEMPLATE?${validationtype}?g" |  sed "s?TRACKTYPETEMPLATE?${tracktype}?g" | sed "s?OUTFILETEMPLATE?${outfile}?g" >! ${jobname}_cfg.py
else
   echo "Multi Vertex Validation"
   cp ${CMSSW_DIR}/MultiPVValidation_TEMPL_cfg.py .
   cat MultiPVValidation_TEMPL_cfg.py | sed "s?DATASETTEMPLATE?${datasetpath}?g" | sed "s?GLOBALTAGTEMPLATE?${globaltag}?g" | sed "s?ALIGNOBJTEMPLATE?${alignobj}?g" | sed "s?GEOMTAGTEMPLATE?${taggeom}?g" | sed "s?APEOBJTEMPLATE?${apeobj}?g" | sed "s?ERRORTAGTEMPLATE?${tagape}?g"  | sed "s?VALIDATIONMODULETEMPLATE?${validationtype}?g" |  sed "s?TRACKTYPETEMPLATE?${tracktype}?g" | sed "s?OUTFILETEMPLATE?${outfile}?g" >! ${jobname}_cfg.py
endif

cp ${CMSSW_DIR}/PVValidation_TEMPL.lsf .
cat PVValidation_TEMPL.lsf | sed  "s?JOBNAMETEMPLATE?${jobname}?g" | sed "s?OUTFILETEMPLATE?${outfile}?g" >! ${jobname}.lsf

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

 if (! -d  ${CMSSW_DIR}/test/PVValResults) then
     mkdir ${CMSSW_DIR}/PVValResults
     cp ${outfile}   ${CMSSW_DIR}/PVValResults
     cp ${jobname}.out ${CMSSW_DIR}/PVValResults 
 else
     cp ${outfile}   ${CMSSW_DIR}/PVValResults
     cp ${jobname}.out ${CMSSW_DIR}/PVValResults 
 endif
endif

echo  "Job ended at "
\date;

exit 0
