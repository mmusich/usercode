#!/bin/tcsh

set option=$1

eval `scramv1 runtime -csh`

setenv CMSSW_DIR ${CMSSW_BASE}/src/ZbbAnalysis/AnalysisStep/

if(${option} == "--acc") then
    cp -pr ${CMSSW_DIR}/test/AcceptanceStudies/scripts $TMPDIR/
    cp -pr ${CMSSW_DIR}/scripts/SubmitThemAll.csh $TMPDIR/scripts/
    cp -pr ${CMSSW_DIR}/scripts/checkFilesOnCastor.csh $TMPDIR/scripts/
    cp -pr ${CMSSW_DIR}/macros/ScanTheTreeForCorrectionFactors.C  $TMPDIR/scripts/
else
    cp -pr ${CMSSW_DIR}/test/scripts $TMPDIR/
    cp -pr ${CMSSW_DIR}/scripts/SubmitThemAll.csh $TMPDIR/scripts/
    cp -pr ${CMSSW_DIR}/scripts/checkFilesOnCastor.csh $TMPDIR/scripts/
endif

cd $TMPDIR/scripts
chmod +x $TMPDIR/scripts/*.*
./SubmitThemAll.csh
