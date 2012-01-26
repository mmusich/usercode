#! /bin/tcsh

cmsenv
setenv CMSSW_DIR ${CMSSW_BASE}/src/Alignment/OfflineValidation/test

foreach inputfile (`ls ./datFiles | grep Input`)
    bsub -o tmp.tmp -q cmscaf1nd PVValidationSubmitter.csh ${CMSSW_DIR}/datFiles/$inputfile
end
