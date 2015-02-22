#!/bin/bash

CMSSW_DIR=${CMSSW_BASE}/src/AuxCode/SLHCSimPhase2/test/GenSimStep
cd $CMSSW_DIR

#for i in {0..159}; do   
for i in {0..9}; do   
    echo "------ submitting job with seed = $i"
    bsub -o tmp.tmp -q 1nd MinBias_TuneZ2star_14TeV_pythia6_cff_py_GEN_SIM_forPhaseII_TEMPL.lsf $i 300
done


