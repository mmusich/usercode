#!/bin/csh
set version=$1
scramv1 project CMSSW $version
cd $version/src
eval `scramv1 runtime -csh`
cd $CMSSW_BASE/src
addpkg DataFormats/PatCandidates  V06-04-19-01
addpkg PhysicsTools/PatAlgos      V08-06-42
addpkg RecoJets/Configuration     V02-04-17
cvs co -r V00-04-11 RecoBTag/PerformanceDB
cvs co -A -d ZbbAnalysis UserCode/musich/ZbbAnalysis
cvs update -P ZbbAnalysis
#scramv1 b


