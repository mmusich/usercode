#!/bin/sh

function clean_and_resubmit() {
#    grep -E 'MSG-e | MSG-s' $1
#    grep -E 'MSG-s' $1
    grep -E 'exit' $1
    if [ $? = 0 ]; then
	tmp=${1%.*}
	job=${tmp#*/}
	echo $job
	rm log/${job}.log
	rfrm ${CASTOR_HOME}/ZbbAnalysis/test_2012_02_02_11_09_50_MCforJEC/${job}.root
	bsub < LSF/${job}.lsf
    fi
}

export -f clean_and_resubmit

find log -exec bash -c "clean_and_resubmit {}"  \;
