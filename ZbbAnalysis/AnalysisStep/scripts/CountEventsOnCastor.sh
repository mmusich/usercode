#!/bin/sh

CASTORDIR=$CASTOR_HOME/zbb/ZbbSkimSummer11_PAT42X_Sep11_V1/

for DIR in $(nsls $CASTORDIR); do
    for file in $(nsls $CASTORDIR/$DIR); do 	
#        stager_qry -M $CASTORDIR/$DIR/$file 
	stager_get -M $CASTORDIR/$DIR/$file  &
    done
    echo $DIR
done


for DIR in $(nsls $CASTORDIR); do
    let TotEvents=0
# reset counter    
    for file in $(nsls $CASTORDIR/$DIR); do 	
	set -- `edmEventSize -v rfio:$CASTORDIR/$DIR/$file |grep Events`
	let TotEvents+=$4
    done
    echo $DIR $TotEvents
done
