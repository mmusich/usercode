#!/bin/csh

set CASTORDIR=$1

# reset counter    
set sum = 0

foreach file (`nsls ${CASTORDIR}`) 	
    set evts=`edmEventSize -v -a $CASTORDIR/$file | grep Events`
    set num=`echo ${evts} | awk '{split($0,b," "); print b[4]}'`
    echo "+" ${num} "events ("$file")" 
    @ sum += $num
    echo "=> partial sum:" ${sum}
end
echo "======> total sum:" $sum

