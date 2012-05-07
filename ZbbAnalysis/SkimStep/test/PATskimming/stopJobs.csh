#!/bin/tcsh

set status=$1
echo "killing job in ${status} status"

foreach job ( `bjobs | grep ${status} | awk '{split($0,a," "); print a[1]}'`) 
    bkill $job 
end
