#!/bin/tcsh

set folder=$1
foreach filetocopy (`nsls ${folder}`)
  
    set log=`echo ${filetocopy} | awk '{split($0,a,"_"); print a[4]}'`
    set num=`echo ${log} | awk '{split($0,b,"."); print b[1]}'`
    echo $num
    
    if(${num}>102) then
	echo ${filetocopy}
	rfrm ${folder}/$filetocopy
    endif
end

