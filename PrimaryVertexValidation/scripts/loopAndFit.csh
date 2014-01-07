#!/bin/tcsh

set COUNT=0
set myDir=$1

foreach inputfile (`cmsLs /store/caf/user/musich/Alignment/PVValidation/$myDir`)
    #echo $inputfile
    set namebase=`echo $inputfile |awk '{split($0,b,"/"); print b[9]}'`
    if ("$namebase" =~ *"$myDir"*) then
	echo "copying: /store/caf/user/musich/Alignment/PVValidation/$myDir/PVValidation_$namebase" 
	cmsStage /store/caf/user/musich/Alignment/PVValidation/$myDir/$namebase .
    @ COUNT+=1
    endif  
end

 echo "copied $COUNT files"  
