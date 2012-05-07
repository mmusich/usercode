#!/bin/tcsh

if ($#argv < 2) then
    echo "************** Argument Error: at least 2 arg. required **************"
    echo "*   Usage:                                                           *"
    echo "*     ./checkFilesOnCastor.csh castorfolder <options>                *"
    echo "*                                                                    *"
    echo "*   Options:                                                         *"
    echo "*   --copy                locally copy all files in folder           *"
    echo "*   --qry                 checks which files are staged              *"
    echo "*   --get                 stager_get of all files in folder          *"
    echo "*   --rm                  stager_rm -a -M of all files in folder     *"
    echo "*   --rfrm                rfrm all files in folder                   *"
    echo "*   --ls                  lists only not staged files                *"
    echo "**********************************************************************"
    exit 1
endif

set folder=$1
set option=$2

foreach filetocopy (`nsls ${folder}`)
    set string=`stager_qry -M ${folder}/$filetocopy`
    set isStaged=`echo $string | awk '{split($0,a," "); print a[3]}'`
    if(${option} == "--copy") then
        if(${isStaged} == "STAGED") then
            echo "--- copying $filetocopy"
	    xrdcp root://castorcms/${folder}/$filetocopy
        else echo "--- sorry not yet staged"
	endif
    else if(${option} == "--qry") then
	 echo "$filetocopy : $isStaged"
    else if(${option} == "--get") then
      	stager_get -M ${folder}/$filetocopy
    else if(${option} == "--rm") then
      	stager_rm -a -M ${folder}/$filetocopy
    else if(${option} == "--rfrm") then
        rfrm ${folder}/$filetocopy
    else if(${option} == "--ls") then
	 if(${isStaged} != "STAGED") then
	 echo "$filetocopy : $isStaged"
	 endif
    else 
	echo "Unrecognized option. please select --copy, --qry, --get or --rm"
    endif
    if($filetocopy == $2) then 
    echo "reached an end"
        break
    endif
end

