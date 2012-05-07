#! /bin/tcsh

@ i = 1
foreach inputfile (`ls | grep _ `)
    while ($i <= 1)
	set string=`head -n 1 $inputfile`
	set folder=`echo $string | awk '{split($0,a,"/acc"); print a[1]}' | awk '{split($0,b," "); print b[2]}'`
	echo $folder   
	./checkFilesOnCastor.csh $folder --get 
	@ i += 1
    end
    chmod +x $inputfile
    echo $inputfile
    ./$inputfile
end

rm -fr acceptanceCalculator_*_PU_*.root
rm -fr analyzePAT_DATA2011_*_*.root
rm -fr analyzePAT_MC_*_*.root
