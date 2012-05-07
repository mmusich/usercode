#! /bin/tcsh

foreach inputfile (`ls | grep _ `)
    chmod +x $inputfile
    echo $inputfile
    ./$inputfile
end

rm -fr analyzePAT_MC_*_*_*.root
rm -fr analyzePAT_DATA2011_*_*.root
