#!/bin/sh

sample=$1
blank=" "
stack=""

rm -fr dbs.txt
rm -fr dataset.txt
python $DBSCMD_HOME/dbsCommandLine.py --url http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet --query="find file where dataset=$sample and dataset.status like VALID*" > dbs.txt
ls | grep .root dbs.txt > dataset.txt
rm -fr dbs.txt
for file in `less dataset.txt`; do
    #echo "$file"
    #echo $file | awk '{split($0,a,"/zbbPATSkim"); print a[2]}'
    suffix=`echo $file | awk '{split($0,a,"/zbbPATSkim"); print a[2]}'`
    #echo $suffix
    filename="analyzePatBasics$suffix"
    stack="$stack$filename$blank"
done

#echo $stack
hadd analyzePatBasics_merged.root `echo $stack`