#!/bin/bash

loaddata=$1
cellfile=$2
platename=$(echo $loaddata | awk -F'[_.]' '{print $1}')

echo *AGP* > illumAGP.txt
echo *DNA* > illumDNA.txt
echo *ER* > illumER.txt
echo *Mito* > illumMito.txt

pwd > pwd.txt
cp ${loaddata} .

linecount=$(awk 'END { print NR -1 }' load_data.csv)

for location in $(cat ../${cellfile}); do for i in $(cat illum${location}.txt); do seq 1 ${linecount} | xargs -i -- echo $i > filename_illum${location}.txt; done; done

for i in $(cat pwd.txt); do seq 1 ${linecount} | xargs -i -- echo $i > ${platename}_pwd.txt; done

echo FileName_IllumER,PathName_IllumER,FileName_IllumAGP,PathName_IllumAGP,FileName_IllumMito,PathName_IllumMito,FileName_IllumDNA,PathName_IllumDNA > temp.csv

paste -d, filename_illumER.txt ${platename}_pwd.txt filename_illumAGP.txt ${platename}_pwd.txt filename_illumMito.txt ${platename}_pwd.txt filename_illumDNA.txt ${platename}_pwd.txt >> temp.csv

paste -d, ${loaddata} temp.csv > ${platename}_load_data_with_illum.csv