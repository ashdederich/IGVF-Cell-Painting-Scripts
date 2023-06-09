#!/bin/bash

dir=$1
loaddata=$2

cp ${loaddata} ${dir}

cd ${dir}

echo *AGP* > illumAGP.txt
echo *DNA* > illumDNA.txt
echo *ER* > illumER.txt
echo *Mito* > illumMito.txt
echo *RNA* > illumRNA.txt
echo *Brightfield* > illumBrightfield.txt
echo *LowZBF* > illumLowZBF.txt
echo *HighZBF* > illumHighZBF.txt

pwd > pwd.txt

linecount=$(awk 'END { print NR -1 }' load_data.csv)

cell_locations="ER AGP Mito DNA RNA Brightfield LowZBF HighZBF"

for location in $(echo $cell_locations); do for i in $(cat illum${location}.txt); do seq 1 ${linecount} | xargs -i -- echo $i > filename_illum${location}.txt; done; done

for i in $(cat pwd.txt); do seq 1 ${linecount} | xargs -i -- echo $i > full_pwd.txt; done

echo FileName_IllumER,PathName_IllumER,FileName_IllumAGP,PathName_IllumAGP,FileName_IllumMito,PathName_IllumMito,FileName_IllumDNA,PathName_IllumDNA,FileName_IllumRNA,PathName_IllumRNA,FileName_IllumBrightfield,PathName_IllumBrightfield,FileName_IllumLowZBF,PathName_IllumLowZBF,FileName_IllumHighZBF,PathName_IllumHighZBF > temp.csv

paste -d, filename_illumER.txt full_pwd.txt filename_illumAGP.txt full_pwd.txt filename_illumMito.txt full_pwd.txt filename_illumDNA.txt full_pwd.txt filename_illumRNA.txt full_pwd.txt filename_illumBrightfield.txt full_pwd.txt filename_illumLowZBF.txt full_pwd.txt filename_illumHighZBF.txt full_pwd.txt >> temp.csv

paste -d, load_data.csv temp.csv > load_data_with_illum.csv

rm *.txt temp.csv
cd ../../