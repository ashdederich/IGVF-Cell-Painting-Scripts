#!/bin/bash
dir=$1

cd $dir

platename=${dir%_JobArray*}

echo Plate_illumAGP.npy > illumAGP.txt
echo Plate_illumDNA.npy > illumDNA.txt
echo Plate_illumER.npy > illumER.txt
echo Plate_illumMito.npy > illumMito.txt

pwd > pwd.txt

linecount=$(( $(ls /project/shared/gcrb_igvf/ashley/2021_Chandrasekaran_submitted/2020_11_04_CPJUMP1/images/${platename}*/Images/*.tiff | wc -l) / (8) ))

for location in $(cat ../cell_locations.txt); do for i in $(cat illum${location}.txt); do seq 1 ${linecount} | xargs -i -- echo $i > name_illum${location}.txt; done; done

for i in $(cat pwd.txt); do seq 1 ${linecount} | xargs -i -- echo $i > ${platename}_pwd.txt; done

echo FileName_IllumER,PathName_IllumER,FileName_IllumAGP,PathName_IllumAGP,FileName_IllumMito,PathName_IllumMito,FileName_IllumDNA,PathName_IllumDNA > temp.csv
	
paste -d, name_illumER.txt ${platename}_pwd.txt name_illumAGP.txt ${platename}_pwd.txt name_illumMito.txt ${platename}_pwd.txt name_illumDNA.txt ${platename}_pwd.txt >> temp.csv

paste -d, /project/shared/gcrb_igvf/ashley/2021_Chandrasekaran_submitted/2020_11_04_CPJUMP1/images/${platename}*/load_data.csv temp.csv > load_data_with_illum.csv

cd .. 
