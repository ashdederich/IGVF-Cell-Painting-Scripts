#!/bin/bash
dir=$1

cd $dir

plate=${dir%__*}

find . -type f -name "*ch5*" > FileName_OrigDNA.txt
find . -type f -name "*ch4*" > FileName_OrigER.txt
find . -type f -name "*ch1*" > FileName_OrigMito.txt
find . -type f -name "*ch2*" > FileName_OrigAGP.txt

cd Images/

pwd > pwd.txt
mv pwd.txt ../

cd ../

for location in $(cat ../cell_locations.txt); do linecount=$(cat FileName_Orig${location}.txt | wc -l); for i in $(cat pwd.txt); do seq 1 ${linecount} | xargs -i -- echo $i > ${location}_pwd.txt; done; for j in $(echo ${plate}); do seq 1 ${linecount} | xargs -i -- echo $j > metadata_plate.txt; done; done	

echo FileName_OrigER,PathName_OrigER,FileName_OrigAGP,PathName_OrigAGP,FileName_OrigMito,PathName_OrigMito,FileName_OrigDNA,PathName_OrigDNA,Metadata_Plate,Metadata_Well,Metadata_Site,Metadata_Col,Metadata_Row > load_data.csv 

paste -d, FileName_OrigER.txt ER_pwd.txt FileName_OrigAGP.txt AGP_pwd.txt FileName_OrigMito.txt Mito_pwd.txt FileName_OrigDNA.txt DNA_pwd.txt metadata_plate.txt ../metadata_well.txt ../metadata_site.txt ../metadata_col.txt ../metadata_row.txt >> load_data.csv

sed -i 's#./Images/##g' load_data.csv 
cd .. 
