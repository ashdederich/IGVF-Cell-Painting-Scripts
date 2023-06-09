#!/bin/bash
dir=$1

cd $dir/Images

plate=$(cat platename.txt)

find . -type f -name "*ch3*" > FileName_OrigRNA.txt
find . -type f -name "*ch5*" > FileName_OrigDNA.txt
find . -type f -name "*ch4*" > FileName_OrigER.txt
find . -type f -name "*ch1*" > FileName_OrigMito.txt
find . -type f -name "*ch2*" > FileName_OrigAGP.txt
find . -type f -name "*ch8*" > FileName_OrigBrightfield.txt
find . -type f -name "*ch6*" > FileName_OrigHighZBF.txt
find . -type f -name "*ch7*" > FileName_OrigLowZBF.txt

sed -i 's#./##g' FileName_OrigDNA.txt
sed -i 's#./##g' FileName_OrigER.txt
sed -i 's#./##g' FileName_OrigMito.txt
sed -i 's#./##g' FileName_OrigAGP.txt
sed -i 's#./##g' FileName_OrigRNA.txt
sed -i 's#./##g' FileName_OrigBrightfield.txt
sed -i 's#./##g' FileName_OrigHighZBF.txt
sed -i 's#./##g' FileName_OrigLowZBF.txt

sort -o FileName_OrigDNA.txt FileName_OrigDNA.txt
sort -o FileName_OrigER.txt FileName_OrigER.txt
sort -o FileName_OrigMito.txt FileName_OrigMito.txt
sort -o FileName_OrigAGP.txt FileName_OrigAGP.txt
sort -o FileName_OrigAGP.txt FileName_OrigRNA.txt
sort -o FileName_OrigAGP.txt FileName_OrigBrightfield.txt
sort -o FileName_OrigAGP.txt FileName_OrigHighZBF.txt
sort -o FileName_OrigAGP.txt FileName_OrigLowZBF.txt

cat FileName_OrigDNA.txt | awk '{print substr($0,5,2)}' > metadata_col.txt
cat FileName_OrigDNA.txt | awk '{print substr($0,2,2)}' > metadata_row.txt

number_to_letter.r metadata_row.txt
paste metadata_rowforcol.txt metadata_col.txt > tempwell.txt
cat tempwell.txt | tr -d '[:blank:]' > metadata_well.txt
sed -i 's/"//g' metadata_well.txt

pwd=$(cat ${plate}_pwd.txt)

count=$(expr $(cat FileName_OrigDNA.txt | wc -l) / 9)
for (( c=1; c<=$count; c++)); do seq 9 >> metadata_site.txt; done

cell_locations="ER AGP Mito DNA RNA Brightfield HighZBF LowZBF"

for location in $(echo $cell_locations); do linecount=$(cat FileName_Orig${location}.txt | wc -l); for i in $(echo $pwd); do seq 1 ${linecount} | xargs -i -- echo $i > ${location}_pwd.txt; done; for j in $(echo ${plate}); do seq 1 ${linecount} | xargs -i -- echo $j > metadata_plate.txt; done; done	

echo FileName_OrigER,PathName_OrigER,FileName_OrigAGP,PathName_OrigAGP,FileName_OrigMito,PathName_OrigMito,FileName_OrigDNA,PathName_OrigDNA,FileName_OrigRNA,PathName_OrigRNA,FileName_OrigBrightfield,PathName_OrigBrightfield,FileName_OrigHighZBF,PathName_OrigHighZBF,FileName_OrigLowZBF,PathName_OrigLowZBF,Metadata_Plate,Metadata_Well,Metadata_Site,Metadata_Col,Metadata_Row > load_data.csv 

paste -d, FileName_OrigER.txt ER_pwd.txt FileName_OrigAGP.txt AGP_pwd.txt FileName_OrigMito.txt Mito_pwd.txt FileName_OrigDNA.txt DNA_pwd.txt FileName_OrigRNA.txt RNA_pwd.txt FileName_OrigBrightfield.txt Brightfield_pwd.txt FileName_OrigHighZBF.txt HighZBF_pwd.txt FileName_OrigLowZBF.txt LowZBF_pwd.txt metadata_plate.txt metadata_well.txt metadata_site.txt metadata_col.txt metadata_row.txt >> load_data.csv

sed -i 's#./Images/##g' load_data.csv 
rm *.txt
cd ../../

#this is the for loop to use 
for plate in BR0011701{0..1}*; do echo ${plate} > platename.txt; sed -i 's/\__.*//' platename.txt; mv platename.txt ${plate}/Images; cp ~/nextflow_docs/bin/number_to_letter.r ${plate}/Images; cd ${plate}/Images; pwd > ${plate}_pwd.txt; cd ../../; ~/./generate_load_data_broad.sh ${plate}; done
