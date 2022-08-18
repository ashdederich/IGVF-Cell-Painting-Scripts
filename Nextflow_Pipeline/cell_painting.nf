#!/usr/bin/env nextflow

/*this is the directory that leads to all of the plates*/
image_dir_ch=Channel.fromPath('/project/shared/gcrb_igvf/ashley/*')
/*image files*/
images_ch=Channel.fromPath('/path/to/images/files/plateid/*.tif')
/* cell profiler pipeline for illumination correction */
images_ch.into {images_illumcorrection; images_analysis }
illum_correct_pipe=Channel.fromPath('/path/to/cellprofiler/illum/pipeline/*.cppipe')
/* cell profiler pipeline for analysis */
analysis_pipe=Channel.fromPath('/path/to/cellprofiler/analysis/pipeline/*.cppipe')
/*a text file with the list of cellular locations, one per line. These must match how they are written in the .npy files*/
cellfile_ch=Channel.fromPath('/path/to/cellfile/cell_locations.txt')

process createLoadDataCsvs{
    tag "${plateid}"

    input:
    path image_dir_ch
    
    ouput:
    tuple val (plateid), file("{plateid}_load_data.csv") into load_data_csv_ch

    script:
   ./generate_load_data.sh ${image_dir_ch}
   mv ${image_dir_ch}/load_data.csv ${plateid}_load_data.csv
}

process illuminationMeasurement{
    tag "${plateid}"
    
    input:
    /*load_data.csv file*/
    tuple val(plateid), file("${plateid}_load_data.csv") from load_data_csv_ch
    /*path to the images*/
    path(image_files) from images_illumcorrection
    /*file for cellprofiler pipeline*/
    file illumpipe from illum_correct_pipe

    output:
    tuple val(plateid), path("${plateid}_out/${plateid}_illumAGP.npy"), path("${plateid}_out/${plateid}_illumDNA.npy"), path("${plateid}_out/${plateid}_illumER.npy"), path("${plateid}_out/${plateid}_illumMito.npy") into illum_npy_files
    tuple val(plateid), file("${plateid}_load_data.csv") into load_data_IllumCSV_generation

    script:
    cellprofiler -c -r -p ${illumpipe} --data-file ${plateid}_load_data.csv -i $image_files -o ${plateid}_out
/*I need to download cellprofiler and pycytominer into a singularity or docker container*/   
}

process createIllumLoadDataCsvs{
    tag "${plateid}"
    
    input:
    tuple val(plateid), file("${plateid}_illumAGP.npy"), file("${plateid}_illumDNA.npy"), file("${plateid}_illumER.npy"), file("${plateid}_illumMito.npy") from illum_npy_files
    tuple val(plateid), file("${plateid}_load_data.csv") from load_data_IllumCSV_generation
    file cellfile from cellfile_ch

    output:
    tuple val(plateid), file("${plateid}_load_data_with_illum.csv") into illum_loaddata_ch


    script:
    ./generate_illumCsv.sh ${plateid}_load_data.csv $cellfile

}

process cpAnalysis{
    tag "${plateid}"
    
    input:
    tuple val(plateid), file("${plateid}_load_data_with_illum.csv") into illum_loaddata_ch
    file analysispipe from analysis_pipe
    path(imagefiles) from images_analysis

    output:
    tuple val(plateid), path("${plateid}_out/*") into analysis_output_ch

    script:
    cellprofiler -c -r -p $analysispipe --data-file ${plateid}_load_data_with_illum.csv -i $imagefiles -o ${plateid}_out
}

process pycytoAnalysis{
    tag "${plateid}"
    

}