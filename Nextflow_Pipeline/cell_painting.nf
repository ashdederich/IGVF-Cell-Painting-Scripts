#!/usr/bin/env nextflow

/*this is the directory that leads to all of the plates*/
image_dir_ch=Channel.fromPath('/project/shared/gcrb_igvf/ashley/*')
/* cell profiler pipeline for illumination correction */
illum_correct_pipe=Channel.fromPath('/path/to/cellprofiler/illum/pipeline/*.cppipe')
/* cell profiler pipeline for analysis */
analysis_pipe=Channel.fromPath('/path/to/cellprofiler/analysis/pipeline/*.cppipe')

process createLoadDataCsvs{
    tag "${plateid}"

    input:
    path image_dir_ch
    
    ouput:
    tuple val (plateid), file("{plateid}_load_data.csv") into load_data_csv_ch

    script:
   ./generate_load_data.sh ${image_dir_ch}
}

process illuminationMeasurement{
    input:
    /*load_data.csv file*/
    tuple val (plateid), file("{plateid}_load_data.csv") from load_data_csv_ch
    /*path to the images*/
    path image_dir_ch
    /*file for cellprofiler pipeline*/
    file cellpipe from illum_correct_pipe

    output:

    script:
    cellprofiler -c -r -p Cell_Painting_Illum_8x_forbroad.cppipe --data-file load_data_${plate_id}.csv -i /*need to specify a directory to the images*/

/*I need to download cellprofiler and pycytominer into a singularity or docker container*/
    
}

process createIllumLoadDataCsvs{

}

process cpAnalysis{

}

process pycytoAnalysis{

}