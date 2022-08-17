#!/usr/bin/env nextflow

/* images are input files, each with a specific plateID */
image_dir_ch=Channel.fromPath('/project/shared/gcrb_igvf/ashley/${plateid}')
/* cell profiler pipeline for illumination correction */
/* cell profiler pipeline for analysis */

/*process breakupfiles{
    break up image files into equal sized directories - maybe don't need this, as I would have max 13,824 images per plate.
} */

process createLoadDataCsvs{
    /*use script I created to create the load data csv*/
    input:
    tuple val(plate_id), path(plate_dir) from input_images_ch
    /*images separated by plateid*/
    
    ouput:
    tuple val(plate_id), file("load_data_${plate_id}.csv") into load_data_ch

    script:
    ./generate_load_data.sh ${plate_dir}
    /*create list of plates*/
}

process illuminationMeasurement{
    input:
    tuple val(plate_id), file("load_data_${plate_id}.csv") from load_data_ch
    file("Cell_Painting_Illum_8x_forbroad.cppipe") from illum_cppipe_ch
    path images from '/some/data/file.txt'

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