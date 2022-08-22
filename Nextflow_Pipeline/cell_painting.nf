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
config_ch=Channel.fromPath('/path/to/pycyto/configFiles/*_config.yml')
params.batchID=''
barcode_platemap_ch=Channel.fromPath('/path/to/platemap/barcode_platemap.csv')
platemap_ch=Channel.fromPath('/same/path/as/above/platemap.txt')
metadata_ch=Channel.fromPath('/path/to/external_metadata/metadata.tsv')


process createLoadDataCsvs{
    tag "${plateid}"

    input:
    path image_dir_ch
    
    ouput:
    tuple val (plateid), file("{plateid}_load_data.csv") into load_data_csv_ch

    script:
    """
   ./generate_load_data.sh ${image_dir_ch}
   mv ${image_dir_ch}/load_data.csv ${plateid}_load_data.csv
   """
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
    """
    cellprofiler -c -r -p ${illumpipe} --data-file ${plateid}_load_data.csv -i $image_files -o ${plateid}_out
    """
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
    """
    ./generate_illumCsv.sh ${plateid}_load_data.csv $cellfile
    """

}

process cpAnalysis{
    tag "${plateid}"
    
    input:
    tuple val(plateid), file("${plateid}_load_data_with_illum.csv") into illum_loaddata_ch
    file analysispipe from analysis_pipe
    path(imagefiles) from images_analysis

    output:
    tuple val(plateid), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFCells.csv"), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFCytoplasm.csv"), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFNuclei.csv") into analysis_output_ch
    tuple val(plateid), file("${plateid}_load_data_with_illum.csv") into load_data_for_PyCyto_ch

    script:
    """
    cellprofiler -c -r -p $analysispipe --data-file ${plateid}_load_data_with_illum.csv -i $imagefiles
    mv ${plateid}_out/IGVF_painting_results/IGVFCells.csv ${plateid}_out/IGVF_painting_results/${plateid}_IGVFCells.csv
    mv ${plateid}_out/IGVF_painting_results/IGVFCytoplasm.csv ${plateid}_out/IGVF_painting_results/${plateid}_IGVFCytoplasm.csv
    mv ${plateid}_out/IGVF_painting_results/IGVFNuclei.csv ${plateid}_out/IGVF_painting_results/${plateid}_IGVFNuclei.csv
    """
}

process prepare_for_PyCyto{
    tag "${plateid}"
    
    input:
    tuple val(plateid), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFCells.csv"), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFCytoplasm.csv"), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFNuclei.csv") from analysis_output_ch

    output:
    tuple val(plateid), file("${plateid}.csv.gz") into pycyto_prepare_ch

    script:
    """
    ./aggregate_cellpainting.r ${plateid}_out/IGVF_painting_results/${plateid}_IGVFCells.csv
    ./aggregate_cellpainting.r ${plateid}_out/IGVF_painting_results/${plateid}_IGVFCytoplasm.csv
    ./aggregate_cellpainting.r ${plateid}_out/IGVF_painting_results/${plateid}_IGVFNuclei.csv
    ./merge_files.r ${plateid}_out/IGVF_painting_results/${plateid}_IGVFCells.csv ${plateid}_out/IGVF_painting_results/${plateid}_IGVFCytoplasm.csv ${plateid}_out/IGVF_painting_results/${plateid}_IGVFNuclei.csv
    bgzip -@ ${task.cpus} ${plateid}.csv
    """
}

process create_PyCyto_Dirs{
    tag "${plateid}"
    tag "${params.batch}"

    input:
    tuple val(plateid), file("${plateid}.csv.gz") from pycyto_prepare_ch
    file(barcode_platemap) from barcode_platemap_ch
    file(platemap) from platemap_ch
    file(metadata) from metadata_ch
    tuple val(plateid), file("${plateid}_load_data_with_illum.csv") from load_data_for_PyCyto_ch

    output:
    /*will need to update this with the actual file names*/
    tuple val(plateid), path("metadata/platemaps/${params.batch}/barcode_platemap.csv"), path("metadata/external_metadata/metadata.tsv"), path("metadata/platemaps/${params.batch}/platemap/platemap.txt"), path("profiles/${params.batch}/${plateid}/${plateid}.csv.gz"), path("load_data_csv/${params.batch}/${plateid}/${plateid}_load_data_with_illum.csv") into for_pycyto_ch 
    
    script:
    """
    mkdir -p profiles/${params.batch}/${plateid}
    mkdir -p metadata/external_metadata
    mkdir -p metadata/platemaps/${params.batch}/platemap
    mkdir -p load_data_csv/${params.batch}/${plateid}

    ln -s $barcode_platemap metadata/platemaps/${params.batch}/
    ln -s $platemap metadata/platemaps/${params.batch}/platemap/
    ln -s $metadata metadata/external_metadata/
    ln -s ${plateid}.csv.gz profiles/${params.batch}/${plateid}/
    ln -s ${plateid}_load_data_with_illum.csv load_data_csv/${params.batch}/${plateid}/
    """
}

process pyCyto{
    tag "${plateid}"
    tag "${params.batch}"

    input:
    tuple val(plateid), file("${plateid}_config.yml") from config_ch
    tuple val(plateid), path("metadata/platemaps/${params.batch}/barcode_platemap.csv"), path("metadata/external_metadata/metadata.tsv"), path("metadata/platemaps/${params.batch}/platemap/platemap.txt"), path("profiles/${params.batch}/${plateid}/${plateid}.csv.gz"), path("load_data_csv/${params.batch}/${plateid}/${plateid}_load_data_with_illum.csv") from for_pycyto_ch

    output:
    tuple val(plateid), path("profiles/${params.batch}/${plateid}/${plateid}.csv.gz"), path("profiles/${params.batch}/${plateid}/${plateid}_normalized_feature_select_batch.csv.gz"), path("profiles/${params.batch}/${plateid}/${plateid}_normalized_feature_select_negcon_batch.csv.gz") into pycyto_output_ch
    tuple val(plateid), path("quality_control/heatmap/${params.batch}/${plateid}/*.png"), path("quality_control/summary/summary.tsv") into quality_control_ch
    tuple val(plateid), path("metadata/external_metadata/metadata.tsv"), path("metadata/platemaps/${params.batch}/platemap/platemap.txt") into metadata_ch


    script:
    """
    python profiling_pipeline.py --config  ${plateid}_config.yml
    """

}

process qc_figs{
    tag "${plateid}"
    tag "${params.batch}"

    input:
    tuple val(plateid), path("profiles/${params.batch}/${plateid}/${plateid}.csv.gz"), path("profiles/${params.batch}/${plateid}/${plateid}_normalized_feature_select_batch.csv.gz"), path("profiles/${params.batch}/${plateid}/${plateid}_normalized_feature_select_negcon_batch.csv.gz") from pycyto_output_ch
    tuple val(plateid), path("metadata/external_metadata/metadata.tsv"), path("metadata/platemaps/${params.batch}/platemap/platemap.txt") from metadata_ch

    output:
    tuple val(plateid), file("${plateid}_CellProfiler_Output_violin_plot.pdf"), file("${plateid}_Compound_Normalization_violin_plot.pdf"), file("${plateid}_Negcon_Normalization_violin_plot.pdf") into qc_images

    script:
    """
    ./create_violin_plots_cellprof.r profiles/${params.batch}/${plateid}/${plateid}.csv.gz metadata/external_metadata/metadata.tsv metadata/platemaps/${params.batch}/platemap/platemap.txt CellProfiler_Output
    ./create_violin_plots_pycyto.r profiles/${params.batch}/${plateid}/${plateid}_normalized_feature_select_batch.csv.gz Compound_Normalization
    ./create_violin_plots_pycyto.r profiles/${params.batch}/${plateid}/${plateid}_normalized_feature_select_negcon_batch.csv.gz Negcon_Normalization
    """
}