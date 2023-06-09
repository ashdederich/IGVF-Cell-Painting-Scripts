#!/usr/bin/env nextflow

/*NOTE: All plate IDs within a batch must end in a '_1' to be recognized by Nextflow in this script*/
params.userid=''
params.batch=''
image_dir_ch=Channel.fromPath("/path/to/images/${params.batch}/*_1/",type: 'dir').map { [ it.name, it ] }.into {image_dir_loaddata; image_dir_illumcorrection; image_dir_analysis }
params.cellprofanalysis=''

process pwdloaddata{
    cache 'lenient'
    label 'coreutils'
    tag "${plateid}"

    input:
    tuple val(plateid), path("imagedir") from image_dir_loaddata

    output:
    tuple val(plateid), file("${plateid}_pwd.txt"), path("imagedir") into loaddata_ch

    script:
    """
    readlink ${imagedir} > ${plateid}_pwd.txt
    """
}

process createLoadDataCsvs{
    cache 'lenient'
    publishDir "/path/to/load_data_csv/${params.batch}/${plateid}/", mode: 'copy'
    label 'pycytoandr'
    tag "${plateid}"

    input:
    tuple val(plateid), file("${plateid}_pwd.txt"), path("imagedir") from loaddata_ch
    file lettertonum from "/path/to/nextflow_docs/bin/letter_to_number.r"

    output:
    tuple val(plateid), file("${plateid}_load_data.csv"), path("imagedir") into illum_correction_ch

    script:
    """
    echo ${plateid} > platename.txt
    mv platename.txt ${imagedir}
    mv ${plateid}_pwd.txt ${imagedir}
    cp ${lettertonum} ${imagedir}/letter_to_number.r
    generate_load_data.sh $imagedir
    mv ${imagedir}/load_data.csv ${plateid}_load_data.csv
   """
}

process illuminationMeasurement{
    cache 'lenient'
    publishDir "/path/to/illum/${params.batch}/${plateid}/", mode: 'copy'
    tag "${plateid}"
    label 'cellprof'

    input:
    tuple val(plateid), file("${plateid}_load_data.csv"), path("imagedir") from illum_correction_ch
    file illumpipe from '/path/to/pipelines/Cell_Painting_Illum_8x_remade.cppipe'

    output:
    tuple val(plateid), path("${plateid}_out/${plateid}_illumAGP.npy"), path("${plateid}_out/${plateid}_illumDNA.npy"), path("${plateid}_out/${plateid}_illumER.npy"), path("${plateid}_out/${plateid}_illumMito.npy"), file("${plateid}_load_data.csv"), path("imagedir") into illum_npy_files

    script:
    """
    cp \$(cat ${illumpipe}) .
    cellprofiler -c -r -p Cell_Painting_Illum_8x_remade.cppipe --data-file ${plateid}_load_data.csv -i $imagedir -o ${plateid}_out
    mv ${plateid}_out/IGVF/Plate_illumAGP.npy ${plateid}_out/${plateid}_illumAGP.npy
    mv ${plateid}_out/IGVF/Plate_illumDNA.npy ${plateid}_out/${plateid}_illumDNA.npy
    mv ${plateid}_out/IGVF/Plate_illumER.npy ${plateid}_out/${plateid}_illumER.npy
    mv ${plateid}_out/IGVF/Plate_illumMito.npy ${plateid}_out/${plateid}_illumMito.npy
    """
}

process createIllumLoadDataCsvs{
    cache 'lenient'
    publishDir "/path/to/load_data_csv/${params.batch}/${plateid}/", mode: 'copy'
    tag "${plateid}"
    label 'pycytoandr'

    input:
    tuple val(plateid), path("${plateid}_out/${plateid}_illumAGP.npy"), path("${plateid}_out/${plateid}_illumDNA.npy"), path("${plateid}_out/${plateid}_illumER.npy"), path("${plateid}_out/${plateid}_illumMito.npy"), file("${plateid}_load_data.csv"), path("imagedir") from illum_npy_files

    output:
    tuple val(plateid), file("${plateid}_load_data_with_illum.csv"), path("imagedir"), file("${plateid}_load_data.csv") into illum_loaddata_ch

    script:
    """
    echo ${plateid} > platename.txt
    generate_illumCsv.sh ${plateid}_out ${plateid}_load_data.csv
    mv ${plateid}_out/${plateid}_load_data_with_illum.csv .
    """

}

process cpAnalysis{
    cache 'lenient'
    publishDir "/path/to/profiles/${params.batch}/${plateid}/", mode: 'copy'
    tag "${plateid}"
    label 'cellprof'

    input:
    tuple val(plateid), file("${plateid}_load_data_with_illum.csv"), path("imagedir"), file("${plateid}_load_data.csv") from illum_loaddata_ch
    file analysispipe from params.cellprofanalysis

    output:
    tuple val(plateid), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFCells.csv"), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFCytoplasm.csv"), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFNuclei.csv"), path("${plateid}_load_data_with_illum.csv"), file("${plateid}_load_data.csv") into analysis_output_ch

    script:
    """
    cp \$(cat ${analysispipe}) .
    cellprofiler -c -r -p Cell_Painting_Analysis_IGVF_loaddata.cppipe --data-file ${plateid}_load_data_with_illum.csv -i $imagedir -o ${plateid}_out
    mv ${plateid}_out/IGVF_painting_results/IGVFCells.csv ${plateid}_out/IGVF_painting_results/${plateid}_IGVFCells.csv
    mv ${plateid}_out/IGVF_painting_results/IGVFCytoplasm.csv ${plateid}_out/IGVF_painting_results/${plateid}_IGVFCytoplasm.csv
    mv ${plateid}_out/IGVF_painting_results/IGVFNuclei.csv ${plateid}_out/IGVF_painting_results/${plateid}_IGVFNuclei.csv
    """
}

process zipfiles{
    cache 'lenient'
    publishDir "/path/to/profiles/${params.batch}/${plateid}/", mode: 'copy'
    tag "${plateid}"
    label 'tabix'

    input:
    tuple val(plateid), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFCells.csv"), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFCytoplasm.csv"), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFNuclei.csv"), path("${plateid}_load_data_with_illum.csv"), file("${plateid}_load_data.csv") from analysis_output_ch

    output:
    tuple val(plateid), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFCells.csv.gz"), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFCytoplasm.csv.gz"), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFNuclei.csv.gz"), path("${plateid}_load_data.csv.gz") into analysis_gzip_ch

    script:
    """
    bgzip -@ ${task.cpus} ${plateid}_out/IGVF_painting_results/${plateid}_IGVFCells.csv
    bgzip -@ ${task.cpus} ${plateid}_out/IGVF_painting_results/${plateid}_IGVFCytoplasm.csv
    bgzip -@ ${task.cpus} ${plateid}_out/IGVF_painting_results/${plateid}_IGVFNuclei.csv
    bgzip -@ ${task.cpus} ${plateid}_load_data.csv
    """
}

process prepare_for_PyCyto{
    cache 'lenient'
    publishDir "/path/to/profiles/${params.batch}/${plateid}/", mode: 'copy'
    tag "${plateid}"
    label 'pycytoandr'

    input:
    tuple val(plateid), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFCells.csv.gz"), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFCytoplasm.csv.gz"), path("${plateid}_out/IGVF_painting_results/${plateid}_IGVFNuclei.csv.gz"), path("${plateid}_load_data.csv.gz") from analysis_gzip_ch

    output:
    tuple val(plateid), path("${plateid}.csv.gz"), path("${plateid}_load_data.csv.gz") into pycyto_prepare_ch

    script:
    """
    aggregate.r ${plateid}_out/IGVF_painting_results/${plateid}_IGVFCells.csv.gz ${plateid}_IGVF ${plateid}
    aggregate.r ${plateid}_out/IGVF_painting_results/${plateid}_IGVFCytoplasm.csv.gz ${plateid}_IGVF ${plateid}
    aggregate.r ${plateid}_out/IGVF_painting_results/${plateid}_IGVFNuclei.csv.gz ${plateid}_IGVF ${plateid}
    merge_aggregated.r ${plateid}_out/IGVF_painting_results/${plateid}_IGVFCells_aggregated.csv.gz ${plateid}_out/IGVF_painting_results/${plateid}_IGVFCytoplasm_aggregated.csv.gz ${plateid}_out/IGVF_painting_results/${plateid}_IGVFNuclei_aggregated.csv.gz
    """
}

process createdirs{
    cache 'lenient'
    tag "${plateid}"
    label 'pycytoandr'
    publishDir "/path/to/workspace", mode: 'copy'

    input:
    tuple val(plateid), path("${plateid}.csv.gz"), path("${plateid}_load_data.csv.gz") from pycyto_prepare_ch
    file(barcode_platemap) from "/path/to/metadata/platemaps/${params.batch}/barcode_platemap.csv"
    file(metadata) from "/path/to/metadata/external_metadata/${params.batch}_metadata.txt"
    file(platemap) from "/path/to/metadata/platemaps/${params.batch}/platemap/${params.batch}_platemap.txt"

    output:
    tuple val(plateid), path("metadata/platemaps/${params.batch}/barcode_platemap.csv"), path("metadata/external_metadata/${params.batch}_metadata.txt"), path("metadata/platemaps/${params.batch}/platemap/${params.batch}_platemap.txt"), path("profiles/${params.batch}/${plateid}/${plateid}.csv.gz"), path("load_data_csv/${params.batch}/${plateid}/load_data.csv.gz"), path("gct/${params.batch}/${plateid}/${plateid}.csv.gz") into for_pycyto_ch 
    
    script:
    """
    mkdir -p profiles/${params.batch}/${plateid}
    mkdir -p metadata/external_metadata
    mkdir -p metadata/platemaps/${params.batch}/platemap
    mkdir -p load_data_csv/${params.batch}/${plateid}
    mkdir -p gct/${params.batch}/${plateid}

    cp \$(cat ${barcode_platemap}) metadata/platemaps/${params.batch}/
    cp \$(cat ${platemap}) metadata/platemaps/${params.batch}/platemap/
    cp \$(cat ${metadata}) metadata/external_metadata/
    cp ${plateid}.csv.gz profiles/${params.batch}/${plateid}/
    cp ${plateid}.csv.gz gct/${params.batch}/${plateid}/
    cp ${plateid}_load_data.csv.gz load_data_csv/${params.batch}/${plateid}/
    mv load_data_csv/${params.batch}/${plateid}/${plateid}_load_data.csv.gz load_data_csv/${params.batch}/${plateid}/load_data.csv.gz
    """
}