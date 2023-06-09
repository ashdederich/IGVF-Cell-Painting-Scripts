#!/usr/bin/env nextflow

params.batch=''
params.nreplicates='5'
params.outlier=false
/*cellprof_out_ch=Channel.fromPath("/path/to/profiles/${params.batch}/",type: 'dir')
loaddata_ch=Channel.fromPath("/path/to/load_data_csv/${params.batch}/",type: 'dir')*/
gct_ch=Channel.fromPath("/path/to/gct/${params.batch}/",type: 'dir')

process pyCyto{
    cache 'lenient'
    publishDir "/path/to/", mode: 'copy'
    label 'profiling'

    input:
    file platemap from "/path/to/metadata/platemaps/${params.batch}/platemap/${params.batch}_platemap.txt"
    file barcode from "/path/to/metadata/platemaps/${params.batch}/barcode_platemap.csv"
    file metadata from "/path/to/metadata/external_metadata/${params.batch}_metadata.txt"
    file config from "/path/to/config_files/${params.batch}_config.yml"
    file config_outlier from "/path/to/config_files/${params.batch}_outlier_config.yml"

    output:
    path("profiles/${params.batch}") into pycyto_output_ch
    path("gct/${params.batch}") into gct_output_ch

    script:
    """
    mkdir -p metadata/external_metadata/ metadata/platemaps/${params.batch}/platemap/ profiles load_data_csv gct
    cp \$(cat ${platemap}) metadata/platemaps/${params.batch}/platemap/
    cp \$(cat ${barcode}) metadata/platemaps/${params.batch}
    cp \$(cat ${config}) .
    cp \$(cat ${config_outlier}) .
    cp \$(cat ${metadata}) metadata/external_metadata/
    rsync -Lr /path/to/load_data_csv/${params.batch}/ load_data_csv/${params.batch}
    rsync -Lr /path/to/profiles/${params.batch}/ profiles/${params.batch}
    if [[ ${params.outlier} == false ]]
    then
        profiling_pipeline.py --config ${params.batch}_config.yml
    else
        profiling_pipeline.py --config ${params.batch}_outlier_config.yml
    fi
    """
}

process get_platelist{
    cache 'lenient'
    publishDir "/path/to/quality_control/violinplots/${params.batch}/", mode: 'copy'
    label 'pycytoandr'

    input:
    path("profiles/${params.batch}") from pycyto_output_ch

    output:
    file("platelist.txt") into platelist_batch_ch

    script:
    """
    for plate in profiles/${params.batch}/*_1; do echo \$plate >> platelist.txt; done
    sed -i 's#profiles/${params.batch}/##g' platelist.txt
    """
}

process violin_plot{
    cache 'lenient'
    publishDir "/path/to/quality_control/violinplots/${params.batch}/", mode: 'copy'
    tag "${plateid}"
    label 'replicating'

    input:
    file("platelist.txt") from platelist_batch_ch

    output:
    tuple val("${params.batch}"), file("Replicate_ViolinPlot_BAW.png"), file("Replicate_ViolinPlot_Data.png"), file("percent_replicating.txt") into qc_images

    script:
    """
    #!/usr/bin/env python3
import sys
import os
import pandas
import shutil
src='/path/to/profiles/${params.batch}'
dest='${params.batch}'
pull=shutil.copytree(src,dest)
srcfile='/path/to/IGVF-Cell-Painting-Scripts/Plotting_and_QC/plotting/Percent_Replicating/utils.py'
destfile='utils.py'
pullfile=shutil.copy(srcfile,destfile)
from utils import calculate_percent_replicating_Target, calculate_percent_matching_Target, calculate_percent_replicating_Target_for_plotting, plot_simple_comparison, plot_two_comparisons
batch_dir='${params.batch}'
with open("platelist.txt","r") as plates:
    platelist=plates.readlines()

platelist=[s.replace('\\n','') for s in platelist]
plate_df=calculate_percent_replicating_Target(batch_dir,platelist,n_replicates=${params.nreplicates})
file=open("percent_replicating.txt","w")
plate_df=repr(plate_df)
file.write("Percent Replicating = " + plate_df + "%")
file.close
calculate_percent_replicating_Target_for_plotting(batch_dir,platelist,batch_name="${params.batch}",n_replicates=${params.nreplicates},plot_data_vals=True)
    """
}