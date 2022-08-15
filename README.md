These are some scripts that have made it easier to run different analysis steps for the IGVF Cell Painting Consortium.

These steps are configured for compound-based experiments. The lab I am working for is moving into CRISPRi and CRISPRa for cardiac cells. If any of the assumptions change for generating the appropriate files, these scripts will be updated to reflect this.

I have the scripts broken up into a few different sections:
    - "Preparing_For_CellProfiler" 
        * This includes two scripts: a file to generate the load_data.csv file for illumination correction and a file to generate the load_data_with_illum.csv file for analysis.
    - "Preparing_For_Cytominer"
        * This includes a script to aggregate data into the median value for the well position and a script to merge the files together for cytominer analysis.
    - "Plotting_Results"
        * This includes two scripts to plot results of the CellProfiler data and the Cytominer data.

