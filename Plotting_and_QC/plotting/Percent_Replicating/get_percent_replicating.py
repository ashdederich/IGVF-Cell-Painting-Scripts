#!/usr/bin/env python3

#adpated from https://github.com/carpenterlab/2022_Cimini_NatureProtocols/blob/09589b066833849329f7c64f88fecda621460c6d/notebooks/AssessDrugSelection.ipynb

import sys
import pandas
from utils import calculate_percent_replicating_Target, calculate_percent_matching_Target, calculate_percent_replicating_Target_for_plotting, plot_simple_comparison, plot_two_comparisons

batch_dir=sys.argv[1]
platelist = [ sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5] ]

plate_df=calculate_percent_replicating_Target(batch_dir,platelist,n_replicates=5)
plate_df

calculate_percent_replicating_Target_for_plotting(batch_dir,platelist,batch_name="Broad Results, U20S Plates",n_replicates=8,plot_data_vals=True)