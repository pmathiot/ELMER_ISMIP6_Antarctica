Script for processing and visualizing Elmer/Ice NetCDF data.

This script provides functionality to:
1. Compute monitoring data from NetCDF files for specified run IDs.
2. Generate plots for various variables over specified basins.

The script supports two main modes of operation:
- Compute NetCDF data: Processes input NetCDF files to generate monitoring data for fluxes and states.
- Plot results: Generates time series plots and basin maps for the processed data.

Usage:
    python valelmer_cpl_reshape.py -runid <run_id> [options]

Command-line Arguments:
    -runid (required): List of run IDs to process.
    -dir_pattern: Custom pattern for the NetCDF directory (use {runid} as a placeholder). Default: './data/{runid}-S/ist/????/'.
    -file_pattern: Custom pattern for the NetCDF files (use {runid} and {ftype} as placeholders). Default: '{runid}_antarctica_ismip6_{ftype}_????1231.[0-9]*[0-9].nc'.
    -o: Base name for output figures. Default: 'output'.
    -plt: Flag to indicate whether to plot figures.
    -compute_nc: Flag to indicate whether to compute monitoring files.

Key Features:
- Processes NetCDF data for specified run IDs and computes integrated fluxes and states for each basin.
- Generates time series plots for variables such as Surface Mass Balance Flux, Ice Discharge, Floating Ice Area, etc.
- Creates basin maps with geographical features for visualization.

Dependencies:
- Python libraries: os, sys, argparse, xarray, pandas, matplotlib, cartopy, datetime.

Example:
    To compute NetCDF data:
        python valelmer_cpl_reshape.py -runid run1 run2 -compute_nc

    To plot results:
        python valelmer_cpl_reshape.py -runid run1 run2 -plt -o output_name
