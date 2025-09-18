#!/usr/bin/env python
# coding: utf-8

"""
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
"""

import os
import sys
import argparse
import xarray as xr
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime

# Constants
S_TO_Y = 86400 * 365
KG_TO_GT = 1e-9 * 1e-3
M2_TO_MKM2 = 1e-6 * 1e-6
M3_TO_MKM3 = 1e-6 * 1e-9

# Dictionary for plot variables
PLOT_VARIABLES = {
    "SMB_Flux": {
        "mntvar": "SMB_Flux",
        "title": "Surface Mass Balance Flux",
        "ncvar": "acabf",
        "scaling_factor": KG_TO_GT * S_TO_Y,
        "unit": "Gt/y",
        "ncfile": "fluxes"
    },
    "BMB_Flux": {
        "mntvar": "BMB_Flux",
        "title": "Basal Mass Balance Flux",
        "ncvar": "libmassbffl",
        "scaling_factor": KG_TO_GT * S_TO_Y,
        "unit": "Gt/y",
        "ncfile": "fluxes"
    },
    "Ice_Discharge": {
        "mntvar": "Ice_Discharge",
        "title": "Ice Discharge",
        "ncvar": "lifmassbf",
        "scaling_factor": KG_TO_GT * S_TO_Y,
        "unit": "Gt/y",
        "ncfile": "fluxes"
    },
    "Ice_flux_at_Grounding_Line": {
        "mntvar": "Ice_flux_at_Grounding_Line",
        "title": "Ice Flux at Grounding Line",
        "ncvar": "ligroundf",
        "scaling_factor": KG_TO_GT * S_TO_Y,
        "unit": "Gt/y",
        "ncfile": "fluxes"
    },
    "Floating_ice_area": {
        "mntvar": "Floating_ice_area",
        "title": "Floating Ice Area",
        "ncvar": "sftflf",
        "scaling_factor": M2_TO_MKM2,
        "unit": "1e6 km2",
        "ncfile": "states"
    },
    "Volume_Above_Flotation": {
        "mntvar": "Volume_Above_Flotation",
        "title": "Volume Above Flotation",
        "ncvar": "lithkaf",
        "scaling_factor": M3_TO_MKM3,
        "unit": "1e6 km3",
        "ncfile": "states"
    },
    "Volume_rate_of_change": {
        "mntvar": "Volume_rate_of_change",
        "title": "Volume Rate of Change",
        "ncvar": "dlithkdt",
        "scaling_factor": 1e-9 * S_TO_Y,
        "unit": "km3/y",
        "ncfile": "fluxes"
    },
    "Land_Ice_Area": {
        "mntvar": "Land_Ice_Area",
        "title": "Land Ice Area",
        "ncvar": "sftgif",
        "scaling_factor": M2_TO_MKM2,
        "unit": "1e6 km2",
        "ncfile": "states"
    }
}

def parse_arguments():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments including:
            - runid (list): List of run IDs.
            - dir_pattern (str): Custom pattern for the NetCDF directory (use {runid} as a placeholder).
            - file_pattern (str): Custom pattern for the NetCDF files (use {runid} and {ftype} as placeholders).
            - o (str): Output figure name.
            - plt (bool): Flag to indicate whether to plot figures.
            - compute_nc (bool): Flag to indicate whether to compute monitoring files.
    """
    parser = argparse.ArgumentParser(description="Process and visualize Elmer/Ice NetCDF data.")
    parser.add_argument("-runid", required=True, nargs='+', help="List of run IDs to process.")
    parser.add_argument("-dir_pattern", default='./data/{runid}-S/ist/????/', 
                        help="Custom pattern for the NetCDF directory (use {runid} as a placeholder).")
    parser.add_argument("-file_pattern", default='{runid}_antarctica_ismip6_{ftype}_????1231.[0-9]*[0-9].nc', 
                        help="Custom pattern for the NetCDF files (use {runid} and {ftype} as placeholders).")
    parser.add_argument("-o", default='output', help="Base name for output figures.")
    parser.add_argument("-plt", action="store_true", help="Flag to plot figures.")
    parser.add_argument("-compute_nc", action="store_true", help="Flag to compute monitoring files.")
    return parser.parse_args()

def parse_dbfile(runid):
    """
    Parse the style database file to retrieve metadata for a given run ID.

    Args:
        runid (str): The run ID to look up.

    Returns:
        tuple: A tuple containing metadata (name, style, color) for the run ID.

    Raises:
        FileNotFoundError: If the run ID is not found in the database file.
        Exception: If there is an error reading the database file.
    """
    try:
        with open('style_elmer.db') as fid:
            for line in fid:
                att = line.split('|')
                if att[0].strip() == runid:
                    return tuple(map(str.strip, att[:4]))
        raise FileNotFoundError(f"RunID {runid} not found in style_elmer.db")
    except Exception as e:
        sys.exit(f"Error reading style_elmer.db: {e}")

def parse_dbbasin():
    """
    Parse the basin database file to retrieve basin metadata.

    Returns:
        tuple: A tuple containing the basin file path, basin variable name, and a dictionary mapping basin IDs to names.

    Raises:
        Exception: If there is an error reading the basin database file.
    """
    try:
        with open('basin_elmer.db') as fid:
            lines = fid.readlines()
        return lines[0].split(':')[1].strip(), lines[1].split(':')[1].strip(), {l.split('|')[0].strip(): l.split('|')[1].strip() for l in lines[2:]}
    except Exception as e:
        sys.exit(f"Error reading basin_elmer.db: {e}")

def save_netcdf(data_vars, filename):
    """
    Save data variables to a NetCDF file.

    Args:
        data_vars (dict): Dictionary of data variables to save.
        filename (str): Name of the output NetCDF file.
    """
    xr.Dataset(data_vars).to_netcdf(filename, engine="netcdf4")
    print(f"Saved {filename}")

def process_basin_data(ds_fluxes, ds_states, mask_dict, runid):
    """
    Process data for each basin and return fluxes and states.

    Args:
        ds_fluxes (xarray.Dataset): Dataset containing flux variables.
        ds_states (xarray.Dataset): Dataset containing state variables.
        mask_dict (dict): Dictionary of masks for each basin.
        runid (str): The run ID.

    Returns:
        tuple: Two dictionaries containing processed fluxes and states data variables.
    """
    data_vars_fluxes, data_vars_states = {}, {}
    for cbasin, mask in mask_dict.items():
        print(f'Processing data for basin {cbasin}')
        ds_int_flx = ds_fluxes.where(mask.compute(), drop=True)
        ds_int_sts = ds_states.where(mask.compute(), drop=True)
        
        for key, props in PLOT_VARIABLES.items():
            cvar = props["ncvar"]
            ds = ds_int_sts if props["ncfile"] == "states" else ds_int_flx
            da_int = (ds[cvar] * ds['cell_area']).sum(dim='nantarctica_face') * props["scaling_factor"]
            da_int.attrs.update({'units': props["unit"], 'long_name': f'Integrated {key} over basin {cbasin}'})
            (data_vars_states if props["ncfile"] == "states" else data_vars_fluxes)[f'{key}_{cbasin}'] = da_int

    return data_vars_fluxes, data_vars_states

def compute_nc_data(runid_lst, dir_pattern, file_pattern):
    """
    Compute NetCDF data for the given run IDs.

    Args:
        runid_lst (list): List of run IDs to process.
        dir_pattern (str): Directory pattern for locating input files (use {runid} as a placeholder).
        file_pattern (str): File pattern for locating input NetCDF files (use {runid} and {ftype} as placeholders).

    Returns:
        None
    """
    varflx_lst = [v["ncvar"] for v in PLOT_VARIABLES.values() if v["ncfile"] == "fluxes"] + \
                 ['imbie_subbasins', 'time', 'time_centered', 'cell_area']
    varsts_lst = [v["ncvar"] for v in PLOT_VARIABLES.values() if v["ncfile"] == "states"] + \
                 ['imbie_subbasins', 'time', 'time_instant', 'cell_area']
    
    da_basin = None
    mask_dict = {}
    
    for runid in runid_lst:
        print(f'Processing {runid}')

        # Open a states and fluxes dataset
        ds_states = open_dataset(runid, varsts_lst, "states", dir_pattern, file_pattern)
        ds_fluxes = open_dataset(runid, varflx_lst, "fluxes", dir_pattern, file_pattern)
        
        # Compute masks
        if da_basin is None:
            da_basin, mask_dict = create_basin_masks(ds_states)
        
        print(ds_states)

        data_vars_fluxes, data_vars_states = process_basin_data(ds_fluxes, ds_states, mask_dict, runid)
        
        ds_states.close()
        ds_fluxes.close()

        data_vars_fluxes['basins'] = da_basin
        data_vars_states['basins'] = da_basin
        
        save_netcdf(data_vars_fluxes, f'ismip6_fluxes_{runid}_monitoring.nc')
        save_netcdf(data_vars_states, f'ismip6_states_{runid}_monitoring.nc')

def open_dataset(runid, var_list, file_type, dir_pattern, file_pattern):
    """
    Open a dataset (either 'states' or 'fluxes') for a given run ID.

    Args:
        runid (str): The run ID.
        var_list (list): List of variables to preprocess.
        file_type (str): Type of dataset to open ('states' or 'fluxes').
        dir_pattern (str): Directory pattern for locating input files (use {runid} as a placeholder).
        file_pattern (str): File pattern for locating input NetCDF files (use {runid} and {ftype} as placeholders).

    Returns:
        xarray.Dataset: The opened dataset.
    """
    cdir = dir_pattern.format(runid=runid)
    cfiles = file_pattern.format(runid=runid, ftype=file_type)

    print(f'Opening {cdir}/{cfiles}')
    ds = xr.open_mfdataset(
        f'{cdir}/{cfiles}',
        chunks={'time': 1},
        preprocess=lambda ds: ds[var_list],
        engine="netcdf4"
    )
    return ds

def plot_results(runid_lst, output_name):
    """
    Plot results for the given run IDs.

    Args:
        runid_lst (list): List of run IDs to process.
        output_name (str): Base name for output figures.

    Returns:
        None
    """
    basin_file, basin_var, basin_dict = parse_dbbasin()
    basin_df = xr.open_dataset(basin_file)

    # Build the runid data dictionary
    runid_data = build_runid_data_dict(runid_lst)

    # Extract basin list from the first runid's data_states
    first_runid = runid_lst[0]
    basin_list = sorted({f"{int(i):02d}" for i in runid_data[first_runid]["data_states"]['basins'].values}) + ['00']

    # create figures for each basin
    for cbasin in basin_list:
        fig = plt.figure(figsize=(16, 16), dpi=100, facecolor='w', edgecolor='k')
        fig.suptitle(f'Elmer Monitoring (Basin {basin_dict[cbasin]})')

        # plot the time series
        for i, key in enumerate(PLOT_VARIABLES):
            print(f'Plotting {key}')
            ax = fig.add_subplot(3, 3, i + 1)
            for runid in runid_lst:
                plot_basin_data(ax, runid_data[runid], cbasin, PLOT_VARIABLES[key])

        # add plot of the basin map
        plot_basin_map(fig, basin_df[basin_var], cbasin)

        # save figure 
        save_figure(fig, output_name, cbasin)

def create_basin_masks(ds_states):
    """
    Create masks for each basin.

    Args:
        ds_states (xarray.Dataset): Dataset containing state variables.

    Returns:
        tuple: A tuple containing:
            - da_basin (xarray.DataArray): DataArray representing the basin data.
            - mask_dict (dict): Dictionary of masks for each basin.
    """
    da_basin = ds_states['imbie_subbasins'].isel(time=0).drop(['time_instant', 'time'])
    basinlst = [f"{int(i):02d}" for i in set(da_basin.values)] + ['00']
    mask_dict = {cbasin: (da_basin == int(cbasin) if int(cbasin) > 0 else da_basin > 0) for cbasin in basinlst}
    return da_basin, mask_dict

def build_runid_data_dict(runid_lst):
    """
    Build a dictionary with runid as the key and associated data (states, fluxes, styles) as values.

    Args:
        runid_lst (list): List of run IDs.

    Returns:
        dict: A dictionary with runid as the key and a dictionary of data_states, data_fluxes, name, style, and color as values.
    """
    runid_data = {}
    for runid in runid_lst:
        print(f"Loading data for runid: {runid}")
        data_states = xr.open_dataset(f'ismip6_states_{runid}_monitoring.nc', decode_times=True)
        data_fluxes = xr.open_dataset(f'ismip6_fluxes_{runid}_monitoring.nc', decode_times=True)
        name, style, color = parse_dbfile(runid)[1:]  # Extract name, style, and color from the database file
        runid_data[runid] = {
            "data_states": data_states,
            "data_fluxes": data_fluxes,
            "name": name,
            "style": style,
            "color": color
        }
    return runid_data

def plot_basin_data(ax, runid_data, cbasin, prop):
    """
    Plot data for a specific basin and variable.

    Args:
        ax (matplotlib.axes.Axes): The axis to plot on.
        runid_data (dict): Dictionary containing data_states, data_fluxes, name, style, and color for a runid.
        cbasin (str): The basin ID.
        prop (dict): Properties of the variable (e.g., title, unit, ncfile, scaling factor).
    """
    ds = runid_data["data_states"] if prop["ncfile"] == "states" else runid_data["data_fluxes"]
    df = ds[f'{prop["mntvar"]}_{cbasin}'].to_pandas()
    df.index = df.index.strftime('%Y').astype(int)
    if prop["ncfile"] == "states":
        df.index = (df.index - 1)
    df.plot(ax=ax, linestyle=runid_data["style"], color=runid_data["color"], label=runid_data["name"])
    ax.set_title(f'{prop["title"]} [{prop["unit"]}]')  # Use the title from props
    ax.grid(True)

def plot_basin_map(fig, basin_data, cbasin):
    """
    Plot the map for a specific basin.

    Args:
        fig (matplotlib.figure.Figure): The figure object to add the map to.
        basin_data (xarray.DataArray): DataArray representing the basin data.
        cbasin (str): The basin ID to plot.

    Returns:
        None
    """
    proj = cartopy.crs.SouthPolarStereo(central_longitude=0,true_scale_latitude=-71)
    ax_map = fig.add_subplot(3, 3, 9, projection=proj)
    data = basin_data.where(basin_data == int(cbasin), drop=True) if int(cbasin) else basin_data.where(basin_data > 0, drop=True)
    data.plot(x='x', y='y', add_colorbar=False)
    ax_map.set_extent((-180, 180, -90, -65), ccrs.PlateCarree())
    add_map_features(ax_map)

def add_map_features(ax_map):
    """
    Add geographical features to the map.

    Args:
        ax_map (cartopy.mpl.geoaxes.GeoAxesSubplot): The map axis to add features to.

    Returns:
        None
    """
    feature = cartopy.feature.NaturalEarthFeature('physical', 'antarctic_ice_shelves_polys', '50m', facecolor='none')
    ax_map.add_feature(feature, linewidth=0.5, edgecolor='k')
    feature = cartopy.feature.NaturalEarthFeature('physical', 'coastline', '50m', facecolor='none')
    ax_map.add_feature(feature, linewidth=0.5, edgecolor='k')

def save_figure(fig, output_name, cbasin):
    """
    Save the figure to a file.

    Args:
        fig (matplotlib.figure.Figure): The figure object to save.
        output_name (str): The base name for the output file.
        cbasin (str): The basin ID to include in the filename.

    Returns:
        None
    """
    fig.subplots_adjust(left=0.05, right=0.98, bottom=0.08, top=0.93, wspace=0.20, hspace=0.20)
    fig.savefig(f'{output_name}_{cbasin}.png', dpi=150, bbox_inches='tight')
    print(f"Saved figure {output_name}_{cbasin}.png")
    print('')

if __name__ == "__main__":
    args = parse_arguments()
    
    if args.compute_nc:
        compute_nc_data(
            args.runid,
            args.dir_pattern,
            args.file_pattern,
    )    
    elif args.plt:
        plot_results(args.runid, args.o)
