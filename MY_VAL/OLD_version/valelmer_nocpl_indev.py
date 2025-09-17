#!/usr/bin/env python
# coding: utf-8

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
STOY = 86400 * 365
KG_TO_GT = 1e-9 * 1e-3
M2_TO_MKM2 = 1e-6 * 1e-6
M3_TO_MKM3 = 1e-6 * 1e-9

PLOT_KEYS = ['SMB_Flux', 'BMB_Flux', 'Ice_Discharge', 'Ice_flux_at_Grounding_Line',
             'Floating_ice_area', 'Volume_Above_Flotation', 'Volume_rate_of_change', 'Land_Ice_Area']

PLOT_SF    = [KG_TO_GT * STOY] * 4 + [M2_TO_MKM2, M3_TO_MKM3, 1e-9 * STOY, M2_TO_MKM2]
PLOT_UNITS = ['Gt/y'  , 'Gt/y'       , 'Gt/y'     , 'Gt/y'     , '1e6 km2', '1e6 km3', 'km3/y'   , '1e6 km2']
PLOT_NCVAR = ['acabf' , 'libmassbffl', 'lifmassbf', 'ligroundf', 'sftflf' , 'lithkaf', 'dlithkdt', 'sftgif' ]
PLOT_NCFIL = ['fluxes', 'fluxes'     , 'fluxes'   , 'fluxes'   , 'states' , 'states' , 'fluxes'  , 'states' ]


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-runid", required=True, nargs='+', help="Run ID list")
    parser.add_argument("-dir", default='EDDIR', help="Directory of input files")
    parser.add_argument("-o", default='output', help="Output figure name")
    parser.add_argument("-plt", action="store_true", help="Plot figures")
    parser.add_argument("-compute_nc", action="store_true", help="Compute monitoring files")
    return parser.parse_args()


def parse_dbfile(runid):
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
    try:
        with open('basin_elmer.db') as fid:
            lines = fid.readlines()
        return lines[0].split(':')[1].strip(), lines[1].split(':')[1].strip(), {l.split('|')[0].strip(): l.split('|')[1].strip() for l in lines[2:]}
    except Exception as e:
        sys.exit(f"Error reading basin_elmer.db: {e}")

def compute_nc_data(runid_lst, cdir):
    varflx_lst = [PLOT_NCVAR[i] for i, f in enumerate(PLOT_NCFIL) if f == 'fluxes'] + ['imbie_subbasins', 'time', 'time_centered', 'cell_area']
    varsts_lst = [PLOT_NCVAR[i] for i, f in enumerate(PLOT_NCFIL) if f == 'states'] + ['imbie_subbasins', 'time', 'time_instant', 'cell_area']
    
    da_basin = None
    mask_dict = {}
    
    for runid in runid_lst:
        print(f'Processing {runid}')
        
        cfile_states = f'{cdir}/{runid}/{runid}_S/ismip6_states_{runid.lower()}_???.nc'
        cfile_fluxes = f'{cdir}/{runid}/{runid}_S/ismip6_fluxes_{runid.lower()}_???.nc'
        print(f'open {cfile_states}')        
        ds_states = xr.open_mfdataset(cfile_states, concat_dim='time', chunks={'time': 1}, preprocess=lambda ds: ds[varsts_lst], engine="netcdf4")
        print(f'open {cfile_fluxes}')        
        ds_fluxes = xr.open_mfdataset(cfile_fluxes, concat_dim='time', chunks={'time': 1}, preprocess=lambda ds: ds[varflx_lst], engine="netcdf4")
        
        print(f'compute mask')
        if da_basin is None:
            da_basin = ds_states['imbie_subbasins'].isel(time=0).drop(['time_instant', 'time'])
            basinlst = [f"{int(i):02d}" for i in set(da_basin.values)] + ['00']
            
            for cbasin in basinlst:
                ibasin = int(cbasin)
                if ibasin > 0:
                    mask_dict[cbasin] = da_basin == ibasin
                else:
                    mask_dict[cbasin] = da_basin > 0
        
        data_vars_fluxes, data_vars_states = {}, {}
        
        for cbasin, mask in mask_dict.items():
            print(f'process data {cbasin}')
            ds_int_flx = ds_fluxes.where(mask, drop=True)
            ds_int_sts = ds_states.where(mask, drop=True)
            
            for ikey, ckey in enumerate(PLOT_KEYS):
                cvar = PLOT_NCVAR[ikey]
                ds = ds_int_sts if PLOT_NCFIL[ikey] == 'states' else ds_int_flx
                da_int = (ds[cvar] * ds['cell_area']).sum(dim='nmesh2D_face') * PLOT_SF[ikey]
                da_int.attrs.update({'units': PLOT_UNITS[ikey], 'long_name': f'Integrated {ckey} over basin {cbasin}'})
                (data_vars_states if PLOT_NCFIL[ikey] == 'states' else data_vars_fluxes)[f'{ckey}_{cbasin}'] = da_int
        
        data_vars_fluxes['basins'] = da_basin
        data_vars_states['basins'] = da_basin
       
        print('write netcdf')
        print('fluxes ...')
        xr.Dataset(data_vars_fluxes).to_netcdf(f'ismip6_fluxes_{runid}_monitoring.nc', engine="netcdf4")
        print('states ...')
        xr.Dataset(data_vars_states).to_netcdf(f'ismip6_states_{runid}_monitoring.nc', engine="netcdf4")

def plot_results(runid_lst, output_name):

    basin_file, basin_var, basin_dict = parse_dbbasin()
    basin_df = xr.open_dataset(basin_file)

    data_states = [xr.open_dataset(f'ismip6_states_{runid}_monitoring.nc', decode_times=True) for runid in runid_lst]
    data_fluxes = [xr.open_dataset(f'ismip6_fluxes_{runid}_monitoring.nc', decode_times=True) for runid in runid_lst]
    print(data_states[0])
    basin_list = sorted({f"{int(i):02d}" for i in data_states[0]['basins'].values}) + ['00']
    print(basin_list)
    styles = [parse_dbfile(runid)[1:] for runid in runid_lst]
    
    for cbasin in basin_list:
        fig=plt.figure(figsize=(16,16), dpi= 100, facecolor='w', edgecolor='k')
        fig.suptitle(f'Elmer Monitoring (Basin {basin_dict[cbasin]})')
        
        for i, (ckey, unit) in enumerate(zip(PLOT_KEYS, PLOT_UNITS)):
            print(f'plot {ckey}')
            ax = fig.add_subplot(3,3,i+1)
            for ds, (name, style, color) in zip(data_states if PLOT_NCFIL[i] == 'states' else data_fluxes, styles):
                df = ds[f'{ckey}_{cbasin}'].to_pandas()
                df.index = df.index.strftime('%Y')  # Convert datetime index to string format (years only)
                if PLOT_NCFIL[i] == 'states':
                    df.index = (df.index.astype(int) - 1).astype(str)  # Convert to int, subtract 1, then back to string
                lg=df.plot(ax=ax, linestyle=style, color=color, label=name)
                # Use ScalarFormatter and disable scientific notation
                formatter = mticker.ScalarFormatter(useMathText=True)
                formatter.set_scientific(False)  # Force normal numbers
                formatter.set_useOffset(False)  # Remove the +xxxx offset
                ax.yaxis.set_major_formatter(formatter)
            ax.set_title(f'{ckey} [{unit}]')
            ax.grid(True)
 
        
        proj = ccrs.Stereographic(central_latitude=-90)
        ax_map = fig.add_subplot(3, 3, 9, projection=proj)
        if int(cbasin) :
            data=basin_df.where(basin_df[basin_var] == int(cbasin), drop=True)
        else:
            data=basin_df.where(basin_df[basin_var] > 0, drop=True)
        data[basin_var].plot(x='lon', y='lat', transform=ccrs.PlateCarree(),add_colorbar=False)
        ax_map.set_extent((-180, 180, -90, -65), ccrs.PlateCarree())

        feature=cartopy.feature.NaturalEarthFeature('physical', 'antarctic_ice_shelves_polys', '50m', facecolor='none')
        ax_map.add_feature(feature,linewidth=0.5,edgecolor='k')
        feature=cartopy.feature.NaturalEarthFeature('physical', 'coastline'                , '50m', facecolor='none')
        ax_map.add_feature(feature,linewidth=0.5,edgecolor='k')
 
        fig.subplots_adjust(left=0.05,right=0.98, bottom=0.08, top=0.93, wspace=0.20, hspace=0.20)

        lax = plt.axes([0.0, 0.0, 1, 0.05])
        lline, llabel = lg.get_legend_handles_labels()
        leg=plt.legend(lline, llabel, loc='upper left', ncol = 4, frameon=False) #,fontsize=14)
        lax.set_axis_off()
 
        fig.savefig(f'{output_name}_{cbasin}.png', dpi=150, bbox_inches='tight')
    #plt.show()


if __name__ == "__main__":
    args = parse_arguments()
    
    if args.compute_nc:
        compute_nc_data(args.runid, args.dir)
    elif args.plt:
        plot_results(args.runid, args.o)

