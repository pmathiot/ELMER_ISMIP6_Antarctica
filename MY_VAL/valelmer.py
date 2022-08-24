#!/usr/bin/env python
# coding: utf-8
print('load modules')
import re
import os
import glob
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import argparse
import sys
import cartopy.crs as ccrs
import cartopy
import xarray as xr
import pandas as pd

print('load functions')
# load arguments
def load_arguments():
    # deals with argument
    parser = argparse.ArgumentParser()
    parser.add_argument("-runid", metavar='runid list'  , help="used to look information in runid.db", type=str, nargs='+' , required=True )
    parser.add_argument("-basin", metavar='basin number', help="basin number"                       , type=str, nargs='+'  , required=False, default=['00'] )
    parser.add_argument("-dir"  , metavar='directory of input file' , help="directory of input file", type=str, nargs=1    , required=False, default=['EDDIR'])
    parser.add_argument("-o"    , metavar='figure_name', help="output figure name without extension", type=str, nargs=1    , required=False, default=['output'])
    parser.add_argument("-noshow" , help="do not display the figure (only save it)"                                        , required=False, action="store_true")
    return parser.parse_args()

# parse style name
def parse_dbfile(runid):
    try:
        lstyle=False
        with open('style_elmer.db') as fid:
            for cline in fid:
                att=cline.split('|')
                if att[0].strip() == runid:
                    cpltrunid = att[0].strip()
                    cpltname  = att[1].strip()
                    cpltline  = att[2].strip()
                    cpltcolor = att[3].strip()
                    lstyle=True
        if not lstyle:
            print( runid+' not found in style_elmer.db' )
            raise Exception

    except Exception as e:
        print( 'Issue with file : style_elmer.db' )
        print( e )
        sys.exit(42)

    # return value
    return cpltrunid, cpltname, cpltline, cpltcolor

def parse_dbbasin():
    try:
        with open('basin_elmer.db') as fid:
            iline=0
            dict_basin={}
            for cline in fid:
                iline += 1
                if iline == 1:
                    cfile_basin=cline.split(':')[1].strip()
                else :
                    dict_basin[cline.split('|')[0].strip()]=cline.split('|')[1].strip()

    except Exception as e:
        print( 'Issue with file : basin_elmer.db' )
        print( e )
        sys.exit(42)

    return cfile_basin,dict_basin

# inputs
print('load arguments and constants')
args=load_arguments()

cdir=args.dir[0]
if cdir == 'EDDIR':
    cdir=os.environ['EDDIR']+'/'

if args.basin[0] == 'ALL':
    BASINs=["%.2d" % i for i in range(20)]
else :
    BASINs=args.basin[:]

RUNIDs=args.runid[:]

cfile_out=args.o[0]

stoy=86400*365
kgtoGt=1e-9*1e-3
m2toMkm2=1e-6*1e-6
m3toMkm3=1e-6*1e-9

plot_keys =['SMB Flux', 'BMB Flux'   , 'Ice Discharge', 'Ice flux at Grounding Line', 'Floating ice area', 'Volume Above Flotation'  , 'Volume rate of change', 'Land Ice Area']
plot_sf   =[kgtoGt*stoy, kgtoGt*stoy   ,  kgtoGt*stoy    ,  kgtoGt*stoy                 ,  m2toMkm2          ,  m3toMkm3                 ,  1e-9*stoy         , m2toMkm2       ]
plot_units=['Gt/y'    , 'Gt/y'       , 'Gt/y'         , 'Gt/y'                      , '1e6 km2'          , '1e6 km3'                 , 'km3/y'                , '1e6 km2'      ]
plot_ncvar=['acabf'   , 'libmassbffl', 'lifmassbf'    , 'ligroundf'                 , 'sftflf'           , 'lithkaf'                 , 'dlithkdt'             , 'sftgif'        ]
plot_ncfil=['fluxes'  , 'fluxes'     , 'fluxes'       , 'fluxes'                    , 'states'           , 'states'                  , 'fluxes'               , 'states'       ]

# get basin definition
cfile_basin,dict_basin=parse_dbbasin()

# create dictionary for all var
dict_df={}

print('start plotting:')
# load all the cvs file
for cbasin in BASINs:
    ibasin=int(cbasin)
    datadf=[]
    plot_sty=[]
    plot_clr=[]
    line_name=[]

    ds_states=[]
    ds_fluxes=[]
    for runid in RUNIDs:
        CONFIG=runid.split('-')[0]
        cfile_dat=cdir+'/'+CONFIG+'/SIMU/'+runid+'/MY_OUTPUT/ismip6_states_'+runid.lower()+'_[0-9]*.nc'
        ds=xr.open_mfdataset(cfile_dat,concat_dim='time')
        if ibasin > 0:
            ds_states.append(ds.where(ds['basins']==ibasin, drop=True))
        elif ibasin == 0:
            ds_states.append(ds.where(ds['basins']>=0, drop=True))

        cfile_dat=cdir+'/'+CONFIG+'/SIMU/'+runid+'/MY_OUTPUT/ismip6_fluxes_'+runid.lower()+'_[0-9]*.nc'
        ds=xr.open_mfdataset(cfile_dat,concat_dim='time')
        if ibasin > 0:
            ds_fluxes.append(ds.where(ds['basins']==ibasin, drop=True))
        elif ibasin == 0:
            ds_fluxes.append(ds.where(ds['basins']>=0, drop=True))

        _, runid_name, styline, styclr = parse_dbfile(runid)
        line_name.append(runid_name)
        plot_sty.append(styline)
        plot_clr.append(styclr)

    title_ext=''
   
    # transform each run data frame to variable data frame
    for ikey,ckey in enumerate(plot_keys):
        cvar=plot_ncvar[ikey]
        df=[]
        for irunid,runid in enumerate(RUNIDs):
            if plot_ncfil[ikey]=='states':
                ds=ds_states[irunid]
            elif plot_ncfil[ikey]=='fluxes':
                ds=ds_fluxes[irunid]
            da_int=(ds[cvar]*ds['cell_area']).sum(dim='nmesh2D_face')*plot_sf[ikey]
            datadf=da_int.to_pandas()
            datadf.name=line_name[irunid]
            df.append(datadf.sort_index())
        
        dict_df[ckey]=pd.concat(df, axis=1)

        dict_df[ckey].index=dict_df[ckey].index.to_datetimeindex().strftime('%Y')
   
    # plot data
    fig=plt.figure(figsize=(16,14), dpi= 100, facecolor='w', edgecolor='k')
    axes=[None]*9
    fig.suptitle('Elmer monitoring (basin '+dict_basin[cbasin]+')'+title_ext,fontsize=18)
    count=0
    for r in range(3):
        for c in range(3):
            count+=1
            if (r,c) != (2,2):
                axes[count-1] = fig.add_subplot(3,3,count)
                ckey=list(dict_df.keys())[count-1]
                dict_df[ckey]=dict_df[ckey].interpolate(limit_area='inside')

                lg=dict_df[ckey].plot(ax=axes[count-1],legend=False,label=line_name,linewidth=2,fontsize=14,color=plot_clr,style=plot_sty)
                axes[count-1].set_title(ckey+' ['+plot_units[count-1]+']',fontsize=16)
                axes[count-1].set_xlabel('')
                axes[count-1].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useOffset=False)
                axes[count-1].grid(True)
                ymin=dict_df[ckey].quantile([.01]).min().min()
                ymax=dict_df[ckey].quantile([.99]).max().max()
                yrange=ymax-ymin
                ymin=ymin-0.05*yrange
                ymax=ymax+0.05*yrange
                axes[count-1].set_ylim(ymin, ymax)

                if r == 2 or (r,c) == (1,2):
                    axes[count-1].set_xlabel('Time [years]',fontsize=16)
            else:
                proj = ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0) 
                axes[count-1] = fig.add_subplot(3,3,9,projection=proj)
                # read netcdf
                basin_df=xr.open_dataset(cfile_basin)
                basin_df=basin_df.where(basin_df.basins == int(cbasin), drop=True)
                # plot data
                basin_df.basins.plot(x='lon', y='lat', transform=ccrs.PlateCarree(),add_colorbar=False)
                axes[count-1].set_extent((-180, 180, -90, -65), ccrs.PlateCarree())
                feature=cartopy.feature.NaturalEarthFeature('physical', 'antarctic_ice_shelves_polys', '50m', facecolor='none')
                axes[count-1].add_feature(feature,linewidth=0.5,edgecolor='k')
                feature=cartopy.feature.NaturalEarthFeature('physical', 'coastline'                , '50m', facecolor='none')
                axes[count-1].add_feature(feature,linewidth=0.5,edgecolor='k')
    
    fig.subplots_adjust(left=0.01,right=0.98, bottom=0.05, top=0.93, wspace=0.20, hspace=0.20)
    
    lax = plt.axes([0.0, 0.0, 1, 0.00])
    lline, llabel = lg.get_legend_handles_labels()
    leg=plt.legend(lline, llabel, loc='upper left', ncol = 4, frameon=False,fontsize=16)
    lax.set_axis_off()
   
    cfile_out_png=cfile_out+'_'+cbasin+'.png' 
    fig.savefig(cfile_out_png, format='png', dpi=150, bbox_inches="tight")
    
    if args.noshow:
       pass
    else:
       plt.show()
